# G2-B — integración productiva del operador dual

## Estado tras G2-B.2-LIVE

G2-B.2 permanece abierto. La cinemática estructural ya cruza el socket mediante
`StructuralStateMessage`, pero el handshake se detuvo antes de modificar el
acoplador: el `DofManifest` no describe conectividad/soportes Q1 y por ello no
permite construir el `H_ps` exacto de G2-A. La clasificación es **OPERATOR**,
no estabilidad espectral. Producción continúa exclusivamente en
`legacyNearestNeighbour` y G3 no se inició.

Fecha: 2026-08-25  
Commit de entrada: `9130bd35d24dd577ff4ba6d88d504f968dac95c0`  
Estado: **FAIL — contrato estructural bloqueante; producción no modificada**

## Decisión temprana obligatoria

G2-B exige detenerse si el contrato deal.II no permite aplicar fuerzas nodales
generalizadas sin reinterpretarlas como tracción. Ése es exactamente el estado
actual.

OpenFOAM escribe `dealiiSolid/robin-in.csv` con columnas

```text
x,y,z,tx,ty,tz,vx,vy,vz
```

Los valores `t` son tracciones en `N/m2`. `read_samples()` los almacena en
`Sample::traction`. Durante cada Newton,
`PdmsSolid::assemble_newton()` ejecuta, para cada punto de cuadratura de la
superficie estructural:

```cpp
const Sample &s = nearest(samples, face_values.quadrature_point(q));
cell_rhs(i) += (phi_i*(s.traction + Z*s.velocity - Z*solid_velocity))
               *face_values.JxW(q);
```

Por tanto la única API disponible:

1. interpreta la entrada como Pa;
2. hace una búsqueda nearest-neighbour independiente de H;
3. multiplica por la función de forma y por `JxW` en m²;
4. genera recién entonces una fuerza nodal en N.

El vector dual solicitado

```text
f_s = -H^T W_f t_f
```

ya está integrado y tiene unidades N. No existe parámetro, fichero ni overload
de `assemble_newton()` que acepte ese vector. Pasarlo por `Sample::traction`
haría que deal.II lo tratara como N/m² y lo multiplicara otra vez por área:
sería doble integración y dimensionalmente produciría N·m². Esto activa el
criterio de **FAIL inmediato** del gate.

## Infraestructura que también falta

El protocolo actual tampoco contiene los datos necesarios para que OpenFOAM
aplique el operador C++ productivo verificado en G2-A:

- no exporta IDs globales reales de los 12 948 DoFs deal.II;
- no exporta las filas Q1 de `H_ps` ni su ownership;
- `robin-query.csv` sólo contiene coordenadas, sin IDs ni pesos;
- `robin-in.csv` sólo contiene centros de cara, tracción y velocidad, sin
  conectividad punto-cara, áreas o versión/hash de H;
- la cinemática se evalúa dentro del proceso deal.II, mientras la media
  `P_fp` se ejecuta después en OpenFOAM;
- la fuerza actual vuelve a hacer nearest-neighbour en cuadratura y no puede
  derivarse del almacenamiento cinemático.

Aunque el `standAlonePatch` global está replicado en los ranks OpenFOAM y el
proceso deal.II actual es serial, reunir únicamente tracciones en master no
resuelve el problema: falta un contrato de IDs/filas que haga posible acumular
Hᵀ sin doble conteo y aplicar el resultado al vector estructural correcto.

## Cambio mínimo de API requerido

Antes de reanudar G2-B se necesita un commit estructural explícito con este
contrato mínimo:

1. **Manifiesto cinemático deal.II**

   - Para cada punto FVM consultado: ID global estable, lista de
     `(global_dof_id, component, Q1_weight)` y owner.
   - Dimensiones, nonzeros, versión y hash reproducible.
   - La misma instancia/versión debe producir los valores cinemáticos y recibir
     la aplicación transpuesta.

2. **Topología FVM del operador compuesto**

   - IDs globales de cara y punto, conectividad ordenada, `|S_f|`, owner y
     hash.
   - Con ello se compone `H=P_fp H_ps` una vez; H y Hᵀ recorren el mismo
     almacenamiento.

3. **Entrada de fuerza generalizada**

   - Un modo inequívoco, por ejemplo
     `interface load representation = generalized_dof_force`.
   - Fichero/vector con `(global_dof_id, component, force_N)`.
   - `assemble_newton()` debe sumar este vector directamente al residual
     global después del ensamblaje interno, aplicando las restricciones
     homogéneas de `AffineConstraints`, **sin** `FEFaceValues`, `phi_i`,
     nearest-neighbour ni `JxW`.
   - Validar tamaño, IDs únicos/owners, unidades N, versión/hash y componentes.

4. **Término Robin consistente**

   La ruta actual no sólo aplica tracción física: integra

   ```text
   H^T W_f [ -t_f + Z(v_f-H v_s) ]
   ```

   y su tangente estructural contiene el término coherente
   `Z H^T W_f H * dv_s/dd_s`. El nuevo contrato debe aceptar tanto la fuerza
   generalizada como la contribución matricial derivada del mismo H. Convertir
   únicamente la tracción y dejar el término de impedancia en una cuadratura
   nearest-neighbour mantendría dos operadores y violaría el requisito de
   operador productivo único.

5. **Selección runtime coordinada**

   - `legacyNearestNeighbour`: protocolo actual, default para casos antiguos.
   - `dualConservative`: requiere manifiesto y fuerza/tangente generalizadas;
     debe fallar si faltan o si los hashes/versiones no coinciden.
   - El caso canónico sólo podrá seleccionar `dualConservative` después de que
     ambos participantes implementen el contrato.

Este cambio mínimo es una ampliación de API, no un ajuste de tolerancia ni una
modificación de la formulación física.

## Pruebas que no se ejecutaron

No se activó un modo dual falso y, por tanto, no se ejecutaron como evidencia
dual:

- mínimos serial/MPI;
- restart serial/MPI;
- monitor productivo de trabajo;
- baseline fuerte.

Tampoco fue necesario repetir legacy: ningún archivo productivo, caso,
tolerancia, parámetro o binario fue modificado. G0, G1 y G2-A conservan su
estado anterior.

## Criterios

| criterio G2-B | estado | evidencia |
|---|---|---|
| mismo H para cinemática y dinámica | FAIL | dinámica usa nearest en cuadratura |
| fuerza exacta `-H^T W_f t_f` | FAIL | API sólo acepta tracción |
| ausencia de doble integración | BLOCKED | enviar N por API actual multiplicaría por `JxW` |
| signo y unidades | PASS en especificación | no se activó una ruta incorrecta |
| tests C++ productivos | NOT RUN | operador productivo no puede conectarse todavía |
| serial/MPI y restart dual | NOT RUN | bloqueados por contrato |
| legacy sin regresión | PRESERVED | producción sin cambios |
| árbol limpio/publicado | PASS al cierre documental | sólo documentación |

Estado final: **G2-B FAIL por contrato estructural insuficiente**. La siguiente
acción recomendada es **G2-B.1: ampliar y probar la API de cargas generalizadas
y el manifiesto H**, no G3. No se inició G3.

## Resultado de G2-B.1

El preflight de G2-B.1 encontró un bloqueo previo adicional: el sparsity
pattern generado sólo desde celdas FE no contiene 13 124 pares escalares
requeridos por `H^T W_f Z_f H`. Conforme al criterio explícito, G2-B.1 se
detuvo sin densificar ni implementar una API parcial. La ampliación dispersa
mínima y el manifiesto se especifican en
`docs/G2-B.1-generalized-structural-interface.md`.

G2-B.1a construyó y validó el superset disperso en un participante simulado,
pero no lo confundió con integración productiva: falta conectarlo al
`DynamicSparsityPattern` C++ antes de inicializar la matriz. Producción sigue
en modo legacy.

## Estado posterior a G2-B.1-API

El bloqueo estructural queda resuelto: el patrón C++ se amplía antes de
`reinit()` y `GeneralizedInterfaceLoad.H` acepta fuerza N y la contribución de
Jacobiano `J_Gamma=-df_Gamma/dd` en N/m, con constraints, sellos y checksums.
Las pruebas end-to-end deal.II pasan en serial/MPI-2/MPI-4. Aún no existe
transporte desde el fluido real ni se activa dual en casos productivos; esas son
las tareas acotadas de G2-B.2. Producción continúa legacy sin regresiones.

## Preflight G2-B.2

G2-B.2 quedó **BLOCKED antes de modificar el acoplador**: la auditoría mostró
que `GeneralizedLoad` sólo se invoca desde el arnés y aún no desde
`PdmsSolid::assemble_newton()`. El ejecutable aborta deliberadamente el modo
dual. Conforme al criterio obligatorio no se simuló una conexión viva ni se
usaron archivos temporales como atajo. Se requiere primero G2-B.1-RUNTIME;
véase `docs/G2-B.2-live-dual-transfer.md`.

## Resultado de G2-B.2-RETRY

G2-B.1-RUNTIME corrigió la integración dentro de
`PdmsSolid::assemble_newton()`, pero no expuso esa API a través del ejecutable
que invoca el acoplador como proceso separado. El retry confirmó que el arnés
runtime sigue pasando y que `PdmsSolid::run()` aún rechaza incondicionalmente
dual, sin deserializador de fuerza/tangente generalizadas. El preflight del
binario es por ello **FAIL**. En cumplimiento del orden obligatorio no se tocó
el acoplador, no se activó dual y no se reintentaron casos FSI. El cambio
mínimo pendiente es G2-B.1-EXEC-PROTOCOL; la evidencia está en
`docs/G2-B.2-live-dual-transfer.md`.

## Estado posterior a G2-B.1-EXEC-PROTOCOL

El cambio mínimo ya está implementado y aprobado: `dealiiPdmsSolid` expone la
API runtime mediante un socket Unix enmarcado/versionado, exporta el manifiesto
DoF, recibe el operador antes de `SparseMatrix::reinit()` y activa fuerza y
tangente conjuntamente. El participante externo pasa serial/MPI-2/MPI-4 y las
pruebas negativas. G2-B.2 puede reabrirse en un gate explícito posterior, pero
no fue reintentado aquí; producción FSI continúa legacy.
