# G2-B.2 — conexión dual viva FVM–FEM

## G2-B.2-LIVE — resultado dirigido (2026-08-26)

Commit de entrada: `78c30e9da778e2daa1bd75794bb56bed618a880f`.

Estado: **FAIL/BLOCKED en la etapa B (handshake), clase B — OPERATOR**.
G2 global permanece abierto. No se inició G3.

El preflight obligatorio pasó sin modificar el acoplador: el protocolo externo
produjo `||DeltaR||=1`, `||DeltaJ||=100`, repetición `R/J=0/0` en serial,
MPI-2 y MPI-4. Pasaron las 26 corrupciones/fallos vectorizados del script; junto
con el aborto sin participante documentado por EXEC-PROTOCOL, son los 27
controles negativos aprobados. La regresión runtime produjo
error FD `1.02721e-11` y error 2 con el signo deliberadamente invertido en los
tres tamaños MPI.

### Frontera bidireccional añadida antes del acoplador

El protocolo carecía de cinemática estructural. Se añadieron mensajes
`RequestStructuralState` y `StructuralStateMessage` al mismo framing/socket,
con estado `accepted|provisional`, esquema `BE|BDF2`, unidades explícitas y
entradas canónicas `(global_dof_id, componente, d, v, a)`. También se añadieron
`AcceptTimeStep`, `RejectTimeStep` y `RollbackToAcceptedState`. La aceptación
desplaza las historias una sola vez y el rollback restaura el estado aceptado;
la aceleración aceptada se persiste en el estado extendido (los estados antiguos
siguen siendo legibles). Un proceso externo real recibió y validó el mensaje en
serial, MPI-2 y MPI-4. Una fuerza nula continúa siendo válida por flag, no por
norma.

### Primer fallo de la secuencia

La etapa B no puede construir el operador productivo exigido. El `DofManifest`
actual contiene IDs, componentes, coordenadas de soporte, ownership lógico y
restricciones, pero **no contiene la conectividad de las celdas Q1 ni el soporte
de las funciones de forma**. G2-A obtiene `H_ps` localizando cada punto FVM en
un hexaedro estructural y evaluando sus ocho funciones Q1. Ese resultado no se
puede deducir de una nube de DoFs: dos triangulaciones con los mismos puntos
pueden tener conectividad distinta y producir operadores distintos.

Por tanto, `robinRobinCoupling` no puede validar ni reproducir entrada por
entrada el H de G2-A usando el handshake aprobado. Usar nearest-neighbour,
releer `solid.msh` por una ruta personal o reutilizar los CSV legacy violaría
los criterios del gate. Conforme al orden A--Q, se detuvo antes de modificar
`src/robinRobinCoupling`; no se ejecutaron C--Q, casos duales ni baseline fuerte.

### Cambio mínimo para reabrir LIVE

Versionar una ampliación del manifiesto, antes de conectar OpenFOAM, eligiendo
una representación inequívoca:

1. incluir en `DofManifest` la conectividad global de cada celda Q1, el orden
   deal.II de vértices y los tres DoFs/componentes por nodo; o
2. añadir una fase `FluidPointManifest` y hacer que deal.II devuelva las filas
   canónicas de `H_ps` evaluadas con su `DoFHandler`.

La segunda opción reduce duplicación geométrica y garantiza que el operador
proviene del FE real. Debe conservar checksums, hashes, sesión y validación MPI,
y alimentar el mismo `OperatorManifest` usado para ampliar el patrón. Sólo tras
esa prueba puede comenzar la modificación de `robinRobinCoupling`.

Clasificación: `PROTOCOL=PASS` para la frontera ejecutable existente,
`OPERATOR=FAIL` por manifiesto insuficiente; `ALGEBRA`, `MPI` live, `LIFECYCLE`
live y `SPECTRAL` no fueron alcanzados.

Fecha: 2026-08-26

Commit de entrada: `f88684c877ee91a88f14356a1c7ee939c6652a66`

Estado histórico: **BLOCKED en preflight; bloqueo corregido posteriormente por
G2-B.1-RUNTIME. G2-B.2 continúa pendiente y no fue reintentado.**

## Resultado del preflight obligatorio

La API de G2-B.1 existe como clase C++ reutilizable en
`src/dealiiPdmsSolid/GeneralizedInterfaceLoad.H`, pero todavía no forma parte de
la ruta productiva del ejecutable estructural:

```text
rg GeneralizedLoad src/dealiiPdmsSolid
  tests/test_interface_sparsity.cc: ... GeneralizedLoad ...
  GeneralizedInterfaceLoad.H: class GeneralizedLoad
```

No hay ninguna inclusión ni llamada a `GeneralizedLoad::receive()` o
`add_to_newton()` desde `dealiiPdmsSolid.cc`. En particular:

- `PdmsSolid::assemble_newton()` siempre ensambla la carga Robin legacy con
  `FEFaceValues`, `nearest()`, `JxW` y muestras en Pa;
- el vector `rhs` y la matriz productivos nunca reciben `f_Gamma` en N ni
  `J_Gamma` en N/m;
- la matriz ampliada puede prepararse en `setup()`, pero no existe un canal
  productivo que entregue los mensajes sellados;
- `PdmsSolid::run()` aborta explícitamente si `interface transfer` no es
  `legacyNearestNeighbour`.

Por tanto, los siguientes elementos aprobados en G2-B.1 están demostrados sólo
en el arnés deal.II, no en la API/ruta estructural viva:

| requisito de preflight | estado productivo |
|---|---|
| fuerza generalizada directa en N | **ausente** |
| ruta sin `FEFaceValues/JxW` | **ausente** |
| tangente dispersa en Newton productivo | **ausente** |
| protocolo `timeIndex/outerCorrector/hashes` | clase probada, **sin transporte productivo** |
| matriz replicada por rank | arnés probado, **sin recepción productiva** |
| diferencia finita | arnés probado, **no sobre la llamada productiva** |

Esto satisface exactamente la condición de detención indicada en el apartado
0 del gate: «si alguno existe solamente en el arnés y no en la API C++
estructural real, detenerse y declarar G2-B.2 bloqueado».

## Consecuencias

No se modificó el acoplador, no se intentó el handshake, no se activó
`dualConservative`, no se ejecutaron las etapas A–I y no se creó un caso dual.
Tampoco se ejecutó el baseline fuerte. Hacerlo habría obligado a introducir un
canal temporal o a caer en la cuadratura legacy, ambos prohibidos.

G2-B.2, G2-B y G2 global permanecen **no aprobados**. G2-B.1 requiere una
remediación acotada antes de reabrir este gate.

## Cambio mínimo requerido

Crear un subgate **G2-B.1-RUNTIME** que promueva el contrato ya validado a la
ruta estructural real, sin conectar aún OpenFOAM:

1. añadir a `PdmsSolid` un estado de protocolo propietario del corrector;
2. exponer un transporte versionado/no temporal para recibir el manifiesto,
   `ForceMessage` y `TangentMessage` completos;
3. separar `assemble_newton_legacy()` de la rama generalizada;
4. en la rama generalizada llamar realmente a `receive()` y
   `add_to_newton(rhs,matrix)` después del ensamblaje interno;
5. impedir que esa rama cree `FEFaceValues`, consulte `nearest()` o multiplique
   por `JxW`;
6. implementar el ciclo provisional/aceptado y reconstrucción de restart;
7. verificar sobre el ejecutable productivo, no sólo sobre el arnés, fuerza
   unitaria/nula, tangente, diferencias finitas, mensajes inválidos y réplicas
   MPI;
8. mantener el aborto si falta cualquier mensaje o si los modos/hashes no
   coinciden.

Sólo después G2-B.2 puede implementar el handshake en dos fases y conectar el
fluido real. Este cambio no requiere A_f, masa virtual, IQN-ILS ni ajustes de
física/tolerancias.

## Archivos y evidencia

En este cierre bloqueado sólo se añade documentación. El código y los casos no
cambian. Comandos de auditoría:

```bash
git rev-parse HEAD
git status --short
rg -n 'GeneralizedInterfaceLoad|GeneralizedLoad|add_to_newton' src/dealiiPdmsSolid
rg -n 'dualConservative production transport|assemble_newton' \
  src/dealiiPdmsSolid/dealiiPdmsSolid.cc
```

No se inicia G3.

## Resolución del bloqueo

G2-B.1-RUNTIME conectó `GeneralizedLoad` a `PdmsSolid::assemble_newton()` y
validó el camino real en serial/MPI. Por tanto, el bloqueo descrito aquí ya no
impide reabrir G2-B.2 en un gate posterior. Este trabajo no reintentó handshake
ni modificó el acoplador fluido.

## G2-B.2-RETRY — preflight del ejecutable

Fecha: 2026-08-26

Commit de entrada: `1ff38d6c78a3e24326a4af63fbebfaee5ccdab63`

Estado: **FAIL en preflight; acoplador fluido no modificado**

La regresión versionada de G2-B.1-RUNTIME volvió a pasar en serial, MPI-2 y
MPI-4. En los tres casos produjo `rhs=1`, error de repetición cero, error de
diferencia finita `1,02721e-11` y error 2 al invertir deliberadamente el signo.
Esto confirma la integración de `GeneralizedLoad` en el método real
`PdmsSolid::assemble_newton()`.

Sin embargo, el criterio adicional de este retry exige que **el ejecutable**
estructural acepte dual con datos válidos. Esa condición no se cumple. El
único consumidor de los setters runtime es todavía
`tests/test_runtime_newton.cc`, que incluye la clase en el mismo proceso. El
método productivo `PdmsSolid::run()` ejecuta `setup()` y a continuación exige
incondicionalmente `interface transfer == legacyNearestNeighbour`.

El binario no declara parámetros para mensajes de fuerza y tangente, no
deserializa `ForceMessage` ni `TangentMessage`, y no puede invocar
`set_generalized_interface_data()` desde su interfaz de línea de comandos. No
existe, por tanto, una entrada que permita presentar al ejecutable un par
válido `(f_Gamma,J_Gamma)`; cualquier configuración dual aborta antes de
`solve_newton()`.

La separación de procesos es obligatoria porque OpenFOAM y deal.II cargan
versiones distintas de PETSc. No sería correcto enlazar la clase deal.II en el
plugin OpenFOAM. Tampoco son utilizables los CSV legacy: contienen tracción en
Pa y cinemática por coordenada, carecen de esquema, IDs globales, hashes,
unidades, triplets y sellos de corrector, y reintroducirían
nearest-neighbour/JxW.

Conforme a la instrucción “si falla, detenerse sin modificar el acoplador”, no
se ejecutaron las etapas A--K, no se creó un caso dual, no se ejecutó el
baseline fuerte y no se modificó `src/robinRobinCoupling`.

### Cambio mínimo previo requerido

Reabrir primero **G2-B.1-EXEC-PROTOCOL** para:

1. añadir serialización canónica y versionada de `Expected`, `ForceMessage` y
   `TangentMessage`, con checksums y unidades;
2. permitir al ejecutable recibir el manifiesto antes de
   `SparseMatrix::reinit()` y la fuerza/tangente después de inicializarlo;
3. retirar el aborto dual sólo tras validar y activar ambos mensajes;
4. probar el binario real con datos válidos, cero válido, mensajes incompletos,
   hashes incorrectos y restart sin estado provisional;
5. conservar un transporte neutral entre procesos, sin compartir objetos
   PETSc ni reinterpretar fuerzas N como tracciones Pa.

Después de aprobarlo, G2-B.2 podrá implementar el handshake live en el
acoplador. G2-B.2 permanece **FAIL/BLOCKED** y G3 no se inició.

## Resolución por G2-B.1-EXEC-PROTOCOL

El ejecutable ya recibe `DofManifest`, `OperatorManifest`, fuerza, tangente y
activación atómica desde un proceso externo mediante el protocolo versionado.
Las pruebas serial/MPI y negativas son PASS. El bloqueo de preflight queda
resuelto y la ruta estructural está lista para un futuro retry de G2-B.2. Este
subgate no modificó el acoplador, no reintentó G2-B.2 y no activó casos duales.
