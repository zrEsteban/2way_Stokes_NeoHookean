# G2-B.1 — contrato estructural generalizado

Fecha: 2026-08-25  
Commit de entrada: `22551c1c13689a2de254706bb0dd01dbdfefa49e`  
Estado: **FAIL — el sparsity pattern no admite `H^T W_f Z_f H`**

## Resultado del preflight estructural

El gate ordena detenerse si el patrón estructural no permite `K_Gamma`, sin
densificar como atajo. La matriz actual se construye exclusivamente con:

```cpp
DoFTools::make_sparsity_pattern(dofs, dsp, constraints, false);
```

Ese patrón contiene acoplamientos entre DoFs que comparten una celda Q1. El
operador compuesto de G2-A promedia varios puntos FVM por cara. Esos puntos
pueden evaluarse en celdas estructurales vecinas distintas, por lo que una fila
de H contiene DoFs que no comparten ninguna celda. Su producto exterior en

```text
K_Gamma = H^T W_f Z_f H
```

crea acoplamientos adicionales legítimos.

Se reconstruyó el H exacto de G2-A (3 520 filas, 4 316 columnas nodales,
45 522 nonzeros) y se comparó cada par `(i,j)` de cada fila contra todos los
pares de nodos que comparten uno de los hexaedros Q1 del Gmsh productivo:

```text
filas H analizadas:                    3520
pares escalares ordenados ausentes:   13124
```

Ejemplos de pares de IDs estables ausentes:

```text
(4,5) (4,204) (4,212) (4,213) (4,1231)
```

La auditoría reproducible queda en `scripts/audit-g2b1-sparsity.py`; retorna 1
cuando detecta que el patrón actual es insuficiente.

Para impedancia diagonal por componente, cada par se replica en los tres
bloques físicos. El intento de insertar esos triplets en la `SparseMatrix`
actual violaría su patrón. No se implementó una matriz densa, no se descartaron
entradas y no se proyectó K sobre el stencil FE existente.

## Contrato matemático requerido

Las cantidades permanecen inequívocas:

| cantidad | definición | unidad |
|---|---|---|
| `t_f` | tracción fluida por cara | N/m² |
| `W_f` | áreas/pesos de cara | m² |
| `g_f` | `W_f t_f` | N |
| `d_f` | `H d_s` | m |
| `f_s` | `-H^T g_f` | N |

La identidad es

```text
delta_d_f^T g_f + delta_d_s^T f_s = 0.
```

La condición Robin existente usa impedancia de tracción frente a velocidad:

```text
t_R = -t_f + Z_f (v_f - H v_s),
[Z_f] = (N/m²)/(m/s) = N s/m³.
```

Por ello el operador de velocidad es

```text
K_v = H^T W_f Z_f H,       [K_v] = N s/m,
```

y el Jacobiano respecto de desplazamiento debe ser

```text
K_Gamma = K_v c_v,          [K_Gamma] = N/m,
c_v = dv^{n+1}/dd^{n+1}.
```

Según G1:

```text
c_v(BE)   = 1/dt,
c_v(BDF2) = 3/(2 dt).
```

No es lícito extraer `Z_f` de `H^T W_f Z_f H` salvo que un manifiesto la
declare escalar uniforme y esa igualdad se verifique. El contrato debe admitir
`Z_f` variable por cara.

## Ampliación mínima del sparsity pattern

Antes de implementar la API pedida se requiere:

1. leer y validar el manifiesto H **antes** de crear `sparsity` y `matrix`;
2. ejecutar `DoFTools::make_sparsity_pattern` para conservar el stencil de
   volumen;
3. por cada fila H y por cada componente admitida por `Z_f`, añadir al
   `DynamicSparsityPattern` exactamente el clique de DoFs con pesos no nulos;
4. aplicar `AffineConstraints` a esos índices y cerrar/copiar el patrón;
5. registrar por separado nonzeros FE, nonzeros añadidos y hash del patrón;
6. rechazar cualquier triplet que no pertenezca al patrón ampliado o cuya
   versión/hash no coincida.

Esto es una ampliación dispersa proporcional a los nonzeros reales de
`H^T W_f Z_f H`, no densificación global. Debe probarse primero en una clase
neutral y después integrarse en `PdmsSolid::setup()`.

## Manifiesto neutral propuesto

La serialización canónica, sin PETSc compartido, debe contener:

```text
schema_version
mode = dualConservative
time_index, outer_corrector_index, operator_version, z_version
n_fluid_faces, n_structural_dofs, nnz_H
component_order = xyz_interleaved_by_global_dof
face: global_id owner area_m2
dof: global_id component owner interface=true
H: face_id dof_id weight
units: traction=N/m2 area=m2 force=N displacement=m impedance=N*s/m3
sign: fluid_tf=sigma_f*n_f; solid_fs=-transpose(H)*Wf*tf
topology_sha256, geometry_sha256, operator_sha256
```

El texto debe ordenarse numéricamente por tipo/ID, usar precisión de round-trip
y terminar en newline antes de calcular SHA-256. Ambos participantes deben
rechazar versión, modo, tamaños, hashes, unidades o signo incompatibles. No se
permite fallback a legacy cuando se solicita dual.

## API que queda bloqueada

Tras ampliar el patrón, la ruta de fuerza debe recibir explícitamente pares
`(global_dof_id, component, force_N)`, validar interfaz/owner/versión y sumar
duplicados sólo con regla `ADD` declarada. Debe añadirse directamente al vector
global del residual, con restricciones, sin `FEFaceValues`, nearest-neighbour,
`JxW` ni área.

La estrategia elegida para la tangente será **triplets globales dispersos**
`(row_dof_id,column_dof_id,value_N_per_m)`. Es la opción más verificable para
G2-B.1 porque:

- deal.II puede validar cada entrada contra el patrón exacto ampliado;
- no comparte `PETSc Mat` entre PETSc 3.25.3 y 3.15.5;
- permite comparar directamente con la referencia G2-A;
- fuerza y tangente pueden portar el mismo hash/versiones de H y Z.

El residual aprobado en G1 es `R=R_int+R_inercial-f_ext`. El código actual
ensambla el lado derecho de Newton como `-R`, por lo que una fuerza externa
generalizada se suma a `rhs`; el triplet `K_Gamma` se suma a la matriz del
Newton. Esa pareja deberá verificarse por diferencias finitas después de que el
patrón pueda alojarla.

## Pruebas no ejecutadas

No se implementó una API que no pueda ensamblar su propia tangente. Por tanto
quedan bloqueadas, no fingidas como PASS:

- fuerza unitaria/distribuida y validaciones negativas del manifiesto;
- Z uniforme/variable, simetría y semidefinición positiva;
- acción `K*x` y diferencia finita;
- serial/MPI y restart del nuevo contrato;
- participante simulado completo.

Legacy no fue modificado, de modo que no se recompiló ni repitió G0/G1: sus
fuentes, casos y binarios permanecen iguales al commit de entrada.

## Criterios y siguiente acción

| criterio G2-B.1 | estado |
|---|---|
| API inequívoca de fuerza N | NOT IMPLEMENTED: detenida por patrón |
| ruta consistente `K_Gamma` | FAIL: 13 124 pares ausentes |
| orden `H^T W_f Z_f H` | PASS en derivación |
| fuerza/tangente mismo manifiesto | SPECIFIED |
| diferencias finitas | BLOCKED |
| serial/MPI | BLOCKED |
| legacy sin regresión | PRESERVED; producción intacta |

Estado final: **G2-B.1 FAIL** conforme a la rama obligatoria del prompt. El
cambio mínimo siguiente es **G2-B.1a: ampliar de forma dispersa el patrón desde
el manifiesto H y probar su hash/constraints**, antes de añadir las APIs de
fuerza y triplets. No se inicia G2-B.2 ni G3.

## Resultado de G2-B.1a

El grafo y ensamblaje simulado pasan en serial/MPI: después de Dirichlet quedan
3 168 pares escalares efectivos, el superset tensorial añade 28 512 entradas
vectoriales (9,016 %) y faltan cero conexiones. G2-B.1a permanece FAIL parcial
porque esa ampliación aún no se aplica al `DynamicSparsityPattern` C++ antes de
crear la matriz productiva. Véase
`docs/G2-B.1a-interface-sparsity-extension.md`.
