# G2-B.1a — ampliación dispersa C++ del patrón de interfaz

Fecha: 2026-08-26

Commit de entrada: `37bce9b6157194dfe18d4ef74fe1adb759781add`

Estado: **PASS — G2-B.1a completo; G2-B.1 continúa pendiente**

## Backend e inicialización

La matriz estructural real no es PETSc. La cadena exacta es
`dealii::DynamicSparsityPattern`, `dealii::SparsityPattern::copy_from()` y
`dealii::SparseMatrix<double>::reinit()`.

`PdmsSolid::setup()` ahora distribuye DoFs, construye/cierra
`AffineConstraints`, crea el patrón FEM base y, sólo en
`dualConservative`, exige el manifiesto, mapea nodos a DoFs y ejecuta
`augment_interface_sparsity()`. Después copia el patrón e inicializa la matriz.
No es posible añadir conexiones tardías. En `legacyNearestNeighbour` (default)
no se exige manifiesto ni se ejecuta la ampliación.

## Grafo, Robin y constraints

Cada fila r de H aporta el clique disperso `S_r x S_r`. La Robin vigente usa
un `double solid impedance` y ensambla
`impedance*velocity_coefficient*(phi_i*phi_j)`: `Z_f=z_f I`. Por derivación se
añaden los tres bloques por componente (9.504 entradas), no el superset
tensorial de 28.512.

Cada índice se expande con `AffineConstraints::get_constraint_entries()` para
formar el soporte de `C^T K_Gamma C`; Dirichlet homogéneo sin masters se
elimina. Los pares se ordenan, deduplican y simetrizan. Se validan IDs,
componentes, coordenadas y rangos.

El manifiesto exportado desde el mismo H de G2-A declara:

```text
hashGraph = 6ee94c285c4e5d8619f9f08c3a666eb57dd46ddfe47b2b02be25c1248b4abc39
```

El hash post-constraints/DoF es `30b21a67c782183f` (FNV-1a de pares). Pesos y
áreas no determinan el soporte. Cambiar `hashGraph` después de `reinit()` no es
posible en esta etapa; G2-B.1 añadirá el guard de versión durante Newton antes
de activar la tangente.

## Matriz productiva y resultados

Campaña C++ final: `/tmp/g2b1a-cxx-20260826-final`.

| métrica | resultado |
|---|---:|
| nonzeros productivos legacy | 322.992 |
| conexiones requeridas post-constraints | 14.652 |
| ya existentes | 5.148 |
| añadidas | 9.504 |
| nonzeros finales | 332.496 |
| crecimiento real | 2,942488 % |
| faltantes según `exists(i,j)` | 0 |

La estimación Python previa de 316.224 sólo contaba el grafo volumétrico libre;
322.992 es la medida del DSP real, incluidas entradas conservadas por deal.II
al tratar constraints. Esto explica la diferencia con el crecimiento 9,016 %
del superset tensorial; las 9.504 adiciones coinciden.

`K_Gamma` variable se ensambló directamente en la misma
`dealii::SparseMatrix<double>` final, exigiendo `exists(i,j)` antes de cada
`add`. Se comparó contra `H^T W_f Z_f(Hx)`:

| ejecución | almacenamiento | error Kx | simetría | `x^T Kx` |
|---|---|---:|---:|---:|
| serial | replicado | 4,561e-16 | 0 | 0,10014111691516978 |
| MPI-2 | copia por rank | 4,561e-16 | 0 | idéntico |
| MPI-4 | copia por rank | 4,561e-16 | 0 | idéntico |

La matriz es serial y replicada por rank, no distribuida. Grafo, patrón y Kx
son independientes de la partición fluida. No se usó matriz densa ni PETSc.

## Negativos y regresión legacy

El arnés rechaza manifiesto ausente antes de `reinit()`, hashGraph incompatible
e ID/no orden canónico. El patrón legacy demuestra 9.504 conexiones faltantes;
el extendido exige cero.

Comandos principales:

```bash
G2B1A_CXX_WORK=/tmp/g2b1a-cxx-20260826-final scripts/test-g2b1a-cxx.sh
ctest --test-dir /tmp/g2b1a-cxx-20260826-final/build --output-on-failure
src/dealiiPdmsSolid/tests/run-smoke.sh \
  /tmp/g2b1a-cxx-20260826-final/build/dealiiPdmsSolid
```

En el entorno aislado también pasaron `constantHs` y, usando explícitamente el
binario final, `scripts/test-g1-restart.sh`: continuo/partido serial y MPI, equivalencia
serial/MPI `rtol=5e-5`, sin NaN/Inf. SHA-256 continuos:

- serial: `fsiResiduals.dat` `8bc7ad2d...80e0`, `robin-out.csv` `cceaabe1...cc1`;
- MPI: `fsiResiduals.dat` `7bc98c11...a840`, `robin-out.csv` `b6e5fa56...2f4a`.

BE/BDF2, signo, material y smoke Robin también pasan. Único warning: cabecera
Boost deprecada transitiva de deal.II.

## Archivos y decisión

- `InterfaceSparsityExtension.H`: parser, mapeo, constraints y ampliación;
- `dealiiPdmsSolid.cc`: rama dual de setup, inactiva por default;
- `test_interface_sparsity.cc`: matriz real, Kx, simetría/PSD y MPI;
- `audit-g2b1-sparsity.py`: exportador canónico;
- `test-g2b1a-cxx.sh`: build y pruebas positivas/negativas.

**G2-B.1a-CXX PASS**, completando **G2-B.1a**. Producción continúa legacy. No
se ensamblan fuerzas/tangentes duales dentro de Newton; las APIs, guard completo
de ciclo de vida y diferencias finitas quedan para retomar G2-B.1. No se inició
G2-B.2 ni G3.

## Revalidación por G2-B.1-API

El mismo patrón se usó para ensamblar la nueva tangente generalizada: todas las
entradas pasan `exists(i,j)`, `Jx` coincide exactamente con la forma
factorizada y las réplicas serial/MPI-2/MPI-4 coinciden. G2-B.1 implementó los
sellos y la recepción atómica pendientes; no cambió los conteos 322.992/332.496.
