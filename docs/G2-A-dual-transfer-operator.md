# G2-A — operador dual FVM–FEM

## Corrección de alcance tras G2-A-RUNTIME-Q1

`G2-A-OFFLINE=PASS`, `G2-A-RUNTIME-Q1=BLOCKED_GEOMETRY` y
`G2-A-PRODUCTION=INCOMPLETE`. Los resultados Python conservativos permanecen
aprobados como oracle algebraico, pero su búsqueda volumétrica con extrapolación
no representa una evaluación sobre el boundary deal.II real. El constructor C++
productivo detectó la primera separación relativa `3.481e-5 h_face`, superior al
límite `1e-8`; no se añadió proyección ni nearest-neighbour.

`G2-A-GEOMETRY-AUDIT` determinó después la causa primaria `CAD_MISMATCH`. Las
162 consultas incompatibles son exactamente las dos aristas superiores agudas
del patch fluido frente a los filetes `r=7.5e-7 m` del sólido. La discrepancia
máxima fluido→sólido es `4.166776747647578e-7 m`, esencialmente normal, y no es
un defecto de búsqueda, marco, unidades o MappingQ1. G2-A productivo sigue
incompleto hasta corregir la geometría en un subgate separado.

Fecha: 2026-08-25  
Commit de entrada: `4bcef5aca7988ea4aa7380a6a9072500476ddb95`  
Estado: **PASS (especificación y arnés; transferencia productiva intacta)**

## Preflight G1

La protección de `prevAcceleration` usa estado protocolario explícito:
`dealIIFaceAcceleration_[interfaceI].size()` significa que el participante
deal.II ya entregó el campo. `accelerationAvailable` controla la asignación.
No se inspeccionan valores, normas ni tolerancias; un campo recibido y
exactamente nulo es válido. G1 permanece PASS.

## Inventario de espacios y ruta actual

| espacio | tamaño del caso mínimo | representación/IDs |
|---|---:|---|
| DoFs deal.II Q1 vectoriales | 12 948 | 4 316 nodos Gmsh, tres componentes; numeración deal.II local al `DoFHandler` |
| nodos estructurales que soportan la traza | subconjunto de los 4 316 | IDs Gmsh estables; las celdas adyacentes aportan hasta 8 nodos Q1 por punto consultado |
| puntos del patch FVM | 3 645 únicos | índices globales del `standAlonePatch` derivados de `polyMesh/points` |
| caras del patch FVM | 3 520 | caras globales 57 040--60 559 |
| filas escalares de H | 3 520 | una por centro/cara cinemática |
| columnas escalares de H | 4 316 | una por nodo Q1; vectorial = bloque diagonal de tres componentes |

La ruta cinemática productiva es:

1. deal.II almacena `d_s` sobre `FESystem(FE_Q(1),3)`;
2. `VectorTools::point_value` evalúa Q1 en cada punto escrito en
   `robin-query.csv` (`H_ps`);
3. OpenFOAM busca el sample de idéntica coordenada y obtiene velocidad;
4. cada cara promedia uniformemente las velocidades de sus vértices
   (`P_fp`, peso `1/nVertices`);
5. `dealIIFaceVelocity_` es el dato usado por el acoplamiento.

Por tanto:

```text
d_p = H_ps d_s
d_f = P_fp d_p
H   = P_fp H_ps
d_f = H d_s
```

El arnés reconstruye las coordenadas unitarias en la celda hexaédrica Q1 con
el mismo mapa trilineal, usa sus ocho funciones de forma y compone la misma
media punto-cara. No usa una interpolación teórica separada. La distancia
máxima de reconstrucción física fue `6.2071107969972734e-15 m`.

La ruta dinámica productiva actual es diferente y **no fue modificada**:

1. OpenFOAM calcula por cara
   `t_f = faceZoneViscousForce - faceZonePressureForce*n_f`;
2. envía `-t_f` en centros de caras (`robin-in.csv`);
3. deal.II elige de forma independiente el centro FVM nearest-neighbour para
   cada punto de cuadratura estructural;
4. integra la carga con `phi_i*traction*JxW`.

No existe actualmente una aplicación transpuesta del mapa cinemático. Ésta es
precisamente la sustitución reservada para G2-B.

## Operador dual, unidades y signo

Las declaraciones de `fluidModel.H` identifican tanto presión como fuerza
viscosa de patch como `N/m2`. `totalForceOnInterface()` vuelve explícita la
integración al multiplicar por área. deal.II también multiplica la tracción por
`face_values.JxW`. La cantidad transferida es, por tanto, **tracción**, no
fuerza integrada.

| cantidad | unidad |
|---|---|
| `t_f = sigma_f n_f` | Pa = N/m² |
| `W_f = diag(|S_f|)` | m² |
| `W_f t_f` | N |
| H | adimensional |
| `f_s = -H^T W_f t_f` | N |

El signo menos es acción–reacción. `W_f` contiene exactamente el área
poligonal de cada cara global FVM; no se multiplica dos veces por área. La
aplicación H y la aplicación Hᵀ recorren el mismo diccionario de pesos por
fila, de modo que no existe una segunda búsqueda geométrica para fuerzas.

## Trabajo virtual

El arnés verifica

```text
delta_d_f = H delta_d_s
f_s       = -H^T W_f t_f
delta_d_f^T W_f t_f = -delta_d_s^T f_s
```

con vectores deterministas pseudoaleatorios (`seed=20260825`) y con tracción
uniforme. Se normaliza con

```text
e_work = |W_fluido + W_solido| /
         max(|W_fluido|, |W_solido|, 1e-30 J)
```

Las reducciones de trabajo, fuerzas y contribuciones Hᵀ usan
`MPI.COMM_WORLD`. Cada cara tiene ownership único
`globalFaceIndex modulo nRanks`; `Allreduce(SUM)` acumula las contribuciones
al nodo/DoF propietario lógico. El esquema prueba explícitamente ghost
accumulation y ausencia de doble conteo. Los IDs no dependen del rank.

## Resultados cuantitativos

Campaña: `/tmp/g2a-dual-transfer-20260825-3`.

| métrica | serial | MPI 2 (canónico G0) | MPI 4 (segunda partición) |
|---|---:|---:|---:|
| filas escalares / vectoriales | 3 520 / 10 560 | igual | igual |
| columnas escalares / vectoriales | 4 316 / 12 948 | igual | igual |
| nonzeros escalares | 45 522 | igual | igual |
| suma de áreas [m²] | 2.7453599999999997e-6 | igual | igual |
| máximo `|sum(row)-1|` | 4.441e-16 | 4.441e-16 | 4.441e-16 |
| error partición de unidad | 4.441e-16 | 4.441e-16 | 4.441e-16 |
| `e_work` | 2.471e-15 | 1.030e-15 | 2.060e-16 |
| error relativo de resultante | 3.633e-15 | 3.479e-15 | 3.479e-15 |
| error absoluto resultante uniforme [N] | 1.290e-16 | 1.294e-16 | 1.285e-16 |
| error relativo de momento | 6.919e-16 | 1.063e-15 | 1.398e-15 |
| error máximo campo lineal [m] | 1.742e-18 | igual | igual |
| error máximo rotación rígida [m] | 9.698e-19 | igual | igual |

Objetivos de trabajo: serial `<=1e-12`, MPI `<=1e-10`; ambos PASS. H reproduce
traslación, campo lineal y rotación rígida, y su dual conserva resultante y
momento a precisión de máquina.

## Clasificación y diseño G2-B

Clasificación: **A** — H conserva partición de unidad, trabajo, fuerza, momento
y campos lineales. G2-B confirmó posteriormente que su integración requiere
primero ampliar la API estructural: la API actual sólo acepta tracciones y las
integra de nuevo con `JxW`; no acepta fuerzas nodales generalizadas.

Diseño concreto propuesto, no implementado aquí:

1. exportar desde deal.II los IDs globales reales de DoF y las filas Q1 de
   `H_ps` evaluadas en los puntos globales FVM;
2. componer una sola vez con la conectividad punto-cara para obtener H;
3. almacenar cada fila una vez y usar el mismo almacenamiento para `apply()` y
   `apply_transpose()`;
4. asignar ownership por ID global de cara/DoF del runtime (no por orden local),
   acumular ghosts con suma colectiva;
5. reemplazar únicamente en G2-B la búsqueda nearest-neighbour de tracción por
   `-H^T W_f t_f`, manteniendo H como ruta cinemática;
6. repetir este arnés sobre el mapa runtime exportado antes de activar la ruta.

El ownership módulo-rank del arnés es deliberadamente independiente del
particionado geométrico y prueba invariancia algebraica; G2-B aún debe conectar
estos IDs estables con los mapas de ownership reales de OpenFOAM/deal.II.

## Archivos, comandos y riesgos

Archivos nuevos:

- `scripts/g2a_dual_transfer.py`;
- `scripts/test-g2a-dual-transfer.sh`;
- `docs/G2-A-dual-transfer-operator.md`.

Comando:

```bash
G2A_WORK=/tmp/g2a-dual-transfer-20260825-3 \
  scripts/test-g2a-dual-transfer.sh
```

Riesgos abiertos para G2-B:

- la producción sigue usando nearest-neighbour independiente para cargas;
- `assemble_newton()` no dispone de entrada directa para fuerzas nodales N ni
  para la tangente Robin dual;
- el runtime todavía no exporta un manifiesto DoF-ID/owner;
- debe verificarse que el orden global del `standAlonePatch` sea idéntico tras
  descomposiciones OpenFOAM reales, además de la partición adversarial probada;
- la versión de prueba importa NumPy y mpi4py, dependencias diagnósticas que no
  entran al solver.

No se tocó código productivo, Robin, residual FSI, IQN-ILS, integradores ni
solvers. No se ejecutó el baseline fuerte. G2-B y G3 no se iniciaron.
