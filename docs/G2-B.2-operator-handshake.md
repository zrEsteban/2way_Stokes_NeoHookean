# G2-B.2-OPERATOR-HANDSHAKE — evaluación Q1 de H_ps

## Resultado del prerrequisito G2-A-RUNTIME-Q1

El constructor deal.II productivo ya existe, pero quedó `BLOCKED_GEOMETRY`: el
query 20 está separado del boundary 4 por `5.709e-8 m`, equivalente a
`3.481e-5 h_face`. Esto supera el límite `1e-8` y confirma que la extrapolación
del oracle Python era material. El handshake sigue pendiente; no deben añadirse
sus mensajes hasta resolver explícitamente marco o proyección superficial.

La auditoría geométrica dirigida resolvió la ambigüedad: clasificación primaria
`CAD_MISMATCH`. `blockMesh` genera esquinas agudas y la malla sólida usa filetes
circulares de `0.75 microm`. No hay `SEARCH_DEFECT` ni transformación de marco
que lo explique, y una proyección simple no cumple trabajo/momento. La etapa B
permanece bloqueada hasta un subgate de geometrías nominalmente coincidentes;
no debe añadirse todavía `FluidPointManifest`/`HpsRowManifest`.

La geometría maestra filleted ya está especificada y su superficie candidata es
coincidente. G2-A-FILLET-MESH-QUALITY cerró `BLOCKED_MESH_QUALITY` con
`minDet=0.00222376262038637` y `Mesh OK` serial/MPI-2/MPI-4. El handshake sigue
`BLOCKED` únicamente hasta repetir RuntimeQ1 sin reservas; no se reabrió aquí.

Fecha: 2026-08-26

Commit de entrada: `a33dbc305891ab4fd4744be1005532453f730dc2`

Estado: **BLOCKED en preflight; etapa B no aprobada**

G2-B.2 y G2 global permanecen abiertos. No se ejecutaron las etapas C--Q de
G2-B.2-LIVE y no se inició G3.

## Preflight reproducido

| prueba | resultado |
|---|---|
| EXEC-PROTOCOL serial/MPI-2/MPI-4 | PASS; `DeltaR=1`, `DeltaJ=100`, repetición `0/0` |
| controles negativos del protocolo | PASS; 26 casos vectorizados más aborto sin participante |
| G2-B.1-RUNTIME serial/MPI-2/MPI-4 | PASS; FD `1.02721e-11`, signo invertido `2` |
| G2-A serial/MPI-2/MPI-4 | PASS |
| legacy mínimo continuo/5+5 serial y MPI | PASS; End, sin NaN/Inf, `rtol=5e-5` |

Artefactos nuevos de ejecución:

- `/tmp/g2b2-operator-preflight-exec-20260826`;
- `/tmp/g2b2-operator-preflight-runtime-20260826`;
- `/tmp/g2b2-operator-preflight-g2a-20260826`;
- `/tmp/g2b2-operator-preflight-legacy-20260826`.

G2-A reprodujo `H=3520 x 4316`, 45 522 nonzeros, 12 948 DoFs vectoriales,
error de partición de unidad `4.441e-16`, error lineal `1.742e-18 m` y error
rotacional `9.698e-19 m`. Los hashes legacy nuevos fueron:

- serial `fsiResiduals.dat`: `8bc7ad2dd28176344230a6bef559688b6784aab19935a2b415effeea0aa080e0`;
- serial `robin-out.csv`: `cceaabe1995fc8bf45b148a7375f0b632e5afb2c4baf6312c3ad1243c0a14cc1`;
- MPI `fsiResiduals.dat`: `7bc98c11c6423702e2a948c0ac6dea7db2b09ca9665efa3d45885e3322e0a840`;
- MPI `robin-out.csv`: `b6e5fa5648291d58911bca8b6fa0eef08794758f0f486f8b8676f45d03d42f4a`.

## Auditoría de la construcción H_ps aprobada en G2-A

La única implementación que construye H_ps es
`scripts/g2a_dual_transfer.py`. No es una implementación deal.II aislable:

1. `build_operator()` relee `constant/fluid/polyMesh/{points,faces,boundary}`;
2. `gmsh_volume()` relee `dealiiSolid/solid.msh` y extrae hexaedros físicos 10;
3. para cada punto FVM selecciona las 24 celdas volumétricas con centro más
   cercano;
4. `locate_in_cell()` implementa manualmente Newton en coordenadas `(xi,eta,zeta)`
   y funciones trilineales de ocho nodos;
5. elige el menor valor `distancia + extrapolación` y no rechaza explícitamente
   un punto fuera de la interfaz;
6. compone H promediando uniformemente los vértices OpenFOAM de cada cara.

No utiliza `MappingQ1`, `DoFHandler`, `FEFaceValues`, `CellId`, boundary ID
estructural, `AffineConstraints` ni una búsqueda restringida a caras de la
interfaz. Tampoco define una tolerancia relativa a `h_local`, desempate por
`CellId/faceIndex` ni poda explícita de pesos. Sus columnas escalares son IDs de
nodos Gmsh estables, no los `global_dof_id` vectoriales reales de PdmsSolid.

Por tanto existen dos incompatibilidades simultáneas:

- no puede reutilizarse sin releer archivos externos, condición explícita de
  detención del subgate;
- reproducir literalmente su localización volumétrica contradice la obligación
  nueva de localizar sólo sobre el boundary ID estructural y sin extrapolación.

Implementar ahora un localizador nuevo con `MappingQ1` sería un segundo algoritmo
independiente con criterios distintos y no tendría una referencia aprobada con
la cual exigir igualdad de soporte, pesos y hashes. No se improvisó esa ruta.

## Estado del protocolo solicitado

No se añadieron `FluidPointManifest` ni `HpsRowManifest`: un codec sin el
constructor Q1 canónico sólo probaría framing, no el operador requerido. Tampoco
se modificó `robinRobinCoupling`, se avanzó tiempo, se ejecutó un corrector, se
transmitieron cargas reales ni se construyó un patrón nuevo.

La decisión `FluidPointManifest -> evaluación deal.II -> HpsRowManifest` sigue
siendo la arquitectura recomendada frente a exportar conectividad. Antes debe
existir una única implementación canónica evaluable tanto por el arnés como por
el ejecutable.

## Cambio mínimo para desbloquear

Crear un subgate previo y acotado **G2-A-RUNTIME-Q1**:

1. extraer a una clase C++ deal.II la localización sobre caras con el boundary
   ID configurado, `MappingQ1`, `DoFHandler` y `AffineConstraints` reales;
2. fijar marco de referencia, tolerancia respecto de `h_local`, poda y desempate
   determinista de aristas/vértices;
3. sustituir en el arnés G2-A la construcción Python manual por resultados de
   esa clase, sin que OpenFOAM reevalúe Q1;
4. comparar soporte, pesos y acción contra el baseline numérico existente y
   explicar cualquier cambio causado por eliminar extrapolación volumétrica;
5. aprobar una serialización canónica de filas con DoFs vectoriales reales.

Una vez aprobado, OPERATOR-HANDSHAKE puede añadir ambos mensajes, identidad
global OpenFOAM, hashes serial/MPI/restart, composición P_fp/H, OperatorManifest
y Ready. Hasta entonces la etapa B es **FAIL/BLOCKED**, no PASS condicional.

## Archivos modificados

Sólo documentación del repositorio principal. No se modificó código C++, casos,
parámetros físicos ni artefactos de dependencias externas.
