# G2-A-GEOMETRY-AUDIT — compatibilidad de las interfaces

Fecha: 2026-08-27

Commit de entrada: `71c43b9471c6ca70dbe03137112a2acc270532de`

Estado del subgate: **PASS como auditoría**

Clasificación primaria: **CAD_MISMATCH**

La auditoría no modificó `RuntimeQ1InterfaceOperator`, su tolerancia ni
`src/robinRobinCoupling`; tampoco proyectó puntos, regeneró mallas o ejecutó un
caso FSI dual. Los estados que motivaron el gate no cambian:

- `G2-A-RUNTIME-Q1 = BLOCKED_GEOMETRY`;
- `G2-A-PRODUCTION = INCOMPLETE`;
- `G2-B.2 = BLOCKED`;
- `G2 global = OPEN`.

Los datos completos y reproducibles están en
`tests/g2AInterfaceGeometryAudit/geometry-audit.json` y
`point-distances.csv`. `scripts/audit_g2a_interface_geometry.py` reconstruye la
auditoría desde `polyMesh` y `solid.msh`; el wrapper
`scripts/test-g2a-geometry-audit.sh` exige resultados idénticos en serial,
repetición, MPI-2 y MPI-4.

## Conclusión causal

Las mallas no discretizan una misma superficie nominal. El patch fluido se
genera con `blockMesh` como un canal rectangular con esquinas superiores
agudas. La superficie estructural se genera independientemente con OpenCASCADE
en `roundedSolid.geo` y tiene dos filetes circulares superiores de radio
`r=7.5e-7 m`. Las 162 consultas rechazadas son exactamente los 81 puntos de
cada una de esas dos aristas agudas del fluido. El interior del techo y de las
paredes coincide hasta redondeo.

El error es casi completamente normal a la superficie sólida: distancia normal
máxima `4.166776747647578e-7 m`, frente a distancia tangencial máxima
`1.9947587062916708e-13 m`. Ni una transformación de marco, ni una escala, ni un
offset uniforme explican ese patrón localizado. La búsqueda exhaustiva y la
acelerada producen la misma distancia mínima para los 3 645 puntos. Por ello se
descartan como causa primaria `WRONG_PATCH`, `WRONG_FRAME`, `WRONG_UNITS`,
`QUERY_FIXTURE`, `SEARCH_DEFECT`, `EDGE_EXTENT` y `MAPPING_MISMATCH`.

## Procedencia de las consultas fluidas

La fuente es el `polyMesh` real del caso canónico, patch `interface`, índice
OpenFOAM 3, `startFace=57040`, 3 520 caras y 3 645 puntos únicos. Las unidades
son metros y el marco es la configuración inicial/de referencia. El fixture es
exactamente el conjunto de puntos del patch y conserva el orden de primera
aparición de `standAlonePatch::localPoints()`; no es un CSV redondeado ni un
orden global artificial.

`queryPointId=20` procede del punto global OpenFOAM 20251. Su coordenada
almacenada y leída es exactamente `(6.15e-5, 0, 3.3e-5) m`; participa en las
caras globales 57049, 57059, 58640 y 58664. En la descomposición canónica MPI-4
pertenece solamente a `processor0`, punto de malla local 5331, y aparece en las
caras locales de interfaz 14009, 14019, 14400 y 14412. Es una esquina/arista
superior del patch agudo, no un face centre ni un punto de otro patch.

## Boundary estructural

`solid.msh` identifica el physical surface `interface` con boundary ID 4. La
triangulación productiva contiene 324 caras Q1 en ese boundary. La selección es
físicamente la interfaz correcta: comparte las extensiones globales del canal,
pero sus esquinas superiores están redondeadas.

| ID | physical name | caras | área [m²] | bbox min [m] | bbox max [m] | normal media |
|---:|---|---:|---:|---|---|---|
| 1 | ends | 2000 | 1.5580007801e-6 | (0,-5e-4,0) | (4.92e-3,9.92e-4,5.33e-4) | (1,0,0) |
| 2 | outer | 96 | 1.258536e-5 | (0,-5e-4,0) | (4.92e-3,9.92e-4,5.33e-4) | (0,0,1) |
| 3 | base | 48 | 4.92e-6 | (0,-5e-4,0) | (4.92e-3,9.92e-4,0) | (0,0,-1) |
| 4 | interface | 324 | 2.7416797498e-6 | (0,0,0) | (4.92e-3,4.92e-4,3.3e-5) | (0,0,-1) |

El patch fluido tiene área `2.7453600000e-6 m²`, mientras el boundary 4 tiene
`2.7416797498e-6 m²`. Sus centroides son respectivamente
`(2.46e-3,2.46e-4,3.1048387097e-5)` y
`(2.46e-3,2.4600004621e-4,3.1045493630e-5) m`.

## Punto más cercano y escalas locales

La distancia `5.7090350616544896e-8 m` informada por el bloqueo original era la
distancia a la *bounding box* de la cara, no la distancia geométrica a la cara
Q1. La búsqueda exhaustiva encuentra para la consulta 20:

```text
q_s                 (6.15e-5, 3.464548246917359e-7, 3.276850628713691e-5) m
e = x_f-q_s         (1.355252715606881e-20, -3.464548246917359e-7,
                     2.314937128630913e-7) m
n_s                 (0, 0.831469612302543, -0.555570233019606)
signed normal       -4.166776747647044e-7 m
normal              4.166776747647044e-7 m
tangential          1.402578647394519e-20 m
CellId              2952_0:
face                4
surface quad        315 (Gmsh surface element 1436)
reference diagnostic (0.0375, 0.251697118667249)
minimum location    interior de la cara
distancia a aristas (2.0975314e-7, 1.5785e-3, 6.2360221e-7, 6.15e-5) m
```

La cara es muy anisotrópica: `h_diameter=1.64e-3 m`,
`h_minEdge=8.33355e-7 m`, `sqrt(area)=3.697e-5 m` y aspect ratio
`1967.948`. Por tanto `d/h_diameter=2.54e-4`, pero `d/h_minEdge=0.5` y
`d/sqrt(area)=1.1271e-2`. Usar sólo el diámetro de la cara ocultaba la escala
transversal relevante. Respecto a la precisión de coordenadas, la separación
es aproximadamente `3.8e11` epsilon, por lo que no es redondeo.

## Distribución global y localizador

| dirección | min [m] | media [m] | RMS [m] | p95 [m] | p99 [m] | máximo [m] |
|---|---:|---:|---:|---:|---:|---:|
| fluido → sólido | 0 | 1.8519008e-8 | 8.7843367e-8 | 3.9669e-15 | 4.1667767e-7 | 4.1667767e-7 |
| sólido → fluido | 0 | 1.0475294e-9 | 7.7332929e-9 | 5.9292e-21 | 5.7090351e-8 | 5.7090351e-8 |

Hay 162 puntos fuera de cada umbral absoluto `1e-12`, `1e-10` y `1e-8 m`:
81 en la arista superior izquierda y 81 en la superior derecha. Los otros
3 483 puntos coinciden esencialmente. La búsqueda exhaustiva y la acelerada
eligieron la misma cara en 3 643 puntos; los dos casos restantes son empates
geométricamente equivalentes en caras compartidas. Nunca difirieron en la
distancia mínima. No existe evidencia de `SEARCH_DEFECT`.

## Transformaciones diagnósticas

La mejor traslación `(-2.93e-16,2.55e-9,-1.28e-8) m` deja máximo
`4.075e-7 m`. El ajuste rígido deja `4.096e-7 m`; la similitud obtiene escala
`0.9999984512` y deja `4.093e-7 m`; el escalado por ejes deja `3.906e-7 m`.
Un offset normal constante de `-1.8519e-8 m` deja `3.982e-7 m`. Ninguno elimina
los dos filetes localizados ni tiene una procedencia física compatible. No hay
evidencia de marco u unidades incorrectos.

## Geometría nominal y MappingQ1

El fluido se genera con `blockMeshDict` desde ocho vértices por bloque y paredes
planas. El sólido proviene de `solid.geo`/`roundedSolid.geo`, creado con
OpenCASCADE y filetes circulares superiores. No hay una superficie CAD común:
son geometrías nominales distintas.

La malla sólida usa hexaedros Gmsh tipo 5 y quads de superficie tipo 3, ambos
lineales. La no planaridad máxima de los quads de interfaz es
`3.8808664524e-23 m`, del orden del
redondeo, de modo que la búsqueda diagnóstica sobre sus dos triángulos coincide
con la cara MappingQ1 plana. No contiene nodos geométricos de orden superior que se pierdan al
leerla. `MappingQ1` representa coherentemente esos elementos; la linealización
por cuerdas de un arco es una discrepancia secundaria, no la causa de que el
fluido use una esquina aguda. No corresponde clasificar `MAPPING_MISMATCH`.

## Impacto de una proyección hipotética

No se aplicó proyección. Evaluar `q_s` sólo como diagnóstico dio error relativo
de momento `2.3438259937e-5` para fuerzas aleatorias deterministas. Para
rotaciones unitarias, el error máximo de reproducción fue
`(4.1668e-7,3.4645e-7,3.4645e-7) m` y los errores relativos de trabajo fueron
`(8.6168e-5,1.80795e-4,1.03817e-5)`. Esto viola ampliamente `1e-10` en momento
y `1e-12` en trabajo. Una proyección simple no es una corrección admisible de
G2 sin una formulación conservativa adicional.

## Reproducibilidad y regresiones

La auditoría serial, su repetición, MPI-2 y MPI-4 produjeron clasificación
`CAD_MISMATCH` y hash idéntico
`6efb27a84d0129bf85c5469c0cb878d1a0bb5045f3f059956235544766c70faf`.
La reconstrucción usa únicamente entradas versionadas, por lo que equivale al
restart geométrico sin caches persistidos.

Las regresiones de G2-B.1-RUNTIME pasaron en serial/MPI-2/MPI-4:
`||DeltaR||=1`, `||DeltaJ||=100`, repetición cero, FD `1.02721e-11` y prueba de
signo invertido `2`. EXEC-PROTOCOL pasó en los tres tamaños MPI y conservó sus
controles negativos. Las regresiones legacy mínimo/restart se ejecutaron sin
cambios en código de producción; no se observó NaN/Inf.

## Decisión y siguiente gate

La clasificación primaria única es **CAD_MISMATCH**. La siguiente acción mínima
es un subgate de geometría que seleccione una única superficie nominal y regenere
de forma reproducible ambas mallas coincidentes (o corrija una de ellas según la
geometría física autorizada). No se debe escalar ni proyectar dentro de H_ps.
Sólo después se reejecutará G2-A-RUNTIME-Q1 y, si pasa, se podrá reabrir
G2-B.2-OPERATOR-HANDSHAKE. Esta auditoría no autoriza iniciar G3.
