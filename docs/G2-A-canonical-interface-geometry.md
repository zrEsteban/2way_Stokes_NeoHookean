# G2-A-GEOMETRY-CANONICALIZATION

Fecha: 2026-08-27

Commit de entrada: `8e07a7c295af9de317d98a324c9621fad252a3e3`

Clasificación primaria: **BLOCKED_MESH_QUALITY**

El subgate no es PASS. Se consiguió una interfaz discreta exactamente
coincidente y un generador determinista, pero `checkMesh -allGeometry
-allTopology` rechaza la calidad volumétrica cerca de los filetes. La malla
candidata no se instaló en el caso canónico. No se reintentó el handshake, no se
ejecutó FSI dual ni baseline fuerte y no se inició G3.

## Autoridad física y rama elegida

El filete es intencional, no un accidente del generador:

- `cases/pdmsMicrochannelFSI/README.md` declara dos filetes circulares
  longitudinales superiores de `0.75 um` y describe su mallado refinado;
- `AUDITORIA_CASO_PDMS_MICROCHANNEL_FSI.md` lo registra como dimensión del
  canal y del sólido;
- `common/parameters.json` contiene `upperCornerRadius=7.5e-7 m`;
- `roundedFluid.geo` y `roundedSolid.geo` implementan los mismos centros y
  transiciones desde su incorporación en `6af1f9f0`.

Se adoptó **FILLETED_MASTER**. La superficie discreta maestra es boundary 4 de
`solid.msh`, evaluada como caras `MappingQ1`; no es el círculo analítico. El
sólido existente se conserva porque representa el contrato y tiene Jacobianos
positivos. El caso sharp-corner queda como fixture histórico de
`G2-A-interface-geometry-audit`.

## Contrato y generador

`geometry/interface-geometry-v1.json` fija marco `reference`, unidades SI,
ejes, dimensiones, radio, centros, ángulos, extensiones, nombres, boundary 4,
MappingQ1, `projectionMode=none`, discretización y tolerancias. Su hash canónico
es:

```text
b6f0e724f2b12893eb6e49f24904a070ef3e37b44d64ab5c455b6893ff46319b
```

`generate_canonical_interface_geometry.py` genera primero el sólido desde el
contrato. Luego extrae del `solid.msh` las aristas del boundary 4 en el plano de
entrada, elimina sólo vértices intermedios estrictamente collineales y conserva
cada cambio de normal Q1. La sección fluida transfinita usa esas cuerdas,
24×10 celdas y 80 capas axiales. Ninguna fórmula circular independiente entra
en la superficie fluida.

El generador escribe además, sin uso runtime, la identidad de punto, elemento
maestro, coordenada física y versión geométrica. No edita manualmente
`polyMesh`, `solid.msh` ni coordenadas. Herramientas usadas: Gmsh 4.8.4,
`gmshToFoam`/OpenFOAM v2512 y Python 3.

Dos generaciones limpias (`/tmp/g2a-canonical-final-a-20260827` y
`...-final-b-20260827`) fueron idénticas:

| objeto | hash |
|---|---|
| manifiesto | `dc5b1dd28d3bd46d14d2e206db5a3542f21d42525338aacdd5ba98fed8638ae9` |
| solid.msh | `dace8367d5cab024d4136e3b1ffa350cdb479ca888568bba93cdda40a467dd9b` |
| fluid.msh | `9193216a77e97cce8cbf1969af79776320533282a1d94fb628350b0967595cc1` |
| polyMesh semántico | `236d5e67f88eefdd435fc20d97da82109fd2bba0a5a429c5f0eaa1688789209b` |

## Coincidencia de la superficie discreta

La validación incluye todos los vértices, centros y puntos internos de cada cara
fluida, transiciones y distancias bidireccionales.

| métrica | resultado |
|---|---:|
| puntos fluidos de interfaz | 3645 |
| caras fluidas de interfaz | 3520 |
| caras Q1 maestras | 324 |
| máximo fluido→sólido | `0 m` |
| máximo sólido→fluido | `0 m` |
| máximo en muestras interiores | `3.3881317890e-20 m` |
| error relativo de área | `3.6455712154e-14` |
| error de `n_f ≈ -n_s` | `2.2204460493e-16` |
| gaps / overlaps | `0 / 0` |
| caras que cruzan transición no coplanar | `0` |
| no planaridad máxima del master | `3.8808664524e-23 m` |

El área fluida es `2.7416797497924303e-6 m²` y el área maestra
`2.7416797497923304e-6 m²`. Cada cara fluida queda dentro de una cara maestra o
de una unión estrictamente coplanar. No se proyectó, extrapoló ni aumentó una
tolerancia.

## Calidad y bloqueo

La mejor variante conservó 19 200 hexaedros, 22 275 puntos y volúmenes
positivos. `checkMesh` obtuvo:

| métrica | resultado |
|---|---:|
| volumen mínimo | `1.0236947626e-16 m³` |
| aspect ratio máximo | `133.0091` |
| non-orthogonality máxima/media | `66.2841° / 10.1299°` |
| skewness máxima | `2.52003` |
| peso de interpolación mínimo | `0.0642329` |
| determinante mínimo | `6.5128809158e-6` |
| celdas con determinante `<0.001` | `652` |
| resultado | `Failed 1 mesh checks` |

Se probaron por separado tetraedrización, recombinación y malla transfinita
suavizada. La última elimina los fallos de non-orthogonality, skewness, pesos y
volume ratio, pero no el determinante bajo. La causa es la combinación de
cuerdas Q1 submicrométricas (`2.926e-7 m`) con la longitud del canal. Aprobarla
ignorando `checkMesh` contradiría el gate.

La malla estructural conserva 3000 hexaedros, 4316 puntos, boundary IDs y
MappingQ1, sin Jacobianos negativos ni elementos degenerados en Gmsh/deal.II.

## Impacto geométrico

Respecto del canal sharp-corner:

- volumen fluido: `7.988112e-11 → 7.9879200913e-11 m³`, cambio relativo
  `-2.40243e-5`;
- área de interfaz: `2.74536e-6 → 2.7416797498e-6 m²`, cambio relativo
  `-1.34053e-3`;
- sólido, DoFHandler y constraints: sin cambios;
- conteos fluidos: se conservan 19 200 celdas, 3520 caras y 3645 puntos de
  interfaz, aunque cambian sus coordenadas y la topología interior.

Los campos físicos sharp y filleted no deben compararse como una regresión de
solución.

## RuntimeQ1 y pruebas posteriores

La secuencia obligatoria se detuvo en calidad antes de aprobar RuntimeQ1,
MPI/restart y regresiones FSI. Como diagnóstico adicional, el operador original,
sin modificar y con su tolerancia intacta, localizó la superficie pero rechazó
la primera transición compartida con:

```text
incompatible interface rows at edge/vertex for queryPointId 729
```

El inverso MappingQ1 entrega residuos de coordenada de referencia de orden
`1e-13` en las dos celdas anisotrópicas adyacentes. Una prueba temporal confirmó
que canonicalizar coordenadas de referencia ya aceptadas elimina esa diferencia,
pero dicho cambio se retiró: este gate exigía probar el operador sin modificar.
Debe auditarse en un gate separado después de resolver la calidad.

No se ejecutaron mínimos legacy sobre la geometría candidata porque ésta no fue
aceptada ni instalada. Las mallas y casos legacy versionados permanecieron
intactos; `src/robinRobinCoupling` no cambió.

## Decisión y siguiente corrección mínima

Clasificación única: **BLOCKED_MESH_QUALITY**.

El siguiente subgate debe diseñar refinamiento volumétrico local conforme en la
dirección axial de las columnas de los filetes, con transición poliedral/hex
válida, sin cambiar la superficie maestra. Debe conseguir `Mesh OK` antes de
reabrir RuntimeQ1. Después corresponderá resolver separadamente la
canonicalización numérica de trazas Q1 compartidas, repetir RuntimeQ1
serial/MPI/restart y sólo entonces declarar el handshake `READY_TO_REOPEN`.

Estados finales:

- `G2-A-OFFLINE = PASS histórico`;
- `G2-A-GEOMETRY-CANONICALIZATION = BLOCKED_MESH_QUALITY`;
- `G2-A-RUNTIME-Q1 = BLOCKED_GEOMETRY`;
- `G2-A-PRODUCTION = INCOMPLETE`;
- `G2-B.2-HANDSHAKE = BLOCKED`;
- `G2 global = OPEN`.
