# G2-A-FILLET-MESH-QUALITY

Fecha: 2026-08-27

Commit de entrada: `0c39a3783a8134539babdab8192d4cf4cc4b7406`

Clasificación: **PASS_MESH_QUALITY**

El gate se limitó a la calidad volumétrica de `FILLETED_MASTER`. No cambió el
contrato, el radio `0.75 um`, la superficie Q1 maestra, los límites de
`checkMesh`, `RuntimeQ1InterfaceOperator` ni `src/robinRobinCoupling`. No se
ejecutaron RuntimeQ1 completo, FSI dual, G2-B.2, baseline fuerte ni G3.

## Reproducción y causa espacial

La generación limpia original reprodujo exactamente 19 200 hexaedros,
`minDet=6.512880915810036e-6`, 652 celdas bajo `0.001` y `Failed 1 mesh
checks`. Topología, orientación, regiones, volúmenes, pirámides, concavidad,
non-orthogonality, skewness, pesos y volume ratio pasaron.

El mapa completo está en `tests/g2AFilletedMeshQuality/bad-cell-map.csv`. Las
652 celdas se descomponen en ocho columnas axiales completas (640 celdas) y
seis secciones adicionales en cada extremo axial. Los índices transversales
dominantes son `6--9` y `236--239`, cuatro por cada transición superior; todos
los 80 índices axiales aparecen. No coinciden con un orden de vértices
incorrecto ni una interfaz interna defectuosa. La causa primaria es **A+C**:
columnas completas que heredan la escala submicrométrica de las cuerdas Q1 del
filete y una longitud axial inicial de `61.5 um`. Los extremos añaden el efecto
secundario de la transición geométrica.

Histograma acumulado original:

| umbral | celdas |
|---:|---:|
| `1e-5` | 4 |
| `5e-5` | 8 |
| `1e-4` | 168 |
| `5e-4` | 414 |
| `1e-3` | 652 |
| `2e-3` | 980 |
| `5e-3` | 1472 |
| `1e-2` | 1960 |

El máximo cociente de volumen entre vecinos es `26.4646507861` (p99
`4.3750123001`). La geometría usa un único bloque transfinito extruido; por
conformidad, una columna axial no puede subdividirse de forma aislada sin
introducir hanging faces o una topología de transición nueva. Se mantuvo la
topología hex conforme y se refinó la coordenada axial común.

## Candidatos y selección

`scripts/generate_canonical_interface_geometry.py` acepta overrides explícitos
de discretización y registra `effectiveDiscretization` y
`discretizationHash`; el `geometryContractHash` permanece inmutable. El runner
versionado es `scripts/run-g2a-fillet-mesh-candidates.sh`. La tabla completa,
incluidos p01/p05/p50, tiempos y hashes, está en `candidate-summary.csv`.

Resultados decisivos:

| candidato | celdas | minDet | bajo 0.001 | fallos |
|---|---:|---:|---:|---:|
| original axial 80 | 19 200 | `6.513e-6` | 652 | 1 |
| axial 320 | 76 800 | `9.708e-5` | 330 | 1 |
| grading 1.0 | 19 200 | `4.243e-6` | 164 | 1 |
| transversal 48x20 | 76 800 | `2.008e-6` | 5232 | 2 |
| axial 640 | 153 600 | `3.538e-4` | 4 | 1 |
| axial 1280 | 307 200 | `1.179e-3` | 0 | 0 |
| **axial 1920** | **460 800** | **`2.2237626204e-3`** | **0** | **0** |

El grading mejora non-orthogonality/skewness pero no el peor determinante; el
refinamiento transversal lo empeora y crea skewness fallida. Axial 1280 es el
Pareto mínimo que pasa, pero quedó descartado por margen bajo. Axial 1920 es la
primera variante medida que supera el objetivo recomendado `0.002`; no depende
del redondeo del umbral.

## Calidad de la candidata instalada

| métrica | resultado |
|---|---:|
| celdas / puntos | `460800 / 528275` |
| volumen min/max | `4.2653948441e-18 / 2.6360594602e-16 m3` |
| determinant min | `0.00222376262038637` |
| p01 / p05 / p50 determinant | `0.0535927 / 0.0631625 / 0.126812` |
| aspect ratio max | `8.70066499` |
| non-orthogonality max/media | `66.2840773 / 10.1086331 deg` |
| skewness max | `2.52003162` |
| face weight min | `0.0642329183` |
| face volume ratio min | `0.0377862534` |
| concave/warped/illegal | `0 / 0 / 0` |
| topología/geometría fallidas | `0 / 0` |
| `checkMesh` serial/MPI-2/MPI-4 | `Mesh OK / Mesh OK / Mesh OK` |

No hay volúmenes negativos ni caras casi nulas. Quedan sólo cuatro celdas bajo
`0.005`; ninguna está cerca de `0.001`. El aumento de 24 veces en celdas es el
riesgo abierto principal: eleva memoria/coste y debe observarse bajo ALE.

## Superficie, reproducibilidad y MPI

La interfaz refinada tiene 86 445 puntos y 84 480 caras. Todos provienen de la
misma superficie Q1; no hay proyección. La validación obtuvo distancia
fluido--sólido y sólido--fluido cero, máximo interior `3.3881e-20 m`, error de
área `3.4633e-13`, error de normales opuestas `2.2204e-16` y cero
gaps/overlaps/cruces no coplanares.

Dos generaciones limpias produjeron el mismo hash semántico. La candidata fue
instalada exclusivamente con:

```text
openfoam2512 -c 'python3 scripts/generate_canonical_interface_geometry.py \
  --output /tmp/g2a-fillet-selected-install \
  --fluid-axial-cells 1920 --install'
```

Hashes:

- contrato: `b6f0e724f2b12893eb6e49f24904a070ef3e37b44d64ab5c455b6893ff46319b`;
- discretización: `c7c9080ba929e5470f5ff498a07947c3e610a3df11ec4e7f1bb4f67026b6aa33`;
- polyMesh semántico: `8781ced412365874c791a49d8fe8b055a1f25436032b9ad5e425370540994554`;
- manifiesto: `b4104b2b8c5f7bbf84a34524817103ecd36680b55aef2bda96661e836f75ae9c`.

MPI-2 y MPI-4 conservaron 460 800 celdas globales, el mismo mínimo y `Mesh
OK`; no hubo pérdida de addressing ni modificación de la interfaz. La prueba
aislada de motion solver se omitió porque no hay una amplitud aprobada para
este subgate; no se inventó una carga extrema.

## Regresiones y estado

`G2-B.1-RUNTIME` pasó serial/MPI-2/MPI-4: `DeltaR=1`, `DeltaJ=100`, error FD
`1.02721e-11`, repetición `0`. EXEC-PROTOCOL pasó serial/MPI-2/MPI-4 y sus 26
escenarios negativos usando el manifiesto aprobado prevalidado. El oracle
Python offline no se vuelve a calcular sobre las 86445 consultas nuevas: esa
ruta no es el RuntimeQ1 productivo y escala de forma prohibitiva. El protocolo,
el código legacy y `robinRobinCoupling` no cambiaron.

Estados al cierre:

- `G2-A-FILLET-MESH-QUALITY = PASS`;
- `G2-A-GEOMETRY-CANONICALIZATION = PASS_FILLETED_MASTER`;
- `G2-A-RUNTIME-Q1 = PENDING_RETRY`;
- `G2-A-PRODUCTION = INCOMPLETE`;
- `G2-B.2-HANDSHAKE = BLOCKED`;
- `G2 global = OPEN`.

El siguiente paso autorizado es reintentar G2-A-RUNTIME-Q1 sobre la malla
instalada, sin relajar el operador. No se ejecutó en este gate.
