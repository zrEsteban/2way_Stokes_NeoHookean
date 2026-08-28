# G2-A-RUNTIME-Q1-EDGE-CANONICALIZATION

Fecha: 2026-08-28

Commit de entrada: `5a641962862c1a2fd1232dbd35172d383fbffdd9`

Clasificación: **PASS_RUNTIME_Q1_WITH_COST_RISK**

## Causa y corrección

La causa primaria fue `REFERENCE_ROUNDOFF`, no una discontinuidad topológica.
Para el ID histórico 17289, `x=(2.5625e-6,0,3.225e-5) m`, las dos caras
incidentes fueron:

| CellId/cara | xi crudo | xi canónico | soporte final por componente |
|---|---|---|---|
| `2781_0:/0` | `(1.8378817905537546e-18,0.9984375,7.408279391249113e-16)` | `(0,0.9984375,0)` | `(12564,0.0015625),(12567,0.9984375)` |
| `2952_0:/4` | `(-3.5822150435305296e-12,0.9984375,1.1592065929460622e-13)` | `(0,0.9984375,0)` | igual |

Antes de canonicalizar, la segunda evaluación dejaba pesos espurios de orden
`3.58e-12` en nodos fuera de la arista. Después, ambas filas son idénticas,
`rowSum=1` y el residuo físico es `3.6422416731934915e-20 m`.

La corrección actúa sólo en referencia después de aprobar pertenencia física.
Una coordenada se fija en 0/1 si su cambio no supera

```text
physicalTolerance/minEdge + 512 epsilon max(1, diameter/minEdge)
```

y si el remapeo con `MappingQ1` continúa dentro de la tolerancia física
original. No se cambió esa tolerancia, no se proyectó el punto y no se
renormalizaron filas. Las shape functions Q1 reales se reevalúan en `xi`
canónico; la masa eliminada por poda fue exactamente cero.

## Equivalencia y desempate

Para cada componente se construye el funcional disperso y su expansión real
por `AffineConstraints`, `d=Cq+d0`. Antes del desempate se verifican:

```text
(rA-rB)                 = 0 dentro del error derivado del mapping
(rA-rB) C               = 0
(rA-rB) d0              = 0
C^T (rA-rB)^T           = 0
```

Esto cubre cinemática, offset inhomogéneo y carga dual condensada. Sólo entonces
se selecciona menor `CellId` y menor cara. No se promedian filas. Un índice
R-tree de bounding boxes limita candidatos; ya no existe un barrido
punto-por-todas-las-caras.

El barrido completo produjo 69 012 interiores de cara, 17 397 aristas y 36
vértices. Hubo 11 604 puntos multicandidato, 11 628 comparaciones equivalentes
y cero incompatibles. Se canonicalizaron 127 100 coordenadas de candidatos;
el máximo cambio paramétrico fue `1.2307021881907588e-10`. El test analítico
cubre cara, arista, vértice, perturbación de redondeo, FESystem y rechazos por
duplicado, interior volumétrico, boundary, marco y tolerancia inválidos. Los
guardas conservan el fallo ante soporte, componente, constraint o carga
transpuesta incompatibles.

## QueryManifest y operador productivo

`g2aFluidPointManifest` es una utilidad OpenFOAM C++ que lee el `polyMesh` real,
usa `pointProcAddressing`, deduplica por ID global y entrega el mismo manifiesto
en serial, MPI-2 y MPI-4. No lee CSV ni la malla estructural. Se fijó
`writePrecision 17` porque `decomposePar` con precisión 10 alteraba las
coordenadas de los processor meshes; esta corrección preserva, no modifica, la
superficie maestra.

| métrica | resultado |
|---|---:|
| puntos OpenFOAM | 86 445 |
| filas vectoriales H_ps | 259 335 |
| nonzeros | 932 634 |
| `dofManifestHash` | `8b1befa8bcfed1fe` |
| `queryManifestHash` | `6aedd17ba49b68ed` |
| `supportHash` | `594135f2e0ccd49e` |
| `hpsGraphHash` | `a3e23641fde5ff88` |
| `hpsWeightsHash` | `bf73bc060f5e98e0` |
| `geometryHash` | `13da653088ddcec6` |

Serial, MPI-2, MPI-4 y reconstrucción repitieron exactamente hashes, tamaños y
acciones. `constraints.distribute(d_s)` se aplica antes de H_ps; el operador
permanece expresado en DoFs globales completos y la ruta runtime aprobada trata
la transpuesta sin una segunda expansión.

## Resultados numéricos y coste

| control | resultado | límite |
|---|---:|---:|
| máximo `|sum w-1|` | `2.220446049250313e-16` | `1e-14` |
| reproducción lineal | `6.035753494226803e-16 m` | `1e-12 m` |
| rotación rígida | `4.228388472693467e-16 m` | `1e-12 m` |
| trabajo dual | `8.450856094048475e-14` | `1e-12` serial |
| fuerza | `2.244854113527965e-14` | `1e-12` |
| momento | `1.8107744147170954e-14` | `1e-12` |
| residuo geométrico | `6.035759726770217e-16 m` | tolerancia original |

El manifiesto tarda 2.30 s y 149 572 KiB RSS; H_ps completo tarda 1.87 s y
300 744 KiB RSS. El barrido 1/8, 1/4, 1/2 y completo está en
`tests/g2ARuntimeQ1/scalability.csv`. La malla tiene 24 veces más celdas y
23.72 veces más puntos de interfaz que la histórica. El coste del solve
OpenFOAM/ALE queda como riesgo no medido porque este gate no ejecuta FSI.

El framing admite 64 MiB. Las estimaciones son 3.11 MiB para
FluidPointManifest, 34.6 MiB para HpsRowManifest, 0.41 MiB para estado
estructural, 0.21 MiB para fuerza y 0.35 MiB para tangente: no se requiere
chunking en la discretización actual, pero HpsRowManifest conserva margen
limitado para refinamientos futuros.

## Regresiones y estado

- CTest: 5/5 PASS (material, BE/BDF2, signo, codec y aristas Q1).
- G2-B.1-RUNTIME: serial/MPI-2/MPI-4 PASS, FD `1.02721e-11`, repetición 0.
- EXEC-PROTOCOL: serial/MPI-2/MPI-4 y 27 negativos PASS.
- La suite Python offline de sparsity no escala con 460 800 celdas y se detuvo;
  se ejecutó la regresión C++ equivalente con el manifiesto previamente
  aprobado. El oracle offline permanece histórico y no construye este H_ps.
- `src/robinRobinCoupling` no cambió; no se ejecutó FSI ni G2-B.2.

Estado final: `G2-A-RUNTIME-Q1=PASS`, `G2-A-PRODUCTION-HPS=PASS`,
`G2-B.2-OPERATOR-HANDSHAKE=READY_TO_REOPEN`, `G2 global=OPEN`.
