# G2-A-RUNTIME-Q1 — operador de interfaz Q1 productivo

Fecha: 2026-08-26

Commit de entrada: `3e510296513c3e1a38bd35659c3639b35a62699d`

Estado: **BLOCKED_GEOMETRY en la etapa H; no PASS**

Estado corregido de G2-A:

- `G2-A-OFFLINE = PASS`;
- `G2-A-RUNTIME-Q1 = BLOCKED_GEOMETRY`;
- `G2-A-PRODUCTION = INCOMPLETE`.

Los resultados offline siguen siendo válidos como demostración algebraica:
`3520 x 4316`, 45 522 nonzeros, trabajo/fuerza/momento del orden de `1e-15` y
reproducción del orden de `1e-18 m`. No prueban que el soporte coincida con la
interfaz deal.II real.

## Implementación productiva creada

`RuntimeQ1InterfaceOperator.H` define un único `runtime_q1::Builder<dim>` y una
representación canónica de consultas, filas, soportes y hashes. `PdmsSolid`
expone `build_canonical_interface_q1_operator()` después de `read_mesh()` y
`setup()`, por lo que recibe los objetos reales:

- `Triangulation<3>` importada por PdmsSolid;
- `MappingQ1<3>` miembro productivo;
- `DoFHandler<3>` con `FESystem(FE_Q(1),3)`;
- `AffineConstraints<double>` productivas;
- boundary estructural configurado, `4` en el caso canónico.

El constructor indexa exclusivamente caras `at_boundary()` con boundary ID 4.
Para cada consulta usa `MappingQ1::transform_real_to_unit_cell`, verifica la
celda y la coordenada normal de la cara de referencia, vuelve a mapear para
medir el residuo y evalúa `shape_value_component()` del FE real. No contiene
fórmulas trilineales, nearest-neighbour ni proyección. Las alternativas en
aristas/vértices se ordenan por `CellId/faceIndex` y deben producir filas
equivalentes. H_ps permanece en global DoFs completos; no precondensa
`AffineConstraints`.

Los hashes separados cubren consultas, geometría, soporte, grafo y pesos. La
geometría se declara `reference`, versión 1, unidades metro. La tolerancia es
`max(128 epsilon max(1,h_face), relativeTolerance*h_face)` y el API rechaza
tolerancias relativas superiores a `1e-8`.

`testRuntimeQ1` incluye la clase productiva `PdmsSolid`, usa su malla,
DoFHandler, constraints, mapping y boundary ID, y consume el fixture versionado
`dealiiSolid/robin-query.csv`. El fixture contiene únicamente IDs implícitos por
orden canónico y coordenadas; no contiene conectividad, DoFs, CellId ni pesos.

## Primera divergencia geométrica

Compilación de `dealiiPdmsSolid` y `testRuntimeQ1`: PASS. El constructor localizó
las consultas 0--19. La primera consulta que no pertenece a la interfaz real es:

```text
queryPointId       20
x [m]              (6.1500000000000004e-05, 0, 3.3000000000000003e-05)
boundary esperado  4
minBoxDistance     5.7090350616544896e-08 m
closest CellId     2952_0:
face               4
h_face             1.6400002117320407e-03 m
distancia/h_face   3.4811184906037607e-05
tolerancia usada   1.6400002117320407e-13 m
```

La separación relativa `3.48e-5` supera tanto la tolerancia usada (`1e-10 h`)
como el máximo permitido por el gate (`1e-8 h`). No es error de redondeo ni un
punto ambiguo de arista. El operador Python offline lo aceptaba porque buscaba
celdas volumétricas y permitía extrapolación mediante su score. El runtime no
puede reproducir esa conducta bajo `projectionMode=none`.

El estudio explícito `relativeTolerance={1e-12,1e-10,1e-8}` rechazó siempre el
mismo query 20. Las tolerancias físicas efectivas fueron, respectivamente,
`2.842e-14`, `1.640e-13` y `1.640e-11 m`, todas menores que la separación
`5.709e-8 m`; soporte y hashes nunca se publicaron para una fila inválida.

Conforme a la clasificación D, se detuvo la secuencia antes de identidades
globales, MPI, restart, comparación fila a fila y regresiones finales. No se
renormalizaron filas, no se aumentó la tolerancia y no se añadió proyección.

## Resultado y siguiente decisión

G2-A-RUNTIME-Q1 no puede ser PASS con los puntos actuales y
`projectionMode=none`. Se requiere un gate explícito de compatibilidad
geométrica que determine si:

1. el fixture usa un marco/configuración distinto; o
2. las interfaces son no coincidentes y necesitan una proyección superficial
   conservativa explícita.

La segunda alternativa cambia la definición cinemática y debe especificarse y
validarse como tal; no puede ocultarse dentro de RuntimeQ1. Hasta resolverlo,
G2-B.2-OPERATOR-HANDSHAKE permanece pendiente y G2 no está cerrado.

No se modificó `src/robinRobinCoupling`, no se ejecutaron casos FSI duales ni
baseline fuerte, y no se inició G3.
