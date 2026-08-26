# G2-B.1a — ampliación dispersa del patrón de interfaz

Fecha: 2026-08-25  
Commit de entrada: `ac47ff29d5a3270c3789046238250aa6cad2c5bb`  
Estado: **FAIL parcial — grafo/ensamblaje simulado PASS; integración C++ pendiente**

## Derivación del grafo

Para cada fila r de H se construye sólo el soporte disperso
`S_r={i:H_ri!=0}` y se añaden pares ordenados `S_r x S_r`. No se forma una
matriz densa. El hash canónico del grafo escalar requerido es:

```text
hashGraph = 6ee94c285c4e5d8619f9f08c3a666eb57dd46ddfe47b2b02be25c1248b4abc39
```

La condición Robin productiva vigente usa `solid impedance` escalar y el
producto vectorial por componentes, por lo que `Z_f=z_f I`: sólo conecta cada
componente consigo misma. Una impedancia diagonal vectorial conserva esa
estructura. En cambio, un tensor general o `z n tensorProduct n` contiene
términos fuera de diagonal para una normal genérica; el patrón seguro debe
incluir los nueve bloques de componentes por par nodal.

El arnés registra ambos tamaños, pero usa el superset tensorial para dimensionar
la ampliación segura. No extrae Z de `H^T W_f Z_f H`.

## Constraints

El mallado Gmsh importado es fijo y no contiene hanging nodes. Las superficies
físicas 1, 2 y 3 imponen Dirichlet homogéneo sobre 2 256 nodos. Por tanto

```text
d = C d_hat + d0,  d0=0,
```

y C elimina esos componentes; no introduce masters affine adicionales. El
grafo efectivo se forma como `C^T K C`, eliminando simbólicamente cualquier par
con un nodo restringido. No se inserta sobre DoFs Dirichlet.

Para una futura malla con hanging nodes, la implementación C++ deberá expandir
cada índice por `AffineConstraints::get_constraint_entries()` antes de añadir
el producto cartesiano al `DynamicSparsityPattern`; condensar después sin
preparar ese fill-in no es suficiente.

## Resultados del participante simulado

Campaña: `/tmp/g2b1a-sparsity-20260825-1`.

| métrica | valor |
|---|---:|
| nodos Q1 | 4 316 |
| filas H | 3 520 |
| pares escalares faltantes sin constraints | 13 124 |
| pares escalares efectivos faltantes | 3 168 |
| nonzeros vectoriales legacy efectivos | 316 224 |
| añadidos para Z escalar/diagonal | 9 504 |
| añadidos para Z tensorial/`n⊗n` | 28 512 |
| nonzeros finales, superset tensorial | 344 736 |
| crecimiento estimado | 9,0163934426 % |
| conexiones faltantes tras ampliar | 0 |

La estimación de memoria usa el número de entradas del patrón; bytes exactos
dependen del backend y sus índices. El soporte, no los pesos, determina
`hashGraph`; cambios de H/W conservando soporte actualizan `hashWeights` pero no
reconstruyen el patrón. Un cambio de `hashGraph` durante Newton se rechaza.

## Ensamblaje disperso de prueba

Se ensambló K como diccionario disperso de triplets y se comparó con
`H^T(W_f Z_f(Hx))`, usando Z uniforme y variable. Cada rank posee las filas de
cara `globalFaceId modulo nRanks` y las contribuciones se acumulan con suma MPI.

| prueba | serial | MPI-2 | MPI-4 |
|---|---:|---:|---:|
| acción, Z uniforme | 4,894e-16 | 2,826e-16 | 2,705e-16 |
| acción, Z variable | 5,155e-16 | 2,940e-16 | 2,885e-16 |
| simetría, Z uniforme | 5,273e-17 | 5,016e-17 | 5,127e-17 |
| simetría, Z variable | 5,489e-17 | 5,443e-17 | 5,377e-17 |
| `x^T K x`, Z uniforme | 0,031031467429314213 | equivalente | equivalente |
| `x^T K x`, Z variable | 0,027328668240122086 | equivalente | equivalente |

Los modos aleatorios son reproducibles (`seed=20260825`). Tras constraints,
traslaciones/rotaciones que violan Dirichlet no pertenecen al espacio libre;
las pruebas cinemáticas sin constraints permanecen cubiertas por G2-A.

## Pruebas negativas

El participante simulado rechaza:

- manifiesto ausente en modo dual;
- `hashGraph` distinto;
- ID fuera de rango;
- owner negativo/inválido;
- componente fuera de 0--2;
- inserción dual sobre el patrón legacy;
- cambio topológico durante Newton;
- triplet fuera del patrón extendido.

## Por qué el gate todavía no es PASS

El sólido productivo actual ejecuta deal.II serial con
`dealii::SparseMatrix<double>`; no usa una matriz PETSc distribuida. Además:

- `PdmsSolid::setup()` aún crea y finaliza el patrón antes de recibir un
  manifiesto H;
- no existe parser C++ del manifiesto canónico ni mapeo de sus IDs a los IDs
  reales del `DoFHandler`;
- la ampliación demostrada vive en el participante simulado, no en el
  `DynamicSparsityPattern` productivo;
- no puede activarse el error de nuevas asignaciones PETSc en un backend que no
  es PETSc;
- no se ha comprobado todavía `AffineConstraints` con masters/hanging nodes en
  C++.

Por ello no se cumplen dos criterios obligatorios: manifiesto antes de
inicializar la matriz y ausencia de nuevas asignaciones en la matriz productiva.
No se presenta el test de triplets Python como sustituto de esa evidencia.

## Inicialización C++ requerida

El paso pendiente debe separar mínimamente:

1. `dofs.distribute_dofs(fe)`;
2. construir `AffineConstraints`;
3. exportar IDs DoF de interfaz;
4. leer/verificar manifiesto y mapear IDs;
5. crear el DSP FE y expandir constraints de cada clique H;
6. finalizar el patrón;
7. inicializar `SparseMatrix`/vectores;
8. bloquear `hashGraph` durante Newton.

Legacy debe omitir 3--5 y producir exactamente los 316 224 nonzeros efectivos
actuales. Dual debe abortar si el manifiesto no existe antes del paso 6.

## Archivos y decisión

- `scripts/audit-g2b1-sparsity.py`: auditor simbólico y ensamblaje disperso;
- `scripts/test-g2b1a-sparsity.sh`: serial, MPI-2 y MPI-4;
- este informe y actualizaciones documentales.

No se modificó código C++, casos ni binarios; producción sigue legacy. Por ello
no se repitieron G0/G1 ni los casos legacy: no existe refactor productivo que
pueda introducir una regresión.

Estado final: **G2-B.1a FAIL parcial**. Siguiente paso: integrar el mismo grafo
en un `DynamicSparsityPattern` C++ previo a la matriz y probar la prohibición de
nuevas entradas. G2-B.1 no se completa; G2-B.2 y G3 no se inician.
