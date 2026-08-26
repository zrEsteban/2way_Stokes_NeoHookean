# G2-B.1 — protocolo de carga estructural generalizada

Fecha: 2026-08-26

Commit de entrada: `61cc77b4b7e02bfe5ab7aa7bf4c1a03bf37a995a`

Estado: **PASS; producción fluido–sólido permanece legacy**

## Convención residual y tangente

Se eligió la alternativa B. La API recibe dos objetos inequívocos:

```text
f_Gamma [N]
J_Gamma = -df_Gamma/dd [N/m]
R(d) = R_s(d) - f_Gamma(d)
J(d) = J_s(d) + J_Gamma(d)
```

Como `rhs=-R`, `GeneralizedLoad::add_to_newton()` suma `f_Gamma` directamente
a `rhs` y suma `J_Gamma` a la matriz. Mantiene
`f_Gamma=-H^T W_f t_f`. Para la Robin actual `Z_f=z_f I>=0`, la contribución
es `J_Gamma=H^T W_f Z_f H` multiplicada por el coeficiente temporal de la
variable Robin cuando el fluido la suministre. G2-B.1 no fabrica ni activa ese
coeficiente con el fluido real.

La ruta no contiene `FEFaceValues`, `JxW`, áreas ni nearest-neighbour: recibe N
y N/m ya integrados. La ruta legacy en Pa sigue separada e intacta.

## Protocolo y estado

`GeneralizedInterfaceLoad.H` define mensajes neutrales, sin objetos PETSc:

- `schemaVersion`, `transferMode`, `timeIndex`, `outerCorrector`;
- `operatorVersion`, `zVersion`;
- `hashGraph`, `hashWeights`, `dofManifestHash`;
- unidades, flag explícito `valid` y checksum FNV-1a canónico;
- fuerza como `(global_dof_id, component, value_N)`;
- tangente como triplets `(row_dof_id,column_dof_id,value_N_per_m)`.

La recepción es atómica: fuerza y tangente deben tener el mismo sello de
tiempo/corrector/operador. Se rechazan versiones futuras u obsoletas, modo
distinto, hashes, unidades o checksums incompatibles, IDs/componentes fuera de
rango y duplicados. Un vector exactamente nulo con `valid=true` es válido.

Los datos viven sólo en el objeto provisional del corrector. No existe método
de serialización ni promoción automática al estado aceptado; por ello un
restart reconstruye manifiesto/hashes y no puede persistir una carga
provisional. G2-B.2 deberá ligar la creación/destrucción del objeto al protocolo
real y al evento de aceptación temporal de G1.

## Constraints y matriz replicada

La fuerza se transforma como `C^T f`; cada triplet se expande como
`C^T J C` mediante `AffineConstraints::get_constraint_entries()`. Dirichlet
homogéneo se elimina y masters affine reciben sus pesos. Cada entrada expandida
de tangente exige `sparsityPattern.exists(i,j)` antes de `SparseMatrix::add()`.
No se modifica el patrón después de `reinit()`.

Se eligió estrategia MPI A: cada rank recibe el mensaje global completo y lo
ensambla exactamente una vez en su copia de la matriz serial. No hay reducción
de contribuciones. El arnés compara checksums de mensajes, norma de residual y
norma de Jacobiano entre ranks; cualquier diferencia aborta. Esto evita tanto
la multiplicación por `nRanks` como mezclar recepción global y particionada.

## Participante simulado y resultados

Campaña: `/tmp/g2b1-api-20260826`. El participante usa exactamente H/W_f de
G2-A, tracción uniforme en una componente, tracción pseudoaleatoria
determinista en las restantes y `Z_f` variable por cara. También prueba fuerza
nula válida, varios sellos de corrector y mensajes inválidos.

| métrica | serial | MPI-2 | MPI-4 | límite |
|---|---:|---:|---:|---:|
| error fuerza post-constraints | 0 | 0 | 0 | 1e-12 / 1e-10 |
| error acción tangente | 0 | 0 | 0 | 1e-12 / 1e-10 |
| error trabajo virtual | 2,333e-15 | igual | igual | 1e-12 / 1e-10 |
| error FD interfaz | 1,574e-16 | igual | igual | 1e-8 / 1e-7 |
| mejor error FD residual completo | 2,685e-9 | igual | igual | 1e-8 / 1e-7 |
| FD con signo invertido | 2 | igual | igual | debe fallar |

El barrido `epsilon={1e-2,...,1e-12}` del residual completo muestra error de
truncamiento `2,530e-3` en el extremo grande, mínimo `2,685e-9` y contaminación
por redondeo `8,105e-6` en el extremo pequeño. La prueba deliberadamente
invertida da error 2.

La identidad se evalúa en el espacio restringido consistente: `delta_d_s`
satisface constraints, la fuerza estructural es `C^T f` y el trabajo fluido
usa `H C delta_d_hat`. La energía medida con Z no negativa es
`x^T J_Gamma x=0,10014111691516978 >= 0`.

## Regresiones

- G2-A serial/MPI-2/MPI-4: PASS; trabajo máximo `2,471e-15`;
- G2-B.1a Python y C++: PASS; 9.504 entradas, cero faltantes;
- material, BE/BDF2, signo y smoke Robin: PASS;
- `constantHs`: PASS;
- mínimo/restart continuo-partido serial/MPI legacy: PASS;
- serial/MPI legacy `rtol=5e-5`: PASS;
- SHA-256 legacy sin cambios:
  - serial `fsiResiduals`: `8bc7ad2d...80e0`, `robin-out`: `cceaabe1...cc1`;
  - MPI `fsiResiduals`: `7bc98c11...a840`, `robin-out`: `b6e5fa56...2f4a`;
- NaN/Inf: ninguno.

El patrón legacy permanece en 322.992 nonzeros. Ningún caso contiene
`dualConservative`. El ejecutable aborta explícitamente si se intenta ejecutar
ese modo antes de conectar el transporte real; no cae silenciosamente en la
cuadratura legacy. El único warning es la cabecera Boost deprecada transitiva de
deal.II.

## Decisión y G2-B.2

G2-B.1 es **PASS**. G2-B.2 debe implementar exclusivamente el transporte real
del manifiesto y mensajes ya definidos, construir `f_Gamma` y `J_Gamma` en el
participante fluido usando la misma versión de H, aplicar el coeficiente
temporal Robin derivado, y activar explícitamente el caso canónico. Debe
mantener el guard de sellos/checksums y el monitor de trabajo. No requiere
cambiar IQN-ILS, física, tolerancias ni integradores.
