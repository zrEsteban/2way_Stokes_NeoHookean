# G1-RF — causalidad del coeficiente Robin

Fecha: 2026-08-25  
Base: `44ffe699d861fb66e505a3a763c677eb538cd91d`  
Resultado: **PASS; defecto determinista corregido (rama A)**

## Primer suboperando

La instrumentación opt-in descompuso, por cara, la condición
`pdmsElasticWallPressure` y la primera ecuación de presión en `t=6e-7`.
`deltaT`, `constantHs`, `deltaCoeffs`, `weights`, `Sf`, `magSf`, `nf`, `Cf`,
`Cn`, campo interno, `rhoSolidHs` y `prevPressure` coinciden bit a bit. El
primer suboperando diferente es **`prevAcceleration`**:

```text
continuo: 3520 vectores no nulos del estado aceptado
restart:  3520{(0 0 0)}
```

El archivo `5e-07/fluid/p` sí contiene los 3520 vectores. El constructor los
lee: `sumPrevAcceleration=(29.476769898233119, 0.40325589618756358,
202.28858225586043)`. Una llamada posterior los sustituía por cero.

Para el máximo de `pEqn.boundaryCoeffs()[3]`:

| dato | valor |
|---|---|
| índice local de patch | 2235 |
| patch/cara | `interface`, cara global 59275 (`startFace=57040`) |
| continuo | `4.7766321317626014e-17` |
| restart | `2.3475044290358608e-17` |
| diferencia absoluta | `2.4291277027267406e-17` |
| diferencia relativa | `0.5085440192419381` |
| continuo IEEE-754 | `0x1.b89109930eb8bp-55`, bits `0x3c8b89109930eb8b` |
| restart IEEE-754 | `0x1.b109c328c215cp-56`, bits `0x3c7b109c328c215c` |
| distancia de representaciones | `4636040960985647` ULP discretos |
| `magSf` | `1.2607499999999982e-09`, bits `0x3e15a8d64ebbb796` |
| `deltaCoeff` | `606060.6060606063`, bits `0x41227ed9364d9367` |

No es una perturbación de un ULP. En esa cara:

```text
prevAcceleration continuo = (0.07085446353264084 -0.002219997709090627 0.5113493122548148)
prevAcceleration restart  = (0 0 0)
prevDdtUn continuo         = 0.5095581089590642
prevDdtUn restart          = 0
rhs continuo               = -3.374558807830753e-06
rhs restart                = -1.6599271929772835e-06
valueBoundaryCoeffs cont.  = -1.1102872379960343e-06
valueBoundaryCoeffs rest.  = -5.461442764276437e-07
```

`valueFraction`, `refValue` y `refGrad` no existen como miembros de esta
implementación Robin. Sus equivalentes se construyen con:

```text
coeff0 = 1
coeff1 = rhoSolidHs/rhoFluid
rhs = prevPressure - rhoSolidHs*(nf & prevAcceleration)
valueInternalCoeffs = coeff1/(coeff1 + coeff0/deltaCoeffs)
valueBoundaryCoeffs = (rhs/deltaCoeffs)/(coeff1 + coeff0/deltaCoeffs)
```

El laplaciano combina esos coeficientes con la interpolación de `rAtU` y la
geometría del patch para formar `pEqn.boundaryCoeffs()[3]`.

## Orden, estado y caches

En ambos recorridos hubo seis llamadas efectivas a `updateCoeffs()` en el paso
diagnosticado (continuo 30--35; proceso nuevo 0--5), siempre con `updated()==0`
al entrar y con el mismo orden de evaluación. La geometría continua incremental
y la reconstruida desde `polyMesh/points` produjeron listas bitwise idénticas.

La secuencia defectuosa era:

1. el constructor del patch lee `prevAcceleration` correctamente;
2. deal.II aún no ha enviado `dealIIFaceAcceleration` en el proceso nuevo;
3. `updatePdmsElasticWallPressure()` inicializa un campo local a cero;
4. aun sin dato disponible, sobrescribe el estado aceptado restaurado;
5. `updateCoeffs()` ensambla el RHS Robin con aceleración cero.

No hubo cache geométrico obsoleto ni diferencia en `movePoints()`.

## Experimentos causales A--D

| experimento | resultado aislado |
|---|---|
| A: geometría canónica reconstruida | El restart ya realiza esa reconstrucción. Frente al continuo incremental, puntos y todas las cantidades geométricas fueron bitwise iguales; no elimina la divergencia de aceleración. |
| B: misma secuencia de estado Robin | Conservar el aceptado hasta disponer de dato deal.II hace bitwise idénticos todos los suboperandos y `boundaryCoeffs`; residual GAMG y trayectoria coinciden. |
| C: eliminar movimiento de malla | La descomposición por operandos lo vuelve innecesario como discriminante: `Sf`, `magSf`, `nf`, `Cf`, `Cn`, `weights` y `deltaCoeffs` ya coinciden con movimiento activo. No se alteró el caso físico para obtener el PASS. |
| D: incremental frente a puntos finales | Comparación directa bit a bit: ninguna diferencia geométrica; se descarta reconstrucción geométrica como causa. |

No se acumularon cambios en el experimento causal B: la única modificación fue
evitar la asignación sintética de cero cuando el campo externo aún no existe.

## Corrección

`robinRobinCouplingInterface::updatePdmsElasticWallPressure()` ahora marca la
aceleración como disponible sólo cuando deal.II entregó un campo no vacío (o
cuando el sólido OpenFOAM proporcionó el campo). Únicamente entonces actualiza
`prevAcceleration`. En un restart deal.II temprano conserva el estado aceptado
leído. La fórmula Robin y su orden algebraico no cambian.

## Amplificación

Antes de la corrección, en el primer corrector de `t=6e-7`:

| etapa | continuo | restart | perturbación relativa |
|---|---:|---:|---:|
| `prevDdtUn[2235]` | 0.5095581089590642 | 0 | 1.0 |
| `rhs[2235]` | -3.374558807830753e-6 | -1.6599271929772835e-6 | 0.5081 |
| `boundaryCoeffs[2235]` | 4.7766321317626014e-17 | 2.3475044290358608e-17 | 0.5085 |
| residual inicial GAMG | 0.059691256396702426 | 0.059624512273709143 | 0.001118 |
| residuo FSI, corrector externo 1 | 5.1475174054381677e-5 | 5.9826546357700895e-5 | 0.1622 |

El cociente dimensional entre la diferencia del residual GAMG y la del
coeficiente es `2.747658055126e12`. Sólo hubo un corrector externo por paso en
este caso; ambos aceptaban uno, de modo que no existe bifurcación por número de
correctores. La diferencia final provenía de estado reconstruido incorrecto y
su amplificación por el acoplamiento, no de una tolerancia o aceptación
discreta. Después de la corrección, la cadena de perturbación es exactamente
cero en los archivos comparados.

## Regresiones

Campaña: `/tmp/g1rf-restart-regression-20260825-1`.

- serial continuo 10 pasos y serial 5+5: `End`, BDF2 tras restart, 10
  aceptaciones, `robin-out.csv` y `fsiResiduals.dat` combinado idénticos;
- MPI continuo y MPI 5+5: mismos resultados de restart;
- serial/MPI: PASS con `rtol=5e-5`;
- ausencia de NaN/Inf: PASS;
- serial `robin-out.csv`:
  `cceaabe1995fc8bf45b148a7375f0b632e5afb2c4baf6312c3ad1243c0a14cc1`;
- MPI `robin-out.csv`:
  `b6e5fa5648291d58911bca8b6fa0eef08794758f0f486f8b8676f45d03d42f4a`;
- serial `fsiResiduals.dat` continuo:
  `8bc7ad2dd28176344230a6bef559688b6784aab19935a2b415effeea0aa080e0`;
- MPI `fsiResiduals.dat` continuo:
  `7bc98c11c6423702e2a948c0ac6dea7db2b09ca9665efa3d45885e3322e0a840`;
- round-trip Robin y descomposición `constantHs`: PASS;
- build final de `librobinRobinCoupling.so`: exit 0.

La instrumentación pesada fue retirada. El árbol sólo conserva la corrección
causal y documentación. Estado final: **G1 PASS**. Siguiente gate recomendado:
G2, sin iniciarlo aquí.
