# G1-R — primera divergencia de restart

Fecha: 2026-08-25  
Commit de partida: `2a736a1341d7e23e4e28e6b3e348f4fbeb176542`  
Estado del diagnóstico: **PASS (divergencia localizada)**  
Estado de G1: **no se cierra en este gate**

Este documento registra trabajo nuevo de G1-R. No reevalúa ni resume el informe
general de G1.

## Ejecución instrumentada

Se añadió el arnés versionado
`scripts/diagnose-g1-restart-first-divergence.sh`. Cada rama se crea con
`git archive`, usa `writePrecision 17` y no reutiliza resultados. La campaña
válida quedó en:

```text
/tmp/g1r-first-divergence-20260825-4
```

Se instrumentó temporalmente y bajo `G1_RESTART_DIAGNOSTICS=1` el
`pimpleFluid.esi.C` de la revisión fijada de solids4foam. Se registraron los
campos antes de construir la ecuación de presión y, para cada corrección no
ortogonal, `diag`, `lower`, `upper`, `source`, `internalCoeffs` y
`boundaryCoeffs` de `pEqn`. La instrumentación se retiró después de capturar la
evidencia y se recompiló solamente `libsolids4FoamModels.so`. Su SHA-256
restaurado es:

```text
ad3e0764c3bd66d798b5036667c1f1d453c7370d247dadf3cc2efcc684c7228f
```

Ambas ramas parten del mismo estado aceptado en `t=5e-7`:

- continua: prosigue en memoria hasta `t=6e-7`;
- restart: termina en `t=5e-7`, abre una shell/ejecución nueva, relee
  `latestTime` y avanza hasta `t=6e-7`.

Ambas ejecuciones terminaron con código 0 y `End`.

## Antes de escribir y después de releer

La instantánea escrita en `5e-07` se comparó antes de iniciar el restart y de
nuevo después de releer. `diff -qr` no encontró diferencias. Los campos que
entran al primer paso reiniciado son idénticos byte a byte:

| objeto | SHA-256 común |
|---|---|
| `fluid/U` | `14b0b4b04c0e7b7e7eecf25f36e4f7a2f1dfbb245b4a8e08f47d228bd0ae4364` |
| `fluid/U_0` | `caa493361f83a72fd716d15d1d5c8f064acadf0adb86df8af12eb8095f54da7f` |
| `fluid/Uf` | `4036340392f280bf9d994a888e1bf62a55274db876da4ff3b49b74a6638ea097` |
| `fluid/Uf_0` | `5f78452f9e3dcbb92b96e8569f93733dc5dacad023ecc0e98042a7dc6a4bf3fa` |
| `fluid/p` | `b81613173ad8e1881386938e710f9938ee250abdd4e0e8f65e24dad30909f633` |
| `fluid/p_0` | `87e3226922cc39024a730bb039fc7fa06bd428c52328d60ddf2f25366b4dd831` |
| `fluid/phi` | `96603460125ed37e3fa00918f66fd5deada7df5ec1c259ebf1eac3346b0de20d` |
| `fluid/meshPhi` | `cc7409532bb7f0f0386e0d256bd7677a2ff58cea8ee129300d2e9b82fe7fd9c1` |
| `fluid/V0` | `d60a24c4a44f7b029d8891a915508b78cf54ef60ce0ef174c025172b687802c9` |
| `fluid/polyMesh/points` | `6ad2dca5f1207ad7fe1b3781f9813d707b9b4671dcd41ff7280d5727ecfed2e2` |

Esto descarta una pérdida textual de `U`, historias, presión, flujos, volumen o
puntos como primera divergencia.

## Primera divergencia antes del solve de presión

En la primera corrección de `t=6e-7`, antes de `pEqn.solve`, coinciden:

```text
sum(U)       = (0.090685915649390839 1.6426376600340486e-08 0.049003809255311885)
sum(U.old)   = (0.070983272305534917 9.3712694931419133e-09 0.04924112467782403)
sum(U.oldOld)= (0.052261294328790585 3.6608536322184459e-09 0.049284256186904504)
sum(phi)     = 4.5775781393451396e-11
sum(A(U))    = 325490451071.93152
sum(H(U))    = (1321323.1202293925 0.19782806706585263 851465.51759385155)
sum(HbyA)    = (0.07807633995746692 1.165184297424663e-08 0.049024689685730559)
sum(phiHbyA) = 5.1125017138242406e-11
sum(p)       = 19.983843758584925
```

También son idénticos `diag`, `lower`, `upper`, `source` e
`internalCoeffs` de la primera `pEqn`. El **primer operando diferente antes del
solve** es `pEqn.boundaryCoeffs()[3]`, correspondiente al patch físico
`interface` (los patches 0--2 son `inlet`, `outlet` y `glass`). Sus 3520
componentes difieren:

```text
max(abs(delta)) = 2.4291277027267406e-17
```

El máximo relativo no es una métrica útil aquí porque varios coeficientes
cruzan cero y son del orden de `1e-21`; alcanza `1.4242732204915736`. Un ejemplo
del primer coeficiente es:

```text
continua = -3.4497349110955422e-14
restart  = -3.4497349660886193e-14
```

Los agregados ocultaban esta diferencia: `sum(diag)` y `sum(source)` eran
iguales, y `sum(source^2)` sólo variaba en el último bit. Después del primer
solve, el primer campo observable diferente es `U`; en la corrección PIMPLE
siguiente:

```text
sum(U) continua = (0.093223648517699098 9.4925226650602063e-09 0.050123790950080703)
sum(U) restart  = (0.093223657626961873 6.3100426900675629e-09 0.050123450602561244)
```

La primera residual inicial GAMG es `0.059691256396702426` (continua) frente a
`0.059624512273709143` (restart). Por tanto la divergencia no comienza en la escritura
de los campos ni después del solve: comienza al reconstruir los coeficientes de
borde Robin de la matriz de presión a partir del estado releído/geometría
derivada.

## Clasificación A--G

Para hacer inequívoca la clasificación usada por G1-R:

| clase | causa | resultado |
|---|---|---|
| A | campo físico no escrito o valor textual alterado | descartada |
| B | historia `oldTime` incompleta o desplazada | descartada |
| C | parámetro/estado persistente Robin perdido | descartada en los diccionarios escritos |
| D | estado geométrico o coeficiente de borde derivado reconstruido de forma no idempotente | **principal** |
| E | matriz interna o RHS de presión diferente | descartada antes de incorporar el borde |
| F | estado/caché interno del solver lineal | no es la primera divergencia |
| G | indeterminismo o procedencia distinta | descartada |

Clasificación obligatoria: **D**. La hipótesis queda delimitada a la
reconstrucción del coeficiente Robin del patch `interface` (normales,
`deltaCoeffs` y los operandos derivados usados por `updateCoeffs`). Los puntos y
los miembros persistidos son idénticos, pero el coeficiente ensamblado no lo
es. No se modificó la fórmula ni se intentó corregirla en G1-R.

## Archivos y comandos

Archivo nuevo:

- `scripts/diagnose-g1-restart-first-divergence.sh`
- `docs/G1-restart-first-divergence.md`

Comandos principales:

```bash
bash -n scripts/diagnose-g1-restart-first-divergence.sh
G1R_WORK=/tmp/g1r-first-divergence-20260825-4 \
  scripts/diagnose-g1-restart-first-divergence.sh

openfoam2512 -c 'export PETSC_DIR=/home/ezamora/Workspace/petsc-3.25.3 \
PETSC_ARCH=arch-openmpi-cpu-opt-nox; source env/activate-project.sh; \
wmake libso external/solids4foam/source/src/solids4FoamModels'
```

El último comando recompiló el objetivo afectado después de retirar la
instrumentación; no recompiló otros objetivos ni alteró el algoritmo.

## Riesgo y siguiente acción

La perturbación absoluta es pequeña, pero el acoplamiento la amplifica después
del primer solve. Falta instrumentar directamente los operandos por cara dentro
de `pdmsElasticWallPressureFvPatchScalarField::updateCoeffs()` para separar
`nf`, `deltaCoeffs`, `prevAcceleration` y el término Robin. Esa corrección o
decisión de persistencia pertenece a una remediación posterior de G1. **No se
inicia G2.**
