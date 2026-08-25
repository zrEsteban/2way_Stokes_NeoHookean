# G0-MPI-D — diagnóstico de la discrepancia serial/MPI

Fecha: 2026-08-25 (`America/Santiago`)

## Resultado ejecutivo

**Diagnóstico: PASS. Clasificación principal: B — contaminación o error de
ejecución/descomposición.** No se modificó código fuente, formulación FSI,
parámetro físico, tolerancia, discretización, transferencia, IQN-ILS ni
condición Robin.

La procedencia del build aislado es correcta. El fallo no lo causan MPI,
PETSc, la transferencia ni IQN-ILS: durante `decomposePar`, el entorno aislado
carga correctamente `librobinRobinCoupling.so`; al escribir el campo
descompuesto, `pdmsElasticWallPressureFvPatchScalarField::write()` omite
`constantHs`. Cada `processor*/0/fluid/p` queda con `coeff0`, `coeff1`, `rhs`,
`prevPressure`, pero sin `constantHs`.

En runtime, el valor ausente toma el default `-1`; la condición usa
`hs = ap*dt` en vez del espesor explícito `3.469008593e-6 m`. Así duplica
exactamente la masa superficial y la impedancia:

| ejecución | `constantHs` en processor | masa/área [kg/m2] | impedancia [Pa s/m] |
|---|---:|---:|---:|
| serial aislado | sí | 0.003364938335 | 33649.38335 |
| MPI antiguo | sí | 0.003364938335 | 33649.38335 |
| MPI aislado original | **no** | 0.006729876670 | 67298.76670 |
| MPI aislado, `decomposePar -no-libs` | sí | 0.003364938335 | 33649.38335 |

`decomposePar -no-libs` preserva el diccionario como tipo desconocido y
elimina la discrepancia sin cambiar ningún dato físico ni algoritmo. Es una
remediación del procedimiento de descomposición, no una estabilización.

## Identidad y alcance

- Repositorio principal: `/home/ezamora/OpenFOAM/2way_Stokes_NeoHookean`.
- HEAD auditado: `6af1f9f0ca9d21823c9e740a9d863843b898b467`.
- Commit canónico de fuentes: `846a5f0a3f6d3abc9ffc32e83d4c96dc4fe7eb69`.
- solids4foam upstream: `4b254fa5260e0ae94640d7404089bde73907fc2d`.
- Representación fijada: commit base anterior más los parches/overlay
  versionados bajo `external/solids4foam`; checkout reconstruido en
  `external/solids4foam/source`.
- Prefijo aislado:
  `platforms/g0e-846a5f0/linux64GccDPInt32Opt`.
- Caso: `cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS`.
- Control temporal diagnóstico: `endTime=1e-6`, cambiado sólo en copias bajo
  `/tmp`; `deltaT=1e-7` no cambió.
- MPI: 4 rangos, descomposición `(2 2 1)` para fluido y `(4 1 1)` para el
  modelo OpenFOAM sólido; `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`.
- El participante deal.II es siempre un proceso serial separado invocado por
  el rango master. La implementación no soporta sólido deal.II MPI; por ello
  sólo existen las combinaciones soportadas fluido serial+sólido serial y
  fluido MPI+sólido serial.

## Ejecuciones nuevas

Todas las copias se obtuvieron mediante `git archive HEAD`, regeneraron las
mallas y comenzaron sin tiempos posteriores a `0`, `processor*`,
`postProcessing`, CSV ni estados deal.II previos.

| ID | fluido | descomposición | resultado | pasos/correctores | exit |
|---|---:|---|---|---|---:|
| `isolated-a-clean` | MPI 4 | `decomposePar` con libs | fallo en `8e-7` | `1,1,1,1,2,2,5,100` | launcher 1 / MPI_ABORT |
| `isolated-b-clean` | MPI 4 | igual | fallo en `8e-7` | `1,1,1,1,2,2,5,100` | launcher 1 / MPI_ABORT |
| `control-serial` | serial | no aplica | `End` en `1e-6` | `1,1,1,1,1,1,1,3,4,4` | 0 |
| `control-mpi-nolibs` | MPI 4 | `decomposePar -no-libs` | `End` en `1e-6` | `1,1,1,1,1,1,1,3,4,4` | 0 |
| `control-mpi-first` | MPI 4 | `-no-libs`, un paso | `End` en `1e-7` | `1` | 0 |

Hubo un intento descartado (`isolated-b-src`): el arnés llamó al ejecutable
sin hacer `cd` al caso y falló en el primer corrector al no poder crear
`dealiiSolid/robin-in.csv`. No se usó como evidencia y se creó otra copia.

Comandos esenciales (dentro de OpenFOAM v2512 y después de activar el entorno):

```bash
export PETSC_DIR=/home/ezamora/Workspace/petsc-3.25.3
export PETSC_ARCH=arch-openmpi-cpu-opt-nox
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source env/activate-project.sh

CASE=/tmp/g0-mpi-d-20260825/isolated-a-clean
decomposePar -case "$CASE" -region fluid -force
decomposePar -case "$CASE" -region solid -force
cd "$CASE"
mpirun -np 4 "$FOAM_USER_APPBIN/solids4Foam" -case "$CASE" -parallel

CASE=/tmp/g0-mpi-d-20260825/control-mpi-nolibs
decomposePar -no-libs -case "$CASE" -region fluid -force
decomposePar -no-libs -case "$CASE" -region solid -force
cd "$CASE"
mpirun -np 4 "$FOAM_USER_APPBIN/solids4Foam" -case "$CASE" -parallel
```

## Reproducibilidad del fallo MPI aislado

Las dos repeticiones válidas son numéricamente y byte a byte reproducibles:

| métrica final | repetición A | repetición B |
|---|---:|---:|
| tiempo | `8e-7` | `8e-7` |
| correctores | 100 | 100 |
| `max|r_u|` [m/s] | 1.642341764e-6 | 1.642341764e-6 |
| ratio físico | 5.419727823 | 5.419727823 |
| residuo geométrico relativo | 2.869738984e-5 | 2.869738984e-5 |
| residuo geométrico absoluto | 2.965508636e-13 | 2.965508636e-13 |
| defecto relativo de potencia | 2.172840841e-5 | 2.172840841e-5 |
| NaN/Inf | ninguno | ninguno |
| exit | 1 / MPI_ABORT | 1 / MPI_ABORT |

SHA-256:

```text
185e7ea45a51042d531df1d6f54b42f3ef318e1e47f46986e81e54291235aab7  fsiResiduals.dat (A y B)
14e915d474c6a6fe03c0c6bf6147e1b87044b67a7429460b29a56b5b03cdfd6d  robin-out.csv final (A y B)
```

Se archivaron 113 intercambios por repetición, cada uno con `robin-in.csv`,
`robin-out.csv` y `robin-query.csv` (339 archivos). Los conjuntos finales y
los historiales FSI compararon con `cmp=0`.

## Primera divergencia

La primera divergencia es de **entrada/configuración**, al terminar la
descomposición y antes del primer solve:

```text
0/fluid/p:            constantHs 3.469008593e-06
processor0/0/fluid/p: constantHs ausente; coeff0/coeff1/rhs escritos
```

El primer desacuerdo numérico aparece en `t=1e-7`, corrector 1:

| métrica | MPI antiguo/correcto | MPI aislado defectuoso | factor |
|---|---:|---:|---:|
| `max|r_u|` | 1.148873823e-9 | 4.452056315e-9 | 3.8751 |
| ratio físico | 0.003791283617 | 0.01469178584 | 3.8751 |
| residuo geométrico | 3.817299070e-7 | 1.480223271e-6 | 3.8772 |

La captura de los CSV del primer corrector, ordenada por coordenadas, da las
siguientes normas entre el control `-no-libs` y el caso defectuoso:

| cantidad | máximo | RMS |
|---|---:|---:|
| desplazamiento [m] | 3.330544218e-16 | 9.195046411e-17 |
| velocidad [m/s] | 3.330544218e-9 | 9.195046411e-10 |
| aceleración [m/s2] | 3.330544218e-2 | 9.195046411e-3 |
| tracción de salida [Pa] | 12.18507214 | 1.185412927 |
| tracción transferida de entrada [Pa] | 12.18507214 | 1.082110918 |

La velocidad Robin de entrada aún es idéntica (diferencia cero), mientras la
tracción ya difiere. Con las 3,520 caras y el área uniforme de cara
`2.0295e-10 m2`, la fuerza transferida del primer intercambio es:

```text
control:    (6.442564425e-10, -9.268060502e-15, 6.460861721e-8) N
defectuoso: (1.078297317e-9,  -9.093115736e-15, 1.288234916e-7) N
```

La componente dominante casi se duplica (`1.994x`), coherente con la
impedancia duplicada. La primera divergencia precede a la creación de historia
IQN-ILS útil.

## Comparación old/isolated

El árbol externo original y el checkout reconstruido tienen el mismo HEAD
`4b254fa5260e0ae94640d7404089bde73907fc2d`. La comparación recursiva de
`src`, `applications`, `etc` y `Allwmake`, excluyendo plataformas y artefactos,
no produjo diferencias. No hay archivos no rastreados bajo esos inputs de
build. Los seis archivos modificados necesarios son exactamente los seis
representados por los parches versionados; lo demás del árbol externo son
resultados/logs incidentales de tutoriales.

Los `Make/*`, flags efectivos (`linux64GccDPInt32Opt`, DP, Int32, Opt) y líneas
de enlace coinciden. La comparación de símbolos, strings, secciones de datos y
desensamblado normalizado no encontró diferencias semánticas en las funciones
participantes. Las diferencias de SHA en algunas bibliotecas provienen de rutas
de compilación/build-id:

| artefacto | antiguo SHA-256 | aislado SHA-256 | conclusión |
|---|---|---|---|
| `solids4Foam` | `31c342343f...e2e41` | igual | idéntico |
| `libsolids4FoamModels.so` | `3c65a7a5de...34b33` | `ad3e0764c3...7228f` | semánticamente idéntico |
| `librobinRobinCoupling.so` | `366dfdf52c...633bf` | `1196c2f04a...f92df` | sólo ruta de fuente embebida |
| `dealiiPdmsSolid` | `0cde590160...e03f5` | `0513de624a...ae4d` | sólo ruta de build embebida |

`solids4Foam`, `libsolids4FoamModels.so` y `librobinRobinCoupling.so` no tienen
RPATH/RUNPATH. `dealiiPdmsSolid` tiene RUNPATH del sistema para HDF5/OpenMPI.
La activación aislada y los mapas runtime confirman que los objetos propios se
cargan desde el prefijo aislado; no aparece `ezamora-v2512`.

La diferencia histórica decisiva está en los logs de descomposición:

- antiguo: `decomposePar` no pudo cargar los plugins porque
  `libpetsc.so.3.25` no estaba disponible; OpenFOAM trató la condición como
  genérica y preservó `constantHs`;
- aislado: el plugin cargó sin warning, se construyó el tipo real y su método
  `write()` no escribió `constantHs`.

Los meshes iniciales y descompuestos compararon con SHA idéntico; los inputs
`0` y `system` también, salvo el campo `p` reserializado descrito arriba.

## Auditoría de convergencia y reducciones

Definiciones exactas:

- residuo geométrico: desplazamiento deal.II menos desplazamiento de puntos
  fluido; la norma absoluta usa `sqrt(gSum(...))` y, con el floor activo, el
  criterio relativo usa `gMax`;
- ratio físico: máximo de residuo relativo de velocidad, tracción, Robin,
  potencia local y potencia integrada, todos construidos con `gMax/gSum`;
- aceptación: residuo geométrico `<= outerCorrTolerance` **y** ratio físico
  `<= 1`;
- aborto: ambas condiciones siguen fallando al alcanzar `nOuterCorr=100`.

Todas las reducciones revisadas usan el comunicador mundial de OpenFOAM. Los
global patches contienen vectores globales replicados. Por ello `gSum` de una
norma absoluta cuenta la réplica una vez por rango: la norma absoluta MPI es
`sqrt(4)=2` veces la serial. No altera el criterio activo, que usa `gMax`, ni
los cocientes donde el mismo factor aparece en numerador y denominador. Los
cuatro rangos imprimen/reciben los mismos escalares colectivos.

El valor pedido `2.869738984e-5 / 5.419727823` no proviene de una reducción
local: es la respuesta reproducible al `constantHs` perdido y a la impedancia
duplicada. Al preservar el campo, el MPI termina.

## Transferencia e IQN-ILS

- Cada rango tiene 880 caras de interfaz fluida y 640 del patch sólido
  OpenFOAM; el patch global contiene 3,520 caras fluidas. El CSV deal.II tiene
  3,645 puntos y el sólido reporta 12,948 DoF.
- AMI reporta error geométrico `1.738109131e-18` tanto en antiguo como aislado.
- `patchFaceToGlobal/globalFaceToPatch` construyen los vectores globales antes
  de las normas y transferencias.
- IQN-ILS usa `gSum` para productos internos y normas, `gMax` para limitar la
  corrección, mantiene historia por interfaz y reinicia/persiste según la
  configuración existente. El factor uniforme de réplica cancela en QR y en
  los coeficientes.
- En el primer corrector sólo se usa el startup factor; la discrepancia ya
  existe antes de acumular modos. No hay evidencia de categoría D o E.
- No se añadió instrumentación al código. La captura externa de CSV no cambió
  el orden de operaciones ni el algoritmo.

## PETSc y MPI observados en runtime

Se capturaron `/proc/<pid>/maps` durante `isolated-b-clean`:

- cada uno de los cuatro procesos `solids4Foam` cargó únicamente
  `/home/ezamora/Workspace/petsc-3.25.3/arch-openmpi-cpu-opt-nox/lib/libpetsc.so.3.25.3`;
- el proceso separado `dealiiPdmsSolid` cargó únicamente
  `/usr/lib/x86_64-linux-gnu/libpetsc_real.so.3.15.5`;
- todos cargaron `/usr/lib/x86_64-linux-gnu/libmpi.so.40.30.2`; deal.II añadió
  `libmpi_cxx.so.40.30.1` y los rangos OpenFOAM `libmpi_mpifh.so.40.30.0`;
- ninguna dirección de proceso contiene ambas versiones PETSc;
- ninguna biblioteca propia se cargó desde `ezamora-v2512`.

Por tanto no hay mezcla ABI PETSc ni incompatibilidad MPI demostrada.

## Validación de la remediación procedimental

Ambos controles alcanzaron `End`, 10 pasos y la misma secuencia de correctores.
Último corrector:

| métrica | serial | MPI `-no-libs` | diferencia relativa |
|---|---:|---:|---:|
| `max|r_u|` | 1.541915375e-7 | 1.541882474e-7 | 2.13e-5 |
| ratio físico | 0.5088320737 | 0.5088212164 | 2.13e-5 |
| residuo geométrico | 7.807517485e-5 | 7.807540571e-5 | 2.96e-6 |
| defecto relativo de potencia | 9.676089090e-7 | 9.675901179e-7 | 1.94e-5 |
| condición | End, exit 0 | End, exit 0 | equivalente |

Comparando `robin-out.csv` por coordenadas, las diferencias máximas relativas
son `5.61e-6` en desplazamiento, `7.74e-6` en velocidad, `2.87e-5` en
aceleración y `1.94e-5` en tracción. Se adopta tolerancia diagnóstica
`rtol=5e-5`, `atol=1e-12` para cinemática y `atol=5e-3 Pa` para tracción; todas
las magnitudes pasan. No se exige igualdad binaria entre serial y MPI.

SHA-256 de los controles:

```text
e2d3baab1e2844ffe74ef3611c7c66d563dc47aaa9127c594aaa9cf59204a7ee  serial/fsiResiduals.dat
456fcad9861aa7d41f369967798ef1fcc113267b9c81b7868c32e6fae454c472  mpi-nolibs/fsiResiduals.dat
2ad4330e6c0c71d4572ee2b2783e11dce6a7380f8a9e13984c27a98fa8cbf527  serial/robin-out.csv
9fe82c68db1d2e2b22b42d3efe01ae10c49376f06f6be082db4c5134bb14793a  mpi-nolibs/robin-out.csv
```

## Clasificación y decisión recomendada

| categoría | resultado | evidencia |
|---|---|---|
| A procedencia/build | descartada | fuentes compiladas y semántica binaria coinciden |
| **B ejecución/descomposición** | **confirmada** | `write()` pierde `constantHs`; `-no-libs` lo preserva y MPI pasa |
| C reducción/métrica | descartada como causa | métricas colectivas; valor problemático desaparece con input correcto |
| D transferencia | descartada como causa | discrepancia anterior a transferencia; AMI consistente |
| E IQN-ILS | descartada como causa | divergencia en corrector 1 antes de historia útil |
| F sensibilidad | descartada como causa primaria | perturbación determinista de impedancia 2x, no redondeo |

De acuerdo con la regla del gate, **G0 debe permanecer bloqueado mientras la
preparación MPI oficial no fije el procedimiento correcto**. Se recomienda una
remediación G0-E-R estrictamente procedimental: documentar/versionar que ambos
`decomposePar` se ejecutan con `-no-libs` (o, en un gate posterior apropiado,
hacer que `write()` preserve `constantHs`) y repetir el cierre G0. La evidencia
de este informe demuestra que esa remediación recupera el baseline mínimo MPI;
no justifica iniciar G1 ni modificar la física.

## Archivos modificados y riesgos abiertos

Archivo del repositorio modificado por G0-MPI-D:

- `docs/G0-MPI-discrepancy-diagnosis.md` (este informe).

No hubo instrumentación ni cambios de fuente. Riesgos abiertos:

1. cualquier futura descomposición que cargue la condición real volverá a
   perder `constantHs` hasta que el procedimiento use `-no-libs` o se corrija
   la serialización en el gate autorizado;
2. las normas absolutas `gSum` sobre global patches replicados escalan con
   `sqrt(nProcs)`; hoy no controlan la aceptación, pero debe conservarse esta
   advertencia si cambia el floor o el criterio;
3. no existe combinación sólido deal.II MPI en la arquitectura actual; no se
   improvisó una;
4. el baseline fuerte no se volvió a ejecutar en este gate porque la causa
   está confinada a `decomposePar` y el fuerte serial ya es reproducible; su
   inestabilidad física conocida continúa asignada a G7.

**Siguiente acción propuesta, sin ejecutarla aquí:** remediación documental del
procedimiento MPI en G0-E-R y repetición del cierre G0. No iniciar G1 hasta esa
decisión.
