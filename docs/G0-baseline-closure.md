# Cierre G0: baseline reproducible

Fecha: 2026-08-25 (America/Santiago)  
Baseline: `846a5f0a3f6d3abc9ffc32e83d4c96dc4fe7eb69` (`main` y
`origin/main` antes de este informe)  
Resultado global: **FAIL**

G0 no autoriza cambios en el código ni en el algoritmo FSI. El único cambio
del repositorio producido por este gate es este informe. La validación encontró
un mínimo reproducible correcto hasta `1e-6 s`, pero el baseline configurado
hasta `1e-4 s` abortó reproduciblemente en `1.1e-6 s`; por ello no se aprueba
el gate de estabilidad y no se debe comenzar G1.

## Control de versiones y alcance

Antes de ejecutar:

```text
git branch --show-current       -> main
git rev-parse HEAD              -> 846a5f0a3f6d3abc9ffc32e83d4c96dc4fe7eb69
git rev-parse origin/main       -> 846a5f0a3f6d3abc9ffc32e83d4c96dc4fe7eb69
git status --short              -> salida vacía
git remote -v                   -> git@github.com:zrEsteban/2way_Stokes_NeoHookean.git
```

Después de los builds y casos, y antes de crear este informe,
`git status --short` volvió a ser vacío. Los builds se realizaron en
`/tmp/g0-20260825`; ningún artefacto nuevo apareció en el árbol. Los artefactos
OpenFOAM existentes (`Make/*/*.o`, `*.dep`, `lnInclude`, `platforms`), logs,
resultados, `processor*` y `postProcessing` se muestran como ignorados con
`git status --ignored --short`. No fue necesario modificar `.gitignore`.

Limitación de procedencia importante: la regla `.gitignore:11` ignora todo
`cases/`, por lo que
`cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS` existe en la
máquina, pero **no forma parte del commit 846a5f0** (`git ls-files` no devuelve
archivos para ese caso). El caso tiene además rutas absolutas históricas. Cada
copia de ejecución conservó los parámetros físicos/numéricos y sólo sustituyó
esas rutas por las rutas aisladas de `/tmp`.

## Entorno utilizado

| Componente | Versión o identidad efectiva |
|---|---|
| SO | Ubuntu 22.04.5 LTS; Linux 6.8.0-136-generic x86_64 |
| OpenFOAM | OpenCFD OpenFOAM `v2512`; paquete `openfoam2512 2512.0-2`; `linux64GccDPInt32Opt`, DP, Int32 |
| solids4foam | `v2.3`, commit externo `4b254fa5260e0ae94640d7404089bde73907fc2d`; árbol externo **dirty** |
| deal.II | 9.3.2 (`libdeal.ii-dev 9.3.2-1~exp1ubuntu1`) |
| GCC/G++ | 11.4.0 (`11.4.0-1ubuntu1~22.04.3`) |
| CMake | 3.22.1 |
| Open MPI | 4.1.2 (`libopenmpi-dev 4.1.2-2ubuntu1`) |
| PETSc de solids4foam | 3.25.3, repo `fd2ca03442be2a78130b411057f3b922130f6b1b`, tag `v3.25.3` |
| PETSc enlazado por deal.II | 3.15.5 (`libpetsc_real.so.3.15`) |
| Trilinos enlazado por deal.II | 13.2.0 (Teuchos, ML, Ifpack, Amesos, AztecOO, Epetra, Zoltan, Kokkos) |
| Otras dependencias visibles | HDF5 Open MPI, MUMPS 5.4, p4est 2.2, METIS 5, ScaLAPACK 2.1, Boost 1.74 |

El árbol externo dirty de solids4foam tenía modificaciones y resultados ajenos
a G0. No se tocaron. Esto impide reconstruir sólo desde los commits declarados
la biblioteca `libsolids4FoamModels.so` contra la cual se enlazó el acoplador.

Variables/rutas necesarias:

```text
WM_PROJECT=OpenFOAM
WM_PROJECT_VERSION=v2512
WM_OPTIONS=linux64GccDPInt32Opt
WM_PROJECT_DIR=/usr/lib/openfoam/openfoam2512
WM_PROJECT_USER_DIR=/home/ezamora/OpenFOAM/ezamora-v2512
FOAM_USER_APPBIN=/home/ezamora/OpenFOAM/ezamora-v2512/platforms/linux64GccDPInt32Opt/bin
FOAM_USER_LIBBIN=/home/ezamora/OpenFOAM/ezamora-v2512/platforms/linux64GccDPInt32Opt/lib
MPI_ARCH_PATH=/usr/lib/x86_64-linux-gnu/openmpi
S4F_ROOT=/home/ezamora/Workspace/solids4foam
PETSc runtime=/home/ezamora/Workspace/petsc-3.25.3/arch-openmpi-cpu-opt-nox/lib
```

Los comandos OpenFOAM se ejecutaron mediante `openfoam2512 -c`. Para los casos,
`LD_LIBRARY_PATH` antepuso `/tmp/g0-20260825/lib` y el directorio PETSc 3.25.3
anterior al entorno preparado por OpenFOAM.

## Build controlado

Las fuentes de los dos plugins OpenFOAM se copiaron a
`/tmp/g0-20260825/src`; se aplicó `wclean libso`, se dirigió la salida a
`/tmp/g0-20260825/lib` y se añadió en esas copias la ruta de
`libsolids4FoamModels.so`. La forma canónica en el repositorio sigue siendo
`openfoam2512 -c 'cd src/robinRobinCoupling && wmake libso'`.

Comandos efectivos resumidos:

```bash
openfoam2512 -c 'export FOAM_USER_LIBBIN=/tmp/g0-20260825/lib; cd /tmp/g0-20260825/src/robinRobinCoupling; wclean libso; wmake libso'
openfoam2512 -c 'export FOAM_USER_LIBBIN=/tmp/g0-20260825/lib; cd /tmp/g0-20260825/src/fiveParameterMooneyRivlinElastic; wclean libso; wmake libso'
cmake -S src/dealiiPdmsSolid -B /tmp/g0-20260825/build-dealii -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build /tmp/g0-20260825/build-dealii --parallel 2
ctest --test-dir /tmp/g0-20260825/build-dealii --output-on-failure
```

| Objetivo | Exit | Duración pared | Salida |
|---|---:|---:|---|
| `librobinRobinCoupling.so` | 0 | 25.25 s | 707 KiB |
| `libfiveParameterMooneyRivlinElastic.so` (requerida por `controlDict`) | 0 | 20.07 s | 1.3 MiB |
| configurar deal.II | 0 | 0.33 s | build Release |
| `dealiiPdmsSolid` + `testNeoHookean` | 0 | 16.73 s | 7.6 MiB + 873 KiB |
| `neo_hookean_material` CTest | 0 | 0.15 s | 1/1 PASS |

Un primer intento deliberadamente aislado del acoplador terminó con exit 2 al
no encontrar `-lsolids4FoamModels`; se corrigió en la **copia de `/tmp`**
añadiendo su directorio real y el build limpio posterior terminó con exit 0.
Warnings no fatales: `wmkdepend` no encontró 23 cabeceras opcionales/transitivas
de compatibilidad solids4foam/OpenFOAM; GCC mostró el aviso deprecado de
`tbb/task.h` en los plugins y Boost mostró
`boost/detail/no_exceptions_support.hpp` deprecado en deal.II. No hubo warnings
propios del código ni errores en los builds exitosos.

## Casos y comandos

El candidato es apropiado como caso mínimo de 10 pasos: malla pequeña fija,
`deltaT=1e-7 s`, deal.II BDF2, esquema `semiImplicitBeta=1` con IQN--ILS y
tolerancias físicas explícitas. No es evidencia reproducible desde Git mientras
`cases/` siga ignorado.

Se copiaron y limpiaron directorios independientes. `cleanCase.sh` eliminó
tiempos/resultados y estados deal.II; después sólo se reescribieron rutas
absolutas en cada copia. No se usaron archivos previos como evidencia.

```bash
# Serial; el mínimo cambia sólo endTime a 1e-6 en su copia.
solids4Foam -case /tmp/g0-20260825/runs/minimal
solids4Foam -case /tmp/g0-20260825/runs/baseline-a
solids4Foam -case /tmp/g0-20260825/runs/baseline-b

# MPI declarado por el caso: fluid (2 2 1), solid (4 1 1).
decomposePar -case /tmp/g0-20260825/runs/mpi4-minimal -region fluid -force
decomposePar -case /tmp/g0-20260825/runs/mpi4-minimal -region solid -force
mpirun -np 4 solids4Foam -case /tmp/g0-20260825/runs/mpi4-minimal -parallel
```

## Resultados cuantitativos

| Ejecución | Exit/condición | Pasos aceptados | Tiempo alcanzado | Correctores por paso | Pared |
|---|---|---:|---:|---|---:|
| mínimo serial | 0, `End` | 10 | `1e-6 s` | `1,1,1,1,1,1,1,3,4,4` | 230.91 s |
| baseline A | 134, SIGABRT FSI | 10; falla paso 11 | `1.1e-6 s` | igual al mínimo; luego 100 | 1532.51 s |
| baseline B | 134, SIGABRT FSI | 10; falla paso 11 | `1.1e-6 s` | igual al mínimo; luego 100 | 1531.24 s |
| mínimo MPI, 4 ranks | 0, `End` | 10 | `1e-6 s` | `1,1,1,1,1,1,1,3,4,4` | 224.55 s |

El formato de `/usr/bin/time` imprimió `exit=0` al ser la aplicación terminada
por señal, pero el shell/runner devolvió inequívocamente 134 y registró
`Command terminated by signal 6` para A y B.

En el mínimo serial, el último corrector tuvo residuo geométrico relativo
`7.807517485e-5`, `max|r_u|=1.541915375e-7`, defecto de potencia relativo
`9.67608909e-7` y ratio de aceptación `0.5088320737`; todos pasan sus gates.
En ambos baselines completos el último corrector tuvo residuo geométrico
`5.254843908e-5 < 1e-4`, pero `max|r_u|=9.735517512e-7` da ratio físico
`3.212720779 > 1`; al llegar a 100 correctores el solver aborta. No es una
divergencia NaN/Inf sino un incumplimiento reproducible del gate físico.

Normas sobre las 3645 muestras de `robin-out.csv`:

| Ejecución/estado | max desplazamiento (m) | max velocidad (m/s) | max aceleración (m/s2) | max tracción (Pa) |
|---|---:|---:|---:|---:|
| mínimo serial aceptado | `2.101379583e-13` | `9.892006932e-7` | `5.338689470` | `199.1533633` |
| mínimo MPI aceptado | `2.101372875e-13` | `9.891944940e-7` | `5.338643196` | `199.1534989` |
| trial fallido baseline A/B | `4.256686703e-13` | `3.082259442e-6` | `32.51879265` | `632.2268432` |

Los RMS seriales del mínimo son, respectivamente,
`6.297507644e-14 m`, `3.070778798e-7 m/s`, `1.694320606 m/s2` y
`23.12918385 Pa`. El trial fallido se informa para diagnosticar el aborto; no
se trata como estado aceptado.

No se encontraron tokens numéricos NaN/Inf en ninguno de los cuatro logs.

## Repetibilidad y tolerancias

Para dos salidas sobre las mismas coordenadas se ordenan las filas por
`(x,y,z)`. Para cada vector se calcula
`||a-b||_inf` (máximo de la norma euclídea por muestra) y RMS. Se aprueba si
`||a-b||_inf <= atol + rtol*max(||a||_inf,||b||_inf)`, con `rtol=5e-5` y
`atol=1e-15 m` para desplazamiento, `1e-12 m/s` para velocidad,
`1e-7 m/s2` para aceleración y `1e-6 Pa` para tracción. Estas tolerancias son
mucho menores que los valores físicos y permiten diferencias de orden/reducción
MPI sin exigir identidad binaria.

Baseline A/B: `fsiResiduals.dat` y `robin-out.csv` son byte a byte idénticos;
la condición, paso, 100 correctores y métricas finales también coinciden. Los
logs completos difieren sólo por identidad/ruta de caso y tiempos de pared.

Serial/MPI en el mínimo:

| Magnitud | norma infinito de diferencia | RMS diferencia | error relativo al máximo | Resultado |
|---|---:|---:|---:|---|
| desplazamiento | `1.1123e-18 m` | `2.4983e-19 m` | `5.29e-6` | PASS |
| velocidad | `7.3426e-12 m/s` | `1.4445e-12 m/s` | `7.42e-6` | PASS |
| aceleración | `6.0163e-5 m/s2` | `1.5076e-5 m/s2` | `1.13e-5` | PASS |
| tracción | `3.8563e-3 Pa` | `6.6747e-5 Pa` | `1.94e-5` | PASS |
| residuo FSI final | `2.3086e-10` | n/a | `2.96e-6` | PASS |

MPI terminó con `End`, 10 pasos y el mismo vector de correctores. La equivalencia
serial/MPI pasa dentro de la tolerancia definida.

## Tabla del gate

| Criterio obligatorio | Estado | Evidencia/razón |
|---|---|---|
| Control de versiones | PASS | baseline, rama, remoto y árbol inicial confirmados |
| Build del acoplador | PASS | build limpio controlado, exit 0 |
| Build del sólido deal.II | PASS | configure/build exit 0; CTest 1/1 |
| Caso mínimo | PASS | 10 pasos, `End`, gates físicos y geométrico aprobados |
| Baseline estable | **FAIL** | aborta en paso 11, ratio físico 3.212720779 |
| Dos ejecuciones equivalentes | PASS | mismo fallo, residuos y CSV idénticos |
| Sin NaN/Inf | PASS | cero ocurrencias en cuatro logs |
| Artefactos correctamente ignorados | PASS | `git status --short` vacío tras ejecutar |
| Comandos y versiones documentados | PASS | secciones anteriores |
| Serial/MPI | PASS | 4 ranks, descomposición declarada, `End`, errores < `5e-5` relativos |

## Riesgos y limitaciones abiertos

1. El baseline configurado no es estable más allá de 10 pasos: FAIL bloqueante.
2. El caso está fuera de Git por `cases/`; el commit baseline por sí solo no lo
   reproduce.
3. solids4foam es un repositorio externo dirty; su binario/biblioteca no queda
   identificado únicamente por el commit anotado.
4. `Make/options` fija `S4F_ROOT` a una ruta personal y el caso fija rutas
   absolutas; la portabilidad depende de reescrituras controladas.
5. PETSc tiene dos versiones efectivas: solids4foam usa 3.25.3 y deal.II enlaza
   la 3.15.5 del sistema; no produjo fallo aquí, pero debe conservarse explícito.
6. Los avisos `wmkdepend` y deprecaciones TBB/Boost no impiden el build, pero
   reducen la limpieza del diagnóstico.

## Siguiente gate propuesto (no ejecutado)

No comenzar G1. El siguiente paso debe ser un **G0 de remediación y
revalidación**, primero haciendo versionables/portables el caso y el manifiesto
del entorno externo, y luego corrigiendo o justificando la pérdida de
convergencia física en el paso 11. Sólo después deben repetirse exactamente los
ensayos seriales y MPI de este informe y exigirse `End` en `1e-4 s`.
