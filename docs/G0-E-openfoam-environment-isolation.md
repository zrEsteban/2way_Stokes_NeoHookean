# G0-E-R — Remediación de procedencia y aislamiento OpenFOAM/solids4foam

Fecha: 2026-08-25 (America/Santiago)  
Baseline canónico: `846a5f0a3f6d3abc9ffc32e83d4c96dc4fe7eb69`

HEAD inicial: `d41eaaef2b2c59a3daa0e2142f75eaab618dc0b9`
Resultado: **FAIL — el mínimo MPI no completa**

Este gate sólo cambió procedencia, aislamiento, portabilidad y documentación.
No cambió formulación FSI, parámetros físicos, tolerancias ni algoritmos. No se
modificó, limpió ni sobrescribió `~/OpenFOAM/ezamora-v2512` ni el árbol dirty
original `/home/ezamora/Workspace/solids4foam`.

## solids4foam fijado

Se eligió **commit base + parche**, pues no hay un fork autorizado con los
cambios requeridos:

```text
URL             https://github.com/solids4foam/solids4foam.git
commit          4b254fa5260e0ae94640d7404089bde73907fc2d
parche          external/solids4foam/patches/0001-g0-baseline-runtime.patch
SHA-256 parche  9099bbfcc247235ab9d736e5ed3e7587aceae5e1e4ad359e936db04a3ed207d0
SHA-256 diff    9099bbfcc247235ab9d736e5ed3e7587aceae5e1e4ad359e936db04a3ed207d0
overlays        ninguno
```

`external/solids4foam/prepare-checkout.sh` clona, hace checkout detached,
valida/aplica el parche y comprueba commit, hash del diff y ausencia de
untracked. Una reconstrucción adicional en `/tmp/g0er-s4f-verify` reprodujo
los hashes. Durante la auditoría se usó el repositorio original sólo como
mirror local de lectura; el build usó exclusivamente
`external/solids4foam/source`.

La auditoría original registró detached HEAD `4b254fa5`, el remoto anterior,
10 archivos modificados (482 inserciones, 35 eliminaciones) y sus untracked.
Se clasificaron como necesarios y se capturaron seis `.C/.H`:
`newMovingWallVelocityFvPatchVectorField`, `neoHookeanElastic` y
`nonLinGeomUpdatedLagSolid`. Cuatro campos `T` de `hotCylinder`, logs y
mallas/resultados de `sphericalCavity` son incidentales y se excluyeron. No
hubo archivos inciertos ni untracked requeridos. El árbol original conservó
el mismo estado al final.

## Versiones y activación

| Componente | Versión efectiva |
|---|---|
| SO | Ubuntu 22.04.5 LTS; kernel 6.8.0-136-generic x86_64 |
| OpenFOAM | OpenCFD v2512, `linux64GccDPInt32Opt` |
| solids4foam | `4b254fa5` + parche fijado |
| GCC/G++ | 11.4.0 |
| CMake | 3.22.1 |
| Open MPI | 4.1.2; `libmpi.so.40.30.2` |
| deal.II | 9.3.2 |
| PETSc OpenFOAM | 3.25.3, `arch-openmpi-cpu-opt-nox` |
| PETSc deal.II | 3.15.5 del sistema |
| Otras | HDF5/OpenMPI, HYPRE 3.1.0, METIS/ParMETIS, Trilinos |

Activación reproducible:

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
export PETSC_DIR=/path/to/petsc-3.25.3
export PETSC_ARCH=arch-openmpi-cpu-opt-nox
source env/activate-project.sh
```

El script deriva la raíz desde sí mismo, define `S4F_ROOT`,
`WM_PROJECT_USER_DIR`, `FOAM_USER_APPBIN` y `FOAM_USER_LIBBIN`, y usa
`platforms/g0e-846a5f0/$WM_OPTIONS/{bin,lib}`. Elimina de `PATH` y
`LD_LIBRARY_PATH` las plataformas `ezamora-v2512`, antepone prefijo y PETSc,
exige los cinco artefactos y falla ante `not found` o enlace legacy. No toca
`.bashrc`. Desde una shell nueva, `type -a`, `which` y `readlink -f` muestran
únicamente el `solids4Foam` del prefijo aislado.

## Build aislado

```bash
export PETSC_DIR=/home/ezamora/Workspace/petsc-3.25.3
export PETSC_ARCH=arch-openmpi-cpu-opt-nox
export BUILD_JOBS=4
scripts/build-isolated-g0e.sh
```

El script valida commit/diff; ejecuta `Allwmake -j 4`; `wclean/wmake libso`
para `robinRobinCoupling` y `fiveParameterMooneyRivlinElastic`; y CMake Release
para `dealiiPdmsSolid`. Todo se instala en el prefijo aislado. Exit 0;
609.41 s de pared, 2047.97 s usuario, 102.32 s sistema, RSS 1,528,280 KiB.
No hubo errores. Las 397 líneas con `warning:` corresponden principalmente a
392 repeticiones del aviso TBB `task.h` deprecado; el resto son deprecaciones
OpenFOAM y dos `maybe-uninitialized`.

| Artefacto | Bytes | SHA-256 |
|---|---:|---|
| `bin/solids4Foam` | 27,344 | `31c342343fc6185776398a888e235185ad1ca56bd31006f856b5b70ea09e2e41` |
| `bin/dealiiPdmsSolid` | 7,927,392 | `0513de624ad0dc0093c47a8aff09517bf6c0e1c5f40407308eae1dcd0ff5ae4d` |
| `lib/libsolids4FoamModels.so` | 19,015,720 | `ad3e0764c3bd66d798b5036667c1f1d453c7370d247dadf3cc2efcc684c7228f` |
| `lib/librobinRobinCoupling.so` | 723,216 | `1196c2f04a71910e10b366b8a5a6c460d7c4f4bc8fa0a925ddc04d00aa1f92df` |
| `lib/libfiveParameterMooneyRivlinElastic.so` | 1,301,280 | `e4f8a57348d53d26964346494097f029d21ec45932e1196ccc0c371b37afd13e` |

## Enlazado, PETSc y MPI

Con activación, `ldd` no muestra `not found` ni rutas `ezamora-v2512`.
`solids4Foam`, acoplador y ley cargan `libsolids4FoamModels.so` local. Las
bibliotecas OpenFOAM base vienen de `/usr/lib/openfoam/openfoam2512` y las de
sistema de `/usr/lib`, ambas permitidas. `readelf -d` indica que los artefactos
OpenFOAM no tienen RPATH/RUNPATH: la resolución la controla y verifica la
activación. `dealiiPdmsSolid` tiene sólo RUNPATH del paquete hacia
HDF5/OpenMPI, sin rutas personales.

Mapas `/proc/<pid>/maps` capturados durante el fuerte:

| Proceso | PETSc | MPI |
|---|---|---|
| `solids4Foam` | `libpetsc.so.3.25.3` del PETSc fijado | `libmpi.so.40.30.2` |
| `dealiiPdmsSolid` | `libpetsc_real.so.3.15.5` del sistema | `libmpi.so.40.30.2` |

Son procesos separados; cada uno carga una sola versión PETSc y ambos usan la
misma ABI MPI. No existe mezcla PETSc 3.25/3.15 en un espacio de proceso.

## Portabilidad y caso

Los `Make/options` reciben `S4F_ROOT`/`SOLIDS4FOAM_ROOT`. El caso usa rutas
relativas para mallas, CSV, estados y parámetros; `runCase.sh` obtiene ambos
ejecutables del entorno. `prepareCase.sh` regenera las mallas; `solid.msh`
reproducida dio SHA-256
`dace8367d5cab024d4136e3b1ffa350cdb479ca888568bba93cdda40a467dd9b`.

La búsqueda de `/home/`, `ezamora-v2512` y el nombre del repo dejó sólo una
referencia documental y el patrón activo que rechaza la plataforma legacy. No
hay rutas personales en configuraciones activas. `.gitignore` ya no ignora
genéricamente `cases/`: versiona sólo inputs de `semiImplicitBetaIQNILS` y
mantiene ignorados tiempos, `processor*`, postproceso, logs, mallas derivadas,
CSV, estados/restarts, VTK/VTU y binarios.

## Pruebas en copias independientes

No se usaron resultados antiguos. En el mínimo serial, `endTime` se cambió
sólo en `/tmp` a `1e-6`: exit 0, `End`, 231.21 s, 10 pasos y correctores
`1,1,1,1,1,1,1,3,4,4`. Último corrector: residuo geométrico
`7.807517485e-5`, `max|r_u|=1.541915375e-7`, defecto de potencia
`9.67608909e-7`, ratio `0.5088320737`; cero NaN/Inf.

Sobre 3,645 muestras, máximo/RMS: desplazamiento
`2.101379583e-13 / 6.297507644e-14 m`; velocidad
`9.892006932e-7 / 3.070778798e-7 m/s`; aceleración
`5.338689470 / 1.694320606 m/s2`; tracción
`199.1533633 / 23.12918385 Pa`.

El mínimo MPI declara 4 rangos: fluido `(2 2 1)`, sólido `(4 1 1)`. Un ensayo
válido llegó a `8e-7` y abortó a 100 correctores: residuo geométrico
`2.869738984e-5`, ratio físico `5.419727823`, exit MPI 1, 1,426.03 s y sin
NaN/Inf. Una repetición exclusiva, sin otros casos FSI, volvió a atascarse en
`8e-7`; se detuvo adicionalmente en 62 correctores al alcanzar residuo
`1.680353733e-4` y ratio `5.775932799`, porque el criterio ya había fallado.
Se descartó un intento que omitió exportar el ejecutable deal.II por ser error
del arnés, no evidencia numérica.

Dos copias fuertes reprodujeron exactamente `t=1.1e-6`, paso 11, 100
correctores y exit observado 134. Métricas idénticas: ratio `3.212720779`,
residuo geométrico `5.254843908e-5`, `max|r_u|=9.735517512e-7`, defecto de
potencia `1.319919402e-5`; cero NaN/Inf. `robin-in.csv` y `robin-out.csv`
fueron byte a byte idénticos, con SHA-256 `16282a82...3290` y
`9369af55...80a0`. Máximo/RMS del trial: desplazamiento
`4.256686703e-13 / 1.331954459e-13 m`, velocidad
`3.082259442e-6 / 9.791980312e-7 m/s`, aceleración
`32.51879265 / 9.957997670 m/s2`, tracción
`632.2268432 / 109.8604309 Pa`.

La comparación ordena por `(x,y,z)` y exige
`||a-b||inf <= atol + 5e-5 max(||a||inf,||b||inf)`, con `atol` de `1e-15 m`,
`1e-12 m/s`, `1e-7 m/s2` y `1e-6 Pa`. A/B tiene diferencia cero. Alcanzar
`1e-4 s` en el fuerte es criterio de G7, no de G0.

## Gate

| Criterio | Estado |
|---|---|
| URL/commit solids4foam fijados y parche reproducible | PASS |
| Árbol dirty original excluido | PASS |
| Ejecutable y bibliotecas desde prefijo aislado | PASS |
| Sin dependencia runtime `ezamora-v2512` | PASS |
| Un PETSc por proceso y MPI ABI compatible | PASS |
| Rutas activas portables | PASS |
| Mínimo serial | PASS |
| Mínimo MPI | **FAIL**: aborto a 100 correctores en `8e-7` |
| Fuerte reproducible, sin NaN/Inf | PASS |
| Inputs, documentación y scripts versionados | PASS |
| Árbol limpio y push | PASS tras el commit documental (verificado al cierre) |

Por el FAIL MPI, **G0-E-R no se aprueba**. Los mapas completos de ambos
participantes se capturaron en el fuerte serial; en el mínimo MPI se verificó
estáticamente el mismo enlazado, pero no se conservó un mapa de cada rango.
Esto no cambia el diagnóstico del gate y queda como limitación documental.

Riesgos: sensibilidad de convergencia al orden paralelo; PETSc 3.25.3 sigue
siendo una dependencia externa fijada por versión/ruta; avisos de APIs antiguas;
y el fuerte permanece intencionalmente inestable hasta G7. No se inicia G1.
