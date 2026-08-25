# G0-E — Aislamiento y procedencia del entorno OpenFOAM

Fecha: 2026-08-25 (America/Santiago)  
Commit fuente auditado: `846a5f0a3f6d3abc9ffc32e83d4c96dc4fe7eb69`  
HEAD documental durante la auditoría: `0194b5005b7a09356a3603745ff48f46f2f137ce`  
Resultado: **FAIL — entorno compartido; no ejecutar casos FSI**

Los fuentes relevantes del workspace son idénticos a `846a5f0` (`git diff
--exit-code 846a5f0 -- src ...` devolvió 0), pero el ejecutable principal
`solids4Foam` y `libsolids4FoamModels.so` siguen procediendo de
`~/OpenFOAM/ezamora-v2512`. El repositorio fuente externo solids4foam está
dirty. Por tanto no puede atribuirse una futura ejecución FSI sólo al nuevo
repositorio ni aprobarse G0-E.

No se ejecutó ningún caso FSI durante este subgate. No se borró ni sobrescribió
ningún artefacto de `ezamora-v2512`.

## Identidad y entorno inicial

```text
pwd                              /home/ezamora/OpenFOAM/2way_Stokes_NeoHookean
realpath .                       /home/ezamora/OpenFOAM/2way_Stokes_NeoHookean
git rev-parse --show-toplevel    /home/ezamora/OpenFOAM/2way_Stokes_NeoHookean
git rev-parse HEAD               0194b5005b7a09356a3603745ff48f46f2f137ce
```

La shell preparada normalmente por `openfoam2512 -c` produjo:

```text
WM_PROJECT_DIR=/usr/lib/openfoam/openfoam2512
WM_PROJECT_VERSION=v2512
WM_PROJECT_USER_DIR=/home/ezamora/OpenFOAM/ezamora-v2512
WM_OPTIONS=linux64GccDPInt32Opt
FOAM_USER_APPBIN=/home/ezamora/OpenFOAM/ezamora-v2512/platforms/linux64GccDPInt32Opt/bin
FOAM_USER_LIBBIN=/home/ezamora/OpenFOAM/ezamora-v2512/platforms/linux64GccDPInt32Opt/lib
```

`PATH` antepone el `bin` anterior a los binarios de OpenFOAM del sistema y
`LD_LIBRARY_PATH` antepone el `lib` anterior. En consecuencia, el build
canónico sin overrides:

```bash
openfoam2512 -c 'cd src/robinRobinCoupling && wmake libso'
```

instala realmente en
`/home/ezamora/OpenFOAM/ezamora-v2512/platforms/linux64GccDPInt32Opt/lib`, no
en el nuevo repositorio. El entorno por defecto queda clasificado como
**compartido**.

## Inspección de build y ejecución

- `src/robinRobinCoupling/Make/files` define
  `LIB=$(FOAM_USER_LIBBIN)/librobinRobinCoupling`.
- `src/robinRobinCoupling/Make/options` fija
  `S4F_ROOT=/home/ezamora/Workspace/solids4foam`, incluye sus cabeceras y busca
  `-lsolids4FoamModels` exclusivamente mediante `-L$(FOAM_USER_LIBBIN)`.
- La ley adicional usa el mismo patrón y genera
  `$(FOAM_USER_LIBBIN)/libfiveParameterMooneyRivlinElastic.so`.
- `src/dealiiPdmsSolid/CMakeLists.txt` no define `install()`: CMake deja
  `dealiiPdmsSolid` dentro del directorio de build elegido. El smoke script usa
  por defecto `/tmp/dealii-pdms-build`.
- `semiImplicitBetaIQNILS/runCase.sh` invoca `openfoam2512 -c ... solids4Foam`
  y sólo antepone PETSc; por ello vuelve a seleccionar el `solids4Foam` de
  `ezamora-v2512`.
- Sus `fsiProperties` apuntan a `/tmp/dealii-pdms-build/dealiiPdmsSolid` y sus
  CSV/estado/parámetros contienen rutas absolutas bajo
  `/home/ezamora/OpenFOAM/ezamora-v2512/cases/...`.
- Existen muchas referencias equivalentes en otros casos y scripts, además de
  `LD_PRELOAD`/`LD_LIBRARY_PATH` absolutos para PETSc 3.25.3.
- No se encontraron referencias de build útiles a `$WM_PROJECT_USER_DIR` que
  corrijan automáticamente estas rutas absolutas.

## Inventario previo

Metadatos antes del nuevo build:

| Ruta | Tamaño | mtime | SHA-256 |
|---|---:|---|---|
| `ezamora-v2512/.../bin/solids4Foam` | 27,344 | 2026-07-28 12:18:19 -0400 | `31c342343fc6185776398a888e235185ad1ca56bd31006f856b5b70ea09e2e41` |
| `ezamora-v2512/.../lib/libsolids4FoamModels.so` | 19,015,720 | 2026-08-12 13:24:04 -0400 | `3c65a7a5de5a14d29d98db4765bd69abaa9abc91201588d6e93e8daeba534b33` |
| `ezamora-v2512/.../lib/librobinRobinCoupling.so` | 723,216 | 2026-08-25 01:24:24 -0400 | `366dfdf52c3fd7752dcd83c0eea85e1ce0c445c5da8507784b4638fdf3b633bf` |
| `ezamora-v2512/.../lib/libfiveParameterMooneyRivlinElastic.so` | 1,301,280 | 2026-08-25 01:23:59 -0400 | `68fe05522591c92788c4e8e9b57ed3d1e38e2a8fbb90455b3fb98876c4187121` |
| nuevo repo `platforms/.../bin/solids4Foam` | 27,344 | 2026-08-25 01:09:47 -0400 | `31c342343fc6185776398a888e235185ad1ca56bd31006f856b5b70ea09e2e41` |
| nuevo repo `platforms/.../lib/libsolids4FoamModels.so` | 19,015,720 | 2026-08-25 01:09:47 -0400 | `3c65a7a5de5a14d29d98db4765bd69abaa9abc91201588d6e93e8daeba534b33` |
| nuevo repo `platforms/.../lib/librobinRobinCoupling.so` | 723,216 | 2026-08-25 01:09:47 -0400 | `366dfdf52c3fd7752dcd83c0eea85e1ce0c445c5da8507784b4638fdf3b633bf` |
| nuevo repo `platforms/.../lib/libfiveParameterMooneyRivlinElastic.so` | 1,301,280 | 2026-08-25 01:09:47 -0400 | `68fe05522591c92788c4e8e9b57ed3d1e38e2a8fbb90455b3fb98876c4187121` |
| `/tmp/g0-20260825/build-dealii/dealiiPdmsSolid` | 7,927,320 | 2026-08-25 01:32:35 -0400 | `0cde590160bc7fa9f769b5dc3d6cbf43bd1cad8fb593d65da2f440fcd12e03f5` |
| `/tmp/g0-dealii-pdms-build/dealiiPdmsSolid` | 7,927,320 | 2026-08-25 01:23:30 -0400 | `01a053f60d2f063915209c7dd311db69aa166f1d3a1ab4ff809a432657b8f298` |

Las cuatro copias del nuevo repo son bit a bit iguales a las de
`ezamora-v2512`; sin manifiesto de creación no prueban procedencia y se
clasifican como ambiguas. `type -a solids4Foam` y `which solids4Foam` en el
entorno inicial devolvieron únicamente
`/home/ezamora/OpenFOAM/ezamora-v2512/platforms/linux64GccDPInt32Opt/bin/solids4Foam`.
`dealiiPdmsSolid` no estaba en `PATH`.

## Ensayo de prefijo aislado

Se creó un prefijo nuevo, previamente inexistente e ignorado por Git:

```text
/home/ezamora/OpenFOAM/2way_Stokes_NeoHookean/platforms/g0e-846a5f0/
```

Los fuentes se extrajeron directamente con
`git archive 846a5f0 src` a `/tmp/g0e-846a5f0-src`; no se usó el working tree
para compilar. Desde `env -i`, una shell `--noprofile --norc` y
`source /usr/lib/openfoam/openfoam2512/etc/bashrc`, se definió:

```bash
export WM_PROJECT_USER_DIR=/home/ezamora/OpenFOAM/2way_Stokes_NeoHookean
export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/platforms/g0e-846a5f0/$WM_OPTIONS/bin
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/g0e-846a5f0/$WM_OPTIONS/lib
export PATH=$FOAM_USER_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:/home/ezamora/Workspace/petsc-3.25.3/arch-openmpi-cpu-opt-nox/lib:/home/ezamora/OpenFOAM/ezamora-v2512/platforms/$WM_OPTIONS/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=/home/ezamora/OpenFOAM/ezamora-v2512/platforms/$WM_OPTIONS/lib
```

`LIBRARY_PATH` es la alternativa mínima que permite que `wmake` encuentre la
dependencia externa pese a que `Make/options` sólo incluye
`-L$(FOAM_USER_LIBBIN)`. El prefijo propio queda primero en resolución; las
dependencias externas quedan explícitas después.

Builds efectuados, sin casos:

| Artefacto aislado | Exit | Duración | SHA-256 nuevo |
|---|---:|---:|---|
| `.../lib/librobinRobinCoupling.so` | 0 | 25.35 s | `366dfdf52c3fd7752dcd83c0eea85e1ce0c445c5da8507784b4638fdf3b633bf` |
| `.../lib/libfiveParameterMooneyRivlinElastic.so` | 0 | 20.07 s | `68fe05522591c92788c4e8e9b57ed3d1e38e2a8fbb90455b3fb98876c4187121` |
| `.../bin/dealiiPdmsSolid` | 0 | 16.68 s | `b5c607c9b23ae2bc9b251711584a6e19a114cb78ac10a7660f23afd0d085ab12` |

Los dos plugins tienen build determinista y repitieron el hash antiguo. Su
mtime/inodo y logs demuestran que se crearon en el prefijo nuevo, pero el hash
por sí solo no distingue las copias. `dealiiPdmsSolid` sí cambió respecto de
las dos copias `/tmp`. Los hashes y mtimes de `ezamora-v2512` permanecieron
inalterados después del ensayo.

## Resolución comprobada

En una segunda shell limpia con la activación anterior:

```text
which dealiiPdmsSolid
  -> .../2way_Stokes_NeoHookean/platforms/g0e-846a5f0/linux64GccDPInt32Opt/bin/dealiiPdmsSolid

which solids4Foam
  -> .../ezamora-v2512/platforms/linux64GccDPInt32Opt/bin/solids4Foam

ldd .../g0e-846a5f0/.../librobinRobinCoupling.so
  librobinRobinCoupling.so -> archivo inspeccionado del prefijo aislado
  libsolids4FoamModels.so  -> .../ezamora-v2512/.../lib/libsolids4FoamModels.so
  libpetsc.so.3.25         -> .../Workspace/petsc-3.25.3/.../libpetsc.so.3.25

ldd .../g0e-846a5f0/.../libfiveParameterMooneyRivlinElastic.so
  libsolids4FoamModels.so  -> .../ezamora-v2512/.../lib/libsolids4FoamModels.so

ldd .../g0e-846a5f0/.../bin/dealiiPdmsSolid
  deal.II/PETSc/Trilinos   -> bibliotecas del sistema
```

No hubo entradas `not found`. No se carga accidentalmente una copia antigua de
los dos plugins porque el prefijo local tiene precedencia. Sin embargo, el
ejecutable y el modelo base de solids4foam siguen siendo deliberadamente los
de `ezamora-v2512`; esto es precisamente el bloqueo de procedencia.

## Gate

| Criterio | Estado | Evidencia |
|---|---|---|
| Fuentes corresponden a `846a5f0` | PASS | extracción `git archive`; diff de fuentes vacío |
| Destino exacto de cada build conocido | PASS | prefijo `platforms/g0e-846a5f0` y build CMake explícito |
| Ejecutable usado proviene del nuevo build | **FAIL** | deal.II sí; `solids4Foam` resuelve a `ezamora-v2512` |
| Bibliotecas propias cargadas desde rutas esperadas | PASS | plugins aislados primero en `LD_LIBRARY_PATH` |
| Ninguna biblioteca antigua tiene precedencia accidental | PASS parcial | plugins no; solids4foam antiguo es dependencia explícita |
| Repetible desde shell nueva | PASS para builds propios | demostrado dos veces con `env -i`; no para solids4foam externo dirty |
| Rutas y comandos documentados | PASS | este informe |

Resultado G0-E: **FAIL**. G0 permanece no aprobado y no se debe ejecutar ningún
caso FSI.

## Conflictos, dependencia bloqueante y cambio mínimo

Rutas conflictivas/ambiguas:

- `~/OpenFOAM/ezamora-v2512/platforms/$WM_OPTIONS/bin/solids4Foam`;
- `~/OpenFOAM/ezamora-v2512/platforms/$WM_OPTIONS/lib/libsolids4FoamModels.so`;
- copias bit a bit iguales y sin manifiesto en el `platforms/$WM_OPTIONS` del
  nuevo repo;
- `/tmp/dealii-pdms-build` y `/tmp/g0-dealii-pdms-build`;
- rutas absolutas `ezamora-v2512` dentro de `fsiProperties` y `solid.prm`.

Dependencia bloqueante: el repositorio actual no contiene ni construye
`solids4Foam`/`libsolids4FoamModels.so`; `Make/options` consume el árbol externo
`/home/ezamora/Workspace/solids4foam`, que está dirty, mientras sus binarios
instalados viven en `ezamora-v2512`.

Cambio mínimo recomendado para un próximo subgate, no aplicado aquí:

1. fijar un commit limpio de solids4foam y registrar también sus parches, o
   crear un commit limpio que contenga los cambios requeridos;
2. compilar `solids4Foam`, `libsolids4FoamModels.so` y herramientas asociadas
   desde ese estado hacia el mismo prefijo aislado del nuevo repo;
3. hacer que `Make/options` acepte `S4F_ROOT` y un directorio de bibliotecas
   externo configurables, sin rutas personales ni dependencia de
   `FOAM_USER_LIBBIN` para encontrar solids4foam;
4. proporcionar una activación versionada que anteponga sólo el prefijo
   aislado y PETSc fijado;
5. sustituir rutas absolutas de casos por rutas relativas o generadas;
6. repetir `type -a`, `which`, SHA-256 y `ldd`; sólo entonces reabrir los casos
   FSI de G0.
