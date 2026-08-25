# G0-CLOSE — cierre formal de G0

Fecha: 2026-08-25 (`America/Santiago`)

## Decisión

**G0: PASS.** La procedencia, aislamiento, builds existentes y participantes
runtime ya aprobados por G0-E-R no se recompilaron. G0-CLOSE corrigió sólo el
procedimiento de descomposición, añadió su regresión y realizó dos ejecuciones
MPI nuevas. No se modificó C++, física, tolerancias, discretización ni
algoritmos FSI. No se repitió el baseline fuerte.

Estado Git de entrada: `19f62b9cd09c316887a3363a1cf52f16db016480`,
árbol limpio. Commit canónico de fuentes:
`846a5f0a3f6d3abc9ffc32e83d4c96dc4fe7eb69`.

## Procedimiento MPI versionado

El flujo obligatorio es:

```bash
export PETSC_DIR=/path/to/petsc-3.25.3
export PETSC_ARCH=arch-openmpi-cpu-opt-nox
source env/activate-project.sh

case_copy=/path/to/fresh/case-copy
"$case_copy/prepareCase.sh"
"$case_copy/decomposeCase.sh"
export DEALII_PDMS_SOLID_BIN="$FOAM_USER_APPBIN/dealiiPdmsSolid"
"$case_copy/runCase.sh" -parallel
```

`decomposeCase.sh`:

1. aborta ante cualquier `processor*` preexistente;
2. ejecuta `decomposePar -no-libs -region fluid` y
   `decomposePar -no-libs -region solid`;
3. ejecuta el validador antes de devolver éxito.

`runCase.sh -parallel` vuelve a ejecutar el validador justo antes de lanzar
`mpirun`; el número de rangos se obtiene de
`system/fluid/decomposeParDict`, salvo override explícito
`FSI_MPI_PROCS`. La activación exporta `TWO_WAY_FSI_ROOT`, lo que permite usar
copias fuera del repositorio sin rutas personales.

El validador enumera mediante `foamDictionary` los patches de tipo
`pdmsElasticWallPressure`. En cada `processor*/0/fluid/p` exige:

- la misma entrada `constantHs` que en `0/fluid/p`;
- igualdad del valor leído a precisión 17;
- dimensiones de campo idénticas;
- una ocurrencia por patch físico y procesador;
- fallo inmediato ante ausencia, alteración o campo incompleto.

## Regresión automática

Comandos ejecutados:

```bash
bash -n env/activate-project.sh \
  scripts/validate-decomposition-constantHs.sh \
  scripts/test-decomposition-preserves-constantHs.sh \
  cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/decomposeCase.sh \
  cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/runCase.sh

scripts/test-decomposition-preserves-constantHs.sh
```

Resultado: **PASS**. La prueba crea su propio `mktemp`, copia únicamente los
inputs rastreados del caso, regenera mallas, usa `decomposeCase.sh`, valida 4
entradas y dimensiones `[0 2 -2 0 0 0 0]`, comprueba el rechazo de
`processor*` preexistentes, elimina `constantHs` y altera su valor para
demostrar dos fallos esperados, restaura y valida, y elimina sólo su temporal.

## Revalidación MPI dirigida

Se crearon dos copias independientes desde `git archive HEAD`; se copiaron en
ellas los scripts candidatos de este cierre, se fijó sólo `endTime=1e-6` en
`/tmp`, se regeneraron mallas y se ejecutó exclusivamente
`decomposeCase.sh` seguido de `runCase.sh -parallel`. Build aislado existente,
4 rangos, `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`.

| resultado | MPI A | MPI B |
|---|---:|---:|
| exit | 0 | 0 |
| condición | End | End |
| pasos | 10 | 10 |
| correctores | `1,1,1,1,1,1,1,3,4,4` | idénticos |
| ratio físico final | 0.5088212164 | 0.5088212164 |
| residuo geométrico final | 7.807540571e-5 | 7.807540571e-5 |
| `max|desplazamiento|` [m] | 2.10137287524e-13 | idéntico |
| `max|velocidad|` [m/s] | 9.89194494027e-7 | idéntico |
| `max|aceleración|` [m/s2] | 5.33864319568 | idéntico |
| `max|tracción|` [Pa] | 199.153498906 | idéntico |
| NaN/Inf | ninguno | ninguno |

SHA-256, idénticos A/B:

```text
456fcad9861aa7d41f369967798ef1fcc113267b9c81b7868c32e6fae454c472  fsiResiduals.dat
9fe82c68db1d2e2b22b42d3efe01ae10c49376f06f6be082db4c5134bb14793a  robin-out.csv
```

Ambos archivos comparan `cmp=0` entre repeticiones.

Warning no numérico en ambos lanzamientos: Open MPI imprimió dos veces
`Authorization required, but no authorization protocol specified` al sondear
el entorno gráfico. No hubo GUI, fallo de biblioteca ni cambio de exit; ambos
solvers continuaron hasta `End`.

## Equivalencia serial/MPI

Se usó el control serial aislado vigente de G0-MPI-D, ejecutado con el mismo
build y sin cambios de código. SHA-256 serial:

```text
e2d3baab1e2844ffe74ef3611c7c66d563dc47aaa9127c594aaa9cf59204a7ee  fsiResiduals.dat
2ad4330e6c0c71d4572ee2b2783e11dce6a7380f8a9e13984c27a98fa8cbf527  robin-out.csv
```

Comparación por coordenadas, `rtol=5e-5`:

| cantidad | diferencia máxima absoluta | diferencia máxima relativa | estado |
|---|---:|---:|---|
| desplazamiento [m] | 1.091690458e-18 | 5.313541918e-6 | PASS |
| velocidad [m/s] | 7.340644365e-12 | 7.734522463e-6 | PASS |
| aceleración [m/s2] | 5.881275236e-5 | 1.172838648e-5 | PASS |
| tracción [Pa] | 3.856206050e-3 | 1.936802983e-5 | PASS |

Serial y MPI tienen los mismos 10 pasos y correctores. No se exige igualdad
binaria entre configuraciones MPI/serial.

## Defecto diferido obligatoriamente a G1

`pdmsElasticWallPressureFvPatchScalarField::write()` no serializa
`constantHs`. No se corrigió en G0-CLOSE. Hasta aprobar G1,
`decomposePar -no-libs` y la validación post-descomposición son obligatorios.

Criterios añadidos a G1:

- `write()` preserva todos los coeficientes físicos y estados necesarios;
- escribir y releer el campo es idempotente;
- `decomposePar` con el plugin cargado no cambia `constantHs`;
- restart serial y MPI conservan exactamente el estado Robin.

## Baseline fuerte

No se repitió: no cambió código ni binario. Sus dos ejecuciones aisladas
permanecen como baseline negativo reproducible: fallo en `t=1.1e-6`, paso 11,
100 correctores, ratio `3.212720779`, residuo geométrico `5.254843908e-5`, sin
NaN/Inf. Hacerlo converger hasta `1e-4` pertenece a G7.

## Tabla formal

| criterio | estado | evidencia |
|---|---|---|
| Control de versiones y remoto | PASS | HEAD/origin verificados al publicar |
| Entorno y procedencia | PASS | G0-E-R; checkout y prefijo fijados |
| Build aislado | PASS | G0-E-R; no recompilado |
| Mínimo serial | PASS | End, 10 pasos, exit 0 |
| Procedimiento MPI seguro | PASS | `-no-libs`, abort clean-state, validador doble |
| Regresión de `constantHs` | PASS | positivo y negativos pasan |
| Dos mínimos MPI | PASS | End/exit 0; resultados byte a byte |
| Equivalencia serial/MPI | PASS | todas las cantidades bajo `rtol=5e-5` |
| Sin NaN/Inf | PASS | serial y ambas MPI |
| Baseline fuerte reproducible | PASS | fallo conocido preservado; G7 |
| Defecto `write()` contenido | PASS | protección obligatoria; deuda G1 |
| Documentación y scripts versionados | PASS | este cierre |

## Archivos del cierre

- `env/activate-project.sh`
- `cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/decomposeCase.sh`
- `cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/runCase.sh`
- `scripts/validate-decomposition-constantHs.sh`
- `scripts/test-decomposition-preserves-constantHs.sh`
- `docs/G0-E-openfoam-environment-isolation.md`
- `docs/G0-MPI-discrepancy-diagnosis.md`
- `docs/G0-formal-closure.md`

Riesgos abiertos: la protección procedimental sigue siendo necesaria hasta
G1; PETSc 3.25.3 continúa como dependencia externa fijada; la inestabilidad
fuerte continúa asignada a G7.

**Decisión final: G0 PASS. Detenerse; no iniciar G1.**
