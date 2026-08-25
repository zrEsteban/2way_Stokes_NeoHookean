# Semi-implicit beta/AMP + IQN--ILS

Progressive case derived from the clean Dirichlet--Neumann reference. The
deal.II participant receives pure Neumann traction (`dealIIImpedance 0`), while
`pdmsElasticWallPressure` includes the diagonal interface inertia in the fluid
pressure solve. The current validated split is `beta=1`.

La masa local validada es la mitad de la profundidad de onda `ap*dt`:
`rhoSolidHs=0.003364938335 kg/m2`, equivalente a una impedancia inercial
`33649.38335 Pa s/m`. IQN--ILS acelera el desplazamiento, filtra modos casi
dependientes y reinicia la historia en cada paso; la reutilización entre pasos
fue ensayada y descartada por no superar el gate de contracción.

El gate serial aislado de diez pasos alcanzó `End`: los pasos 1--7 necesitaron
un corrector, el paso 8 tres, y los pasos 9--10 cuatro. La razón física final
fue `0.5088320737` y el residual geométrico relativo `7.807517485e-05`.

El caso está configurado hasta `endTime=1e-4 s`, usa el paso validado
`deltaT=1e-7 s` y arranca limpiamente desde cero. Se escribe cada `1e-6 s`,
aproximadamente cada diez pasos.
El mismo paso fijo se configura en deal.II y `adjustTimeStep` se desactiva. La
masa virtual escala con `dt` para conservar la impedancia AMP validada; este
salto temporal de cuatro órdenes de magnitud requiere un nuevo gate antes de
considerarse una configuración de producción.

Run with:

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
export PETSC_DIR=/path/to/petsc-3.25.3
export PETSC_ARCH=arch-openmpi-cpu-opt-nox
source env/activate-project.sh
./cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/prepareCase.sh
./cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/runCase.sh
```

Clean before every independent test with:

```bash
./cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/cleanCase.sh
```

All generated meshes, time directories, processor directories, CSV exchange
files, restart states and logs are ignored. The versioned inputs are sufficient
to regenerate the OpenFOAM meshes from `system/*/blockMeshDict` and the deal.II
mesh from `dealiiSolid/solid.geo`.
