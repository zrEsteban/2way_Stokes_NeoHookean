# Five-parameter Mooney-Rivlin extension

`src/fiveParameterMooneyRivlinElastic` is a runtime-loaded solids4foam v2.3
mechanical law for OpenFOAM OpenCFD v2512. Its isochoric derivatives are

`W1 = C10 + 2*C20*(I1bar-3) + C11*(I2bar-3)` and
`W2 = C01 + C11*(I1bar-3) + 2*C02*(I2bar-3)`.

Following the installed `MooneyRivlinElastic`, it evaluates

`J*sigma = dev(2*W1*Bbar - 2*W2*inv(Bbar)) + 0.5*K*(J^2-1)*I`.

Thus the volumetric energy is the installed-law equivalent
`K/4*(J^2 - 1 - 2*ln(J))`, not `K/2*(J-1)^2`. This choice makes the
`C20=C11=C02=0` response exactly regress to the existing two-parameter law.
The solids4foam interface does not request a fourth-order consistent tangent;
it requests only `impK`, supplied as `K + 4*mu0/3` in the same manner as the
existing law. Both cell and boundary Jacobians are checked before fractional
powers are evaluated.

Coefficients must carry pressure dimensions. `coefficientScale` is explicit:
use `1e6` for numeric inputs in MPa and `1` for inputs in Pa. `bulkModulus` is
always supplied in Pa and is intentionally not affected by that scale.

Build and test:

```bash
openfoam2512 -c 'wmake libso src/fiveParameterMooneyRivlinElastic'
bash tests/fiveParameterMooneyRivlinElastic/Allrun
```

