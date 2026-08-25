# deal.II PDMS solid participant

This executable is the structural participant for the file-based
OpenFOAM--deal.II Robin--Robin coupling. It solves three-dimensional dynamic,
compressible neo-Hookean elasticity with Q1 finite elements, Newton iterations,
and applies the solid Robin condition
weakly on boundary id 3:

```text
sigma_s n + Z_s v_s = robin_traction + Z_s target_velocity.
```

The backward-Euler interface velocity is `(u-u_old)/dt`. Consequently the
Robin impedance contributes `Z_s/dt` to the FEM tangent. This is deliberately
not a traction-only surrogate.

Input and output CSV files contain one row per interface sample:

```text
x,y,z,tx,ty,tz,vx,vy,vz
x,y,z,ux,uy,uz,vx,vy,vz,tx,ty,tz
```

The constitutive energy is the isochoric/volumetric form using the same `mu`
and `K` as the OpenFOAM case. The material tangent is evaluated by centered
finite differences of the analytical first Piola stress; this is slower than a
closed-form tangent but keeps the first implementation independently auditable.
Time integration is selected with `time integration = backwardEuler|bdf2`.
The BDF2 option applies BDF2 to both first-order equations `d_dot=v` and
`v_dot=a`; its first accepted step uses backward Euler.  State files written
by the new solver contain `(d^n,v^n,d^{n-1},v^{n-1},history_depth)`.  Legacy
two-vector states remain readable and deliberately take one BE restart step
before enabling BDF2.
