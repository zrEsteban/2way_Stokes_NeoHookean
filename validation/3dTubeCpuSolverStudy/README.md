# 3dTube CPU solid-solver study

This study compares interchangeable CPU linear solvers for the segregated
solid displacement equation. It does not modify existing solver libraries or
the original `run/solids4foam-validation/3dTube` case.

## Progressive gates

1. One time step, serial: compatibility and numerical equivalence.
2. Ten time steps, serial: repeated-setup cost and trajectory equivalence.
3. Full 3dTube transient: accepted configurations only.
4. MPI scaling: only after the serial full-transient gate.

GPU configurations are intentionally out of scope.

## Acceptance criteria

- Normal exit and solid momentum convergence at every time step.
- Same FSI iteration count as PCG/FDIC, or no more than 5% extra.
- Relative differences against PCG/FDIC below `1e-6` for probe displacement
  and integrated interface force.
- Report wall time, OpenFOAM time, maximum RSS, solid correctors, and HYPRE
  setup/reuse metrics where applicable.

## Stage 1: one step, 6400 solid cells

| Variant | Wall [s] | RSS [MiB] | FSI iterations | Result |
|---|---:|---:|---:|---|
| PCG + FDIC | 3.46 | 152.3 | 6 | reference/pass |
| native GAMG | 4.34 | 158.6 | 6 | pass, 25% slower |
| HYPRE PCG + BoomerAMG, no reuse | 14.28 | 453.4 | 6 | pass, stop here |
| HYPRE PCG + BoomerAMG, reuse | 6.57 | 459.8 | 6 | pass |

All variants produced identical reported probe displacement and integrated
forces. HYPRE reuse reduced AMG setups from 768 to 3 (765 reuses), and setup
time from 6.92815 s to 0.029537 s. It remains slower than PCG/FDIC on this small
serial case, so it is not yet a candidate for the PDMS production case.

The PETScSNES 3dTube setup is preserved. Validation copies use
`implicitSegregated` to exercise the OpenFOAM `lduMatrix` solvers.

## Stage 2: ten steps, 6400 solid cells

| Variant | Wall [s] | RSS [MiB] | Total FSI iterations | Result |
|---|---:|---:|---:|---|
| PCG + FDIC | 26.16 | 163.4 | 88 | reference/pass |
| native GAMG | 30.64 | 161.3 | 88 | pass, 17.1% slower |
| HYPRE PCG + BoomerAMG, reuse | 44.24 | 462.7 | 88 | pass, 69.1% slower |

The three complete trajectories agree at reported precision. At `t=0.00025`,
the displacement magnitude is `1.06837e-09` and the interface normal force is
`-0.0162691` for every variant. HYPRE performed 3 AMG setups and 6222 reuses;
setup cost was only 0.0302009 s, but cumulative Krylov solve time was 13.5654 s.
This confirms that reuse works correctly, while HYPRE is not competitive on
the current small serial mesh.

## Equivalent PDMS-scale gate

The enlarged 3dTube solid contains 491520 cells, only 0.94% fewer than the
496200-cell PDMS solid. Both were decomposed over four MPI ranks, giving about
123000 solid cells per rank (PDMS: 124050). The enlarged fluid region contains
1105920 cells. Both meshes pass `checkMesh`.

| Variant | Wall [s] | RSS [MiB] | FSI iterations | Numerical result |
|---|---:|---:|---:|---|
| PCG + FDIC | 264.87 | 1217 | 6 | reference/pass |
| HYPRE PCG + BoomerAMG, reuse | 401.28 | 1547 | 6 | equivalent, not competitive |

HYPRE was 51.5% slower and used 27.1% more peak resident memory. Probe
displacement, integrated fluid force, solid interface force and FSI history
agree with PCG/FDIC at reported precision. The HYPRE run completed the physical
solve and printed `The momentum equation converged in all time-steps`, but an
MPI rank returned status 1 during shutdown after `MPI_Finalize`; this lifecycle
defect must be fixed before the solver can be considered production-ready.

The current HYPRE profile also emits lifecycle and timing diagnostics for every
component solve. A quiet profiling run is the next optimization gate, followed
by isolated Krylov/AMG tuning only if shutdown is clean.

## CPU lifecycle and profiling implementation

An isolated library, `libhypreFoamSolversParallelSolidCPU.so`, was built from
the existing CPU implementation. The original HYPRE libraries and runtime
solver names were not overwritten. Its source is under `hypreCpuQuiet/`.

The implementation:

- detects OpenFOAM having already called `MPI_Finalize` and avoids late HYPRE
  destructors that would invoke the invalid `MPI_Comm_free`;
- separates per-solve diagnostics from the final accumulated summary;
- records inspection, context creation, matrix update/assembly, vector update,
  AMG setup, Krylov solve, solution copy and total wrapper times;
- records accumulated Krylov iterations for subsequent runs.

The corrected quiet PCG/BoomerAMG run finished with exit status 0. Its log fell
from 9712 to 533 lines. Numerical results and six FSI iterations remained
identical to PCG/FDIC.

| Quiet HYPRE metric (rank-master accumulation) | Time [s] |
|---|---:|
| Krylov solve | 216.084 |
| Matrix inspection | 15.0018 |
| Vector updates | 5.36533 |
| AMG setup | 1.08578 |
| Matrix assembly | 0.198684 |
| Matrix value updates | 0.0142485 |
| Context creation | 0.0275113 |
| Solution copy | 1.05422 |
| Total wrapper | 242.492 |

The quiet wall time was 400.61 s versus 264.87 s for PCG/FDIC. Therefore log
I/O and AMG setup are not the bottlenecks; Krylov dominates.

### Tuning gate

A separate `hypreFlexGMRESParallel` case tested BoomerAMG `relaxType 6`,
`strongThreshold 0.5` and `kDim 50`. It was rejected: after several hours it
had completed only one FSI iteration and was still solving the second solid
update. The run was terminated without altering the accepted PCG case. This
parameter set is explicitly not suitable for the 3dTube elasticity operator.

Current recommendation: keep PCG/FDIC for the PDMS production case. Retain the
corrected quiet HYPRE library as an experimental CPU profile, but do not promote
it until a future Krylov/preconditioner combination beats the native baseline.
