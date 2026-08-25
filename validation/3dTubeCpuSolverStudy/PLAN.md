# Progressive work plan

## Phase A — completed

- Preserve existing libraries and the original 3dTube case.
- Isolate PCG/FDIC, GAMG, HYPRE/no-reuse and HYPRE/reuse variants.
- Validate one serial step with timing, memory, FSI, force and displacement.
- Stop HYPRE/no-reuse because AMG setup dominates.

## Phase B — completed

- Run ten serial steps for PCG/FDIC, GAMG and HYPRE/reuse.
- Compare full displacement, force and FSI-iteration time series.
- Separate HYPRE assembly, AMG setup and Krylov solve costs.

All three candidates were numerically equivalent. PCG/FDIC remains fastest;
GAMG was 17.1% slower and HYPRE/reuse was 69.1% slower in serial.

## Phase C — next

- Run the complete 3dTube transient for candidates passing Phase B.
- Repeat timings at least three times and report median and spread.
- Compare with PETScSNES when its CPU runtime is enabled.

## Phase D

- Test CPU MPI scaling at 1, 2 and 4 ranks on an enlarged 3dTube mesh.
- Promote a solver only after numerical equivalence and speedup are shown.

The first four-rank equivalent-size gate is complete. HYPRE/reuse is
numerically equivalent but 51.5% slower than PCG/FDIC and exits non-zero after
successful physics because of an MPI finalisation lifecycle defect. Before
additional scaling runs, disable per-solve diagnostics and fix the shutdown.

### Follow-up completed

- MPI shutdown fixed in an isolated CPU library; corrected run exits 0.
- Per-solve output disabled; accumulated profiling retained.
- Equivalent 3dTube repeated with identical physics and six FSI iterations.
- Krylov identified as dominant (216.084 s of 242.492 s wrapper time).
- FlexGMRES with stronger AMG smoothing evaluated and rejected after severe
  convergence degradation.
- PCG/FDIC remains the accepted production solver.

## Phase E

- Add the selected CPU solver as an optional PDMS `fvSolution` profile.
- Keep PCG/FDIC as fallback and leave existing implementations untouched.
