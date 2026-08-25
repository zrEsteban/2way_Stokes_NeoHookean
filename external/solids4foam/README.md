# Fixed solids4foam dependency for G0

This directory implements the reproducible base-commit-plus-patch alternative.
It never reads or modifies the developer's dirty solids4foam worktree.

```bash
external/solids4foam/prepare-checkout.sh /path/to/new/solids4foam-checkout
```

The patch contains only six runtime source files required by this baseline:

- `newMovingWallVelocityFvPatchVectorField.C/.H`: ALE interface-velocity fix;
- `neoHookeanElastic.C/.H`: PDMS constitutive/restart diagnostics used by the
  solids4foam runtime library;
- `nonLinGeomUpdatedLagSolid.C/.H`: structural convergence and tangent support
  compiled into `libsolids4FoamModels.so`.

Dirty tutorial fields, logs and generated tutorial results are incidental and
are deliberately excluded. No untracked file from the original worktree is
required. The manifest fixes the upstream URL, commit, patch hash and expected
post-application diff hash.
