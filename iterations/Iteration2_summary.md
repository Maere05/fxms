# Iteration 2 Summary

This iteration started from the Iteration 1 SVBB wall-function implementation and focused on explaining why the Ahmed-body drag response was nearly inert.

## Code Changes

- Added wall-model tuning switches in `src/defines.hpp`:
  - `WALL_MODEL_SVBB_NU_CAP`
  - `WALL_MODEL_SVBB_DELTA_SIGN`
  - `WALL_MODEL_SVBB_FORCE_CORRECTION`
  - `WALL_MODEL_SVBB_OBJECT_ONLY`
  - `WALL_MODEL_SVBB_FLOOR_ONLY`
- Propagated the new switches into OpenCL in `src/lbm.cpp`.
- Changed the SVBB slip inversion in `src/kernel.cpp` to use capped wall-law viscosity instead of unconstrained local Smagorinsky effective viscosity.
- Added object-only and floor-only filter hooks for wall-model isolation tests.
- Added an optional moving-wall force correction path in `update_force_field()`, but left it disabled because it did not materially improve Ahmed drag.

Final default state after the sweep:

```cpp
#define WALL_MODEL_SVBB
#define WALL_MODEL_SVBB_NU_CAP 32.0f
#define WALL_MODEL_SVBB_DELTA_SIGN -1.0f
//#define WALL_MODEL_SVBB_FORCE_CORRECTION
//#define WALL_MODEL_SVBB_OBJECT_ONLY
//#define WALL_MODEL_SVBB_FLOOR_ONLY
//#define WALL_MODEL_DIAGNOSTICS
```

## Main Findings

The original effective-viscosity formulation produced excessive free-slip-like behavior but almost no change in integrated drag. A pure molecular-viscosity denominator collapsed all slip links to zero through reversal clamping. The best tested compromise was a capped viscosity denominator at `32 * nu_molecular`.

Ahmed 25 degree drag remained far too high compared with the expected `Cd ~= 0.285-0.30`, but cap 32 produced the first consistent improvement against same-resolution no-wall baselines.

## Run Results

All runs used the Ahmed 25 degree validation case at `U = 40 m/s`, `A_ref = 0.112032 m2`, `L_ref = 1.044 m`, and `Re ~= 2.82e6`.

```text
case          meanCd       last10Cd     finalCd      delta meanCd   delta last10Cd
12G no_wall   0.762002950  0.733525047  0.699096560  baseline       baseline
12G cap16     0.762425913  0.734067860  0.699727920  +0.000422963   +0.000542813
12G cap32     0.761628968  0.733142508  0.684071280  -0.000373982   -0.000382539
12G cap64     0.762038021  0.733567060  0.686141040  +0.000035071   +0.000042013
14G no_wall   0.788778098  0.802575516  0.802109120  baseline       baseline
14G cap32     0.788249657  0.801850206  0.792223760  -0.000528441   -0.000725310
```

Saved CSVs are in `bin/export/force_validation/`.

## Interpretation

The tiny Cd improvement is not surprising anymore. The current error is not mainly wall shear. The expected force at `Cd = 0.285` is about `31.3 N`; the current simulation reports about `83-87 N`. That excess is pressure and wake dominated.

The likely blockers for accurate Ahmed drag are:

- voxel stair-stepping of the 25 degree rear slant and rounded/front features,
- premature or excessive separation and wake pressure loss,
- coarse Cartesian geometry representation with no interpolated/curved boundary treatment,
- short force averaging for a wake-dominated transient,
- possible force decomposition blind spots: total MEM force is available, but pressure-vs-shear and surface-region contributions are not yet separated.

The next iteration should focus on pressure-drag accuracy and geometry/boundary fidelity, while using 48 GB server runs to determine whether resolution alone reduces the gap.
