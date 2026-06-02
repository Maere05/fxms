# Iteration 1: SVBB Wall Diagnostics

## Scope

This iteration added wall-model diagnostics around the SVBB wall-function path, validated the Ahmed body case at a 10 GB VRAM target, and compared drag with and without the SVBB wall functions enabled.

The main goal was to determine whether the wall model was inactive because the computed slip collapsed to zero, or whether it was active but not materially changing the measured force output.

## Repository Changes

### `src/defines.hpp`

- Added `WALL_MODEL_DIAGNOSTICS` as an optional compile-time flag.
- Left diagnostics disabled by default so the normal release build keeps the SVBB wall model active without the diagnostic overhead.
- Added `#undef WALL_MODEL_DIAGNOSTICS` under `BENCHMARK`, matching the existing pattern for feature flags that should not affect benchmark runs.

Current default wall-model state:

```cpp
#define WALL_MODEL_SVBB
#define WALL_MODEL_POSITIVITY_CLAMP
//#define WALL_MODEL_DIAGNOSTICS
```

### `src/lbm.hpp`

- Added the `Wall_Diagnostics` struct behind:

```cpp
#if defined(WALL_MODEL_SVBB) && defined(WALL_MODEL_DIAGNOSTICS)
```

- Added integer diagnostic counters for:
  - fluid cells checked
  - wall-adjacent cells
  - solid-neighbor links
  - links touching Ahmed object cells
  - links touching plain solid cells
  - zero tangential velocity links
  - nonzero slip links
  - reversal clamp links
  - positivity clamp links
  - nonzero population-delta links

- Added float accumulators for:
  - tangential velocity magnitude
  - wall shear stress
  - effective viscosity
  - raw slip
  - final slip
  - slip ratio
  - absolute population delta
  - signed population delta

- Added domain-level and top-level APIs for resetting and reading diagnostics.

### `src/lbm.cpp`

- Allocated diagnostic buffers on the device when `WALL_MODEL_DIAGNOSTICS` is enabled.
- Added those buffers to the `stream_collide` kernel argument list.
- Implemented:
  - `LBM_Domain::reset_wall_diagnostics()`
  - `LBM_Domain::add_wall_diagnostics(...)`
  - `LBM::reset_wall_diagnostics()`
  - `LBM::wall_diagnostics()`

- Added OpenCL preprocessor defines for diagnostic enum indices in `device_defines()`.

The integer counters are stored as low/high `uint` pairs to avoid overflowing 32-bit counters during large Ahmed runs.

### `src/kernel.cpp`

- Added OpenCL helper functions for diagnostic accumulation:
  - `wall_diag_add_u(...)`
  - `wall_diag_add_f(...)`

- Split the SVBB population and application functions into complete diagnostic and non-diagnostic variants.

This was necessary because embedding conditional macro string fragments inside unbalanced OpenCL function signatures caused a runtime OpenCL compile failure. The corrected implementation now emits complete function definitions for each mode.

- Instrumented the SVBB wall-link path to record:
  - how many wall-adjacent fluid cells are processed
  - how many solid-neighbor links are modified
  - whether links touch the Ahmed object or plain solid boundaries
  - whether tangential velocity is zero
  - whether slip is nonzero
  - whether raw slip is clamped to zero
  - whether the positivity clamp activates
  - whether the corrected population changes

- Kept the non-diagnostic wall model path functionally lightweight for normal builds.

### `src/setup.cpp`

- Added a diagnostics CSV path for force-validation runs:

```text
bin/export/force_validation/{case}_{memory}MB_wall_diag.csv
```

- Added a wall-diagnostics CSV header and per-sample append path.
- Limited wall-diagnostics CSV output to the Ahmed case.
- Reset diagnostics before each force sample interval and appended diagnostics immediately after the force sample.

## Build And Runtime Validation

The project was built with Visual Studio MSBuild using the Release x64 configuration:

```powershell
& 'C:\Program Files\Microsoft Visual Studio\18\Insiders\MSBuild\Current\Bin\MSBuild.exe' 'C:\dev\fxms\FluidX3D.sln' /p:Configuration=Release /p:Platform=x64 /m
```

The final repository state was rebuilt with:

- `WALL_MODEL_SVBB` enabled
- `WALL_MODEL_POSITIVITY_CLAMP` enabled
- `WALL_MODEL_DIAGNOSTICS` disabled

That leaves the production path with the wall model active and without diagnostic atomic overhead.

## Ahmed 10 GB Test Runs

Both validation runs used:

```powershell
& 'C:\dev\fxms\bin\FluidX3D.exe' ahmed 10000 0
```

Shared case details:

- Grid: `476 x 1278 x 257`
- Cells: `156,340,296`
- GPU memory target: about `9988 MB`
- GPU: NVIDIA GeForce RTX 5080
- Reynolds number printed by the run: `2,821,622`
- `nu_lbm`: `0.00000755`
- Ahmed object cells: `923,260`

Two runs were completed:

- Baseline without wall functions.
- SVBB diagnostics run with `WALL_MODEL_SVBB`, `WALL_MODEL_POSITIVITY_CLAMP`, and `WALL_MODEL_DIAGNOSTICS` enabled.

The diagnostic run was significantly slower because each wall-link sample used global atomics. This is expected and is why diagnostics remain disabled by default.

## Drag Comparison

| Metric | Baseline, no wall functions | SVBB diagnostics | Difference |
|---|---:|---:|---:|
| Mean drag, 20 samples | `83.428918738 N` | `83.433510782 N` | `+0.004592044 N` |
| Mean Cd, 20 samples | `0.75988603` | `0.7599278565` | `+0.0000418265` |
| Last-10 mean drag | `83.213882252 N` | `83.20369721 N` | `-0.010185042 N` |
| Last-10 mean Cd | `0.75792744` | `0.757834663` | `-0.000092777` |
| Final drag | `110.48485564 N` | `111.1957016 N` | `+0.71084596 N` |
| Final Cd | `1.00631654` | `1.01279103` | `+0.00647449` |

The mean and last-10-sample comparisons show that SVBB is effectively not changing the Ahmed drag result, even though the diagnostics show the wall model is active.

## Wall Diagnostics Findings

The Ahmed SVBB diagnostic run produced 20 diagnostic samples.

Average per 50-step sample:

| Diagnostic | Value |
|---|---:|
| Wall-adjacent cells | `33,636,600` |
| Solid-neighbor links | `167,371,500` |
| Links touching Ahmed object | `16,192,500` |
| Links touching plain solid | `151,179,000` |
| Object-link share | about `9.67%` |
| Plain-solid-link share | about `90.33%` |
| Slip nonzero links | `167,361,501.6` |
| Slip zero reversal clamp links | `9,998.35` |
| Positivity clamp links | `0` |
| Population delta nonzero links | `69,319,322.1` |
| Average final slip | `0.013761846` |
| Average slip ratio | `0.920814284` |
| Average effective viscosity | `0.001081683` |

Last diagnostic row at step `3000`:

| Diagnostic | Value |
|---|---:|
| Wall-adjacent cells | `33,636,600` |
| Solid-neighbor links | `167,371,500` |
| Links touching Ahmed object | `16,192,500` |
| Links touching plain solid | `151,179,000` |
| Tangential velocity zero links | `0` |
| Slip nonzero links | `167,356,215` |
| Slip zero reversal clamp links | `15,285` |
| Positivity clamp links | `0` |
| Population delta nonzero links | `65,784,012` |
| Average tangential velocity | `0.01534238` |
| Average wall shear stress | `0.00000159` |
| Average effective viscosity | `0.00108219` |
| Average raw slip | `0.01364819` |
| Average final slip | `0.01364855` |
| Average slip ratio | `0.91779496` |

## Main Conclusion

The original suspicion that SVBB is ineffective because the slip velocity collapses to zero is not supported by the diagnostics.

SVBB is active on almost all wall links:

- tangential velocity is nonzero
- final slip is nonzero
- reversal clamp activation is tiny
- positivity clamp activation is zero
- many populations receive nonzero corrections

Despite this, the Ahmed drag output is effectively unchanged. The likely issue is therefore downstream or structural: the correction is being applied but does not materially affect the force measurement, or it affects mostly boundaries that do not drive the reported Ahmed drag.

## Current Risks And Limitations

- The diagnostic float CSV uses fixed decimal formatting. Very small average population deltas may print as zero even when the nonzero-delta counter proves changes occurred.
- The diagnostics currently aggregate object links and plain solid links together for many averages. This makes it harder to isolate Ahmed-body behavior from floor/domain-wall behavior.
- The diagnostic atomics are intentionally expensive and should stay disabled for normal simulation work.
- The implementation proves SVBB activity, but it does not yet prove that the applied population correction is physically signed, indexed, or coupled to the force path correctly.

## Future Research

The next research pass should focus on why active SVBB population changes do not materially change the Ahmed force output.

Priority questions:

1. Verify whether the momentum-exchange force calculation reads the corrected post-SVBB populations or a state that predates the wall-function correction.
2. Audit the SVBB population correction sign and direction indexing, especially `link_i` versus `reflected_i`.
3. Split diagnostics into Ahmed-object links versus floor/domain-wall links for all counters and float averages.
4. Add direct accumulation of the SVBB momentum impulse and compare it to the force accumulated by the existing force-field path.
5. Re-run with object-only SVBB and floor-only SVBB variants to determine whether floor dominance is hiding object-wall effects.
6. Add scientific or higher-precision CSV formatting for population-delta diagnostics.
7. Validate SVBB on a simpler canonical wall-bounded case such as channel flow or Couette flow before using Ahmed drag as the primary indicator.
8. Compare pressure and viscous/shear drag contributions separately, if the existing force path can be decomposed or instrumented.
9. Check whether `nu_eff` from the relaxation-rate calculation is physically appropriate for the Werner-Wengle slip computation or whether molecular viscosity and modeled eddy viscosity should be handled differently.
10. Compare expected Ahmed Cd for this exact setup, domain, ground treatment, Reynolds number, and reference area against the measured values before treating the absolute drag as a validation target.

Recommended next diagnostics:

- Object-only wall diagnostic CSV columns for all averages.
- Per-link correction sign statistics grouped by lattice direction.
- A kernel-side SVBB impulse accumulator.
- A build-time experiment that flips the SVBB correction sign.
- A build-time experiment that swaps `link_i` and `reflected_i` in the correction term.
- A minimal channel-flow validation setup with expected wall-law behavior.

