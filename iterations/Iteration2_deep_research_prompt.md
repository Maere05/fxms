# Deep Research Prompt: Ahmed Drag Accuracy After SVBB Iteration 2

You are a deep research agent assisting with the next implementation iteration for a FluidX3D wall-function and force-validation project.

The implementation agent already has local access to this repository and the following files:

- `iterations/Iteration1.md`
- `iterations/Iteration2_summary.md`
- `context/deep research.md`
- `context/deeo research 2.md`
- `src/defines.hpp`
- `src/kernel.cpp`
- `src/lbm.cpp`
- `src/setup.cpp`

Your task is to research why the current Ahmed-body validation still reports unrealistically high drag and to produce a concrete implementation roadmap for the next coding iteration. The goal is to get more accurate overall drag coefficients for the 25 degree Ahmed body without ridiculous VRAM. The next validation campaign can use up to 48 GB VRAM on a server, so include resolution studies up to that budget.

## Current Project State

FluidX3D is being used for high-Re external aerodynamics with D3Q19 SRT, FP16C storage, Smagorinsky LES, Cartesian voxelized solid geometry, and force extraction through the Momentum Exchange Method.

The Ahmed validation case in `src/setup.cpp` uses:

- 25 degree Ahmed STL: `stl/ahmed_25deg_m.stl`
- `U = 40 m/s`
- `rho = 1.225 kg/m3`
- `nu = 1.48e-5 m2/s`
- body length `L = 1.044 m`
- width `W = 0.389 m`
- height `H = 0.288 m`
- current reference area `A_ref = W * H = 0.112032 m2`
- `Re ~= 2.82e6`

Canonical expectation for the 25 degree Ahmed body is about:

- `Cd ~= 0.285-0.30`
- drag force at this setup: `Fd ~= 31-33 N`

Current measured values are much higher:

```text
12G no_wall   meanCd 0.762002950, last10Cd 0.733525047, finalCd 0.699096560
12G cap32     meanCd 0.761628968, last10Cd 0.733142508, finalCd 0.684071280
14G no_wall   meanCd 0.788778098, last10Cd 0.802575516, finalCd 0.802109120
14G cap32     meanCd 0.788249657, last10Cd 0.801850206, finalCd 0.792223760
```

Iteration 2 changed SVBB wall-law viscosity handling. The final default is:

```cpp
#define WALL_MODEL_SVBB
#define WALL_MODEL_SVBB_NU_CAP 32.0f
#define WALL_MODEL_SVBB_DELTA_SIGN -1.0f
//#define WALL_MODEL_SVBB_FORCE_CORRECTION
//#define WALL_MODEL_SVBB_OBJECT_ONLY
//#define WALL_MODEL_SVBB_FLOOR_ONLY
//#define WALL_MODEL_DIAGNOSTICS
```

This made SVBB physically less degenerate and produced a consistent but tiny Cd improvement. It did not address the much larger error between `Cd ~= 0.76-0.79` and the expected `Cd ~= 0.285-0.30`.

## Core Research Question

Why does the Ahmed-body drag remain about 2.6x too high, and what should the next implementation agent change or measure to reduce total drag toward the experimental range under a practical 48 GB VRAM limit?

Do not focus only on wall shear. The current evidence says the main error is probably pressure and wake dominated. Research and rank the most likely causes, then convert them into implementable code and run-plan recommendations.

## Likely Causes To Investigate

### 1. Voxelized geometry causes artificial pressure drag

FluidX3D currently voxelizes the STL onto a Cartesian lattice. The 25 degree rear slant, rounded front, roof transitions, and any thin supports/features become stair-stepped solids. At high Reynolds number, these stair steps can act as geometric trips and separation triggers. If the flow separates prematurely at or near the slant because the slant is represented by 90 degree steps, the wake pressure drag will dominate. A wall function cannot repair this.

Research what geometric fidelity is needed for a 25 degree Ahmed body:

- how sensitive Cd is to slant-angle surface smoothness,
- how Cartesian cut-cell, interpolated bounce-back, Bouzidi-Firdaouss-Lallemand, multi-reflection, immersed boundary, or signed-distance treatments improve curved/slanted-wall drag in LBM,
- whether FluidX3D already has partial support for sub-grid wall distance, exact triangle intersection, or voxelization metadata that can be reused,
- whether storing a compact per-link wall distance or normal field is feasible under 48 GB.

Implementation output needed:

- a feasible sub-grid boundary treatment compatible with FluidX3D's local OpenCL kernels,
- memory cost estimate per cell or per boundary link,
- expected code touch points,
- fallback if a full interpolated boundary is too invasive.

### 2. Pressure drag and wake topology are not decomposed

The current force CSV only gives integrated force. We need to know where the drag comes from: front, roof, slant, base, underbody, floor-induced effects, supports if present, and skin-friction-like SVBB terms.

Research how to decompose LBM/MEM forces on voxelized objects:

- by local surface normal or link direction,
- by object region labels,
- pressure-like vs shear-like momentum exchange,
- per-link contribution histograms,
- wake pressure proxy from density field,
- separated wake diagnostics relevant to Ahmed body validation.

Implementation output needed:

- a region-labeling strategy for the Ahmed body in FluidX3D,
- a low-overhead force diagnostic that can output front/slant/base/roof/underbody drag contributions,
- a way to compare against reference contribution breakdowns, such as front pressure, slant pressure, base pressure, and skin friction.

### 3. Reference geometry and reference area may not match canonical results

The local case currently uses `A_ref = width * height = 0.112032 m2`. Some references include stilts or use projected area near `0.115 m2`; COMSOL reports a projected area of `0.059 m2` for a half-domain symmetry setup including stilts, which implies about `0.118 m2` full-domain. This difference is small compared with the current factor-2.6 Cd error, but it still matters for validation.

Research:

- canonical 25 degree Ahmed body dimensions, slant geometry, supports/stilts, ground clearance, and reference area,
- whether the local STL matches the canonical body, includes or omits stilts, and has the correct height after scaling,
- whether the current `height = 0.288 m` conflicts with any STL that includes supports or a 0.338 m total height convention,
- how the expected Cd changes depending on body-only vs body-plus-stilts and reference area.

Implementation output needed:

- an exact validation checklist for geometry scale, orientation, ground clearance, and reference area,
- a recommendation for whether to change `A_ref`,
- scripts or code hooks to print mesh extents, voxelized extents, projected area, and occupied boundary-link counts at runtime.

### 4. The domain and boundary conditions may be biasing pressure drag

The current Ahmed force-validation domain is:

- `x = 6 * width = 2.334 m`
- `y = 6 * length = 6.264 m`
- `z = 0.5 * (6 - 1) * width + height = 1.2605 m`

The mesh placement puts the body around the center, and the floor is solid. Need determine if inlet distance, outlet distance, height, side clearance, and top/side boundary conditions are adequate for Ahmed validation.

Research:

- standard Ahmed tunnel domain sizes and blockage recommendations,
- effect of side/top slip vs no-slip vs periodic/open boundaries in FluidX3D,
- required outlet length for wake pressure recovery,
- whether the current 6L streamwise domain is too short,
- whether the top and side boundaries in this setup constrain the wake or pressure field.

Implementation output needed:

- recommended 12 GB, 24 GB, and 48 GB domain/resolution matrices,
- whether to prioritize larger domain or finer body resolution,
- changes to `src/setup.cpp` for server runs.

### 5. Force averaging is too short for wake-dominated drag

The current force-validation pipeline uses:

- `init_steps = 2000`
- `sample_count = 20`
- `sample_interval = 50`

At the reported unit conversion, this is a short physical averaging window for a bluff-body wake. The force traces are noisy and final values jump strongly. Small wall-function deltas are not meaningful until the wake average converges.

Research:

- appropriate non-dimensional averaging time for Ahmed body simulations,
- Strouhal number and wake time scales for Ahmed body at 25 degrees,
- how many convective times `L/U` are needed before sampling and during sampling,
- how to define convergence criteria for Cd mean and standard error.

Implementation output needed:

- revised `init_steps`, `sample_interval`, and `sample_count` for 12 GB, 24 GB, and 48 GB runs,
- a CSV post-processing script that reports mean, last-N mean, standard deviation, standard error, and convergence plots,
- criteria for stopping a run once Cd uncertainty is acceptable.

### 6. Wall model is useful but not sufficient

SVBB mainly affects skin friction. Reference Ahmed drag decompositions indicate pressure drag dominates. COMSOL's example lists a measured total `Cd = 0.285` with skin friction about `0.055`. Even eliminating all skin friction would not close the current `Cd ~= 0.76` gap.

Research:

- expected skin-friction fraction for 25 degree Ahmed body at this Reynolds number,
- whether wall functions can alter slant separation enough to matter,
- whether object-only or floor-only wall modeling should be used,
- whether a floor wall function, no-slip floor, moving ground, or slip floor better matches the target experiment.

Implementation output needed:

- a wall-model test matrix that separates object and floor effects,
- a recommendation on whether to keep cap 32, replace Werner-Wengle, add van Driest damping, use log-law inversion, or make the wall model dependent on reconstructed normals,
- diagnostics that report estimated skin-friction contribution separately from pressure-like force.

### 7. Collision model and lattice choice may affect bluff-body pressure drag

The current configuration is D3Q19 SRT with Smagorinsky. This may be too dissipative or anisotropic for accurate separated external aerodynamics on a voxel grid.

Research:

- evidence for D3Q27 vs D3Q19 in external aerodynamic drag prediction,
- MRT/TRT/cumulant/entropic LBM advantages for bluff-body drag and curved boundaries,
- whether FluidX3D supports any relevant alternatives already,
- expected memory and VRAM impact under 48 GB.

Implementation output needed:

- whether switching to D3Q27 is worth a 48 GB run,
- whether SRT limitations are likely material here,
- practical recommendations constrained by existing FluidX3D architecture.

## Required Research Deliverable

Produce a report aimed at an implementation agent. It must include:

1. Ranked root causes for `Cd ~= 0.76-0.79` instead of `Cd ~= 0.285-0.30`.
2. Evidence for each root cause, including references where possible.
3. A concrete next-iteration implementation plan, split into:
   - low-risk instrumentation,
   - medium-risk setup changes,
   - high-risk boundary/geometry changes.
4. A 48 GB server run matrix with expected VRAM, grid size, domain size, and run duration concerns.
5. Exact metrics to collect:
   - Cd mean, standard deviation, standard error,
   - force by object region,
   - force by link direction,
   - pressure-like vs shear-like proxy,
   - wake diagnostics,
   - boundary-link counts and wall-model activation counts.
6. A decision tree:
   - if higher resolution alone drops Cd substantially, what to do next,
   - if higher resolution does not drop Cd, what boundary/geometry treatment to implement,
   - if force decomposition shows front/slant/base dominance, how to prioritize fixes,
   - if skin-friction contribution is already small, how to stop over-tuning SVBB.
7. Specific code areas likely to change in FluidX3D:
   - `src/setup.cpp` for validation case setup and run length,
   - `src/kernel.cpp` for boundary treatment and force diagnostics,
   - `src/lbm.cpp` for OpenCL define propagation and diagnostic buffers,
   - any mesh/voxelization code paths needed for projected area, wall normals, or wall distance.

## Important Constraints

- Do not propose a method that requires resolving `y+ < 1` over the whole body unless it fits within 48 GB and is part of a comparison run only.
- Prefer changes that preserve FluidX3D's GPU throughput and local-kernel structure.
- Any per-cell storage proposal must include VRAM cost.
- Any boundary treatment proposal must explain how it interacts with MEM force extraction.
- The implementation agent needs actionable formulas and pseudocode, not only high-level CFD advice.
- The final goal is total Ahmed Cd accuracy, not merely a more elegant wall function.

## Useful Reference Targets To Verify

The research agent should independently verify benchmark values rather than trusting this prompt blindly. Start with:

- Ahmed, Ramm, and Faltin, SAE 840300, 1984.
- Lienhart and Becker Ahmed-body validation data.
- COMSOL Ahmed body example: reports `U = 40 m/s`, `Re ~= 2.77e6`, measured total `Cd = 0.285`, skin friction `Cd = 0.055`.
- SimScale Ahmed validation: reports reference `Cd = 0.2875` and simulation values around `0.2835-0.2985`.
- TCFD Ahmed benchmark: reports measured `Cd = 0.2850` and CFD `Cd = 0.2848` at `40 m/s`.

## Expected Style

Be direct and implementation-oriented. Where there are competing explanations, rank them by how likely they are to explain a factor-2.6 drag overprediction. Separate facts, assumptions, and recommended experiments. Include enough mathematical detail for the coding agent to implement diagnostics or boundary changes without redoing the research.
