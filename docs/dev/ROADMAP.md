# Development roadmap and status

Last updated: 2026-11-13 (feature/optimal-flow-components, in progress).

## Goal
Upgrade the steady-state heat-pump model, one feature at a time:
- accurate refrigerant properties without spline mesh errors;
- modular, composable components;
- refrigerant glide (zeotropic mixtures, ε-NTU with glide);
- a robust moving-boundary solution.

Keep what works well: mapping bounded unknowns to ℝ (sigmoid, including pressure) so that
unconstrained solvers can be used.

## Roadmap
| # | Branch | Content | Status |
|---|--------|---------|--------|
| B1 | `fix/sigmoid-inverse` | `sigmoid_01_to_R` was not the inverse of the sigmoid | merged (#18) |
| B3 | `fix/refrigerant-pcrit-sampling` | spline sampling up to `pcrit` (CoolProp error, cp ringing), below triple point; kwargs ignored | merged (#19) |
| B5 | `fix/borehole-exp-cache` | `SingleUPipeBorehole` not constructible with ExponentialUtilities 1.37 | merged (#20) |
| B2 | `fix/borehole-mean-temperature` | borehole coefficient at `(u7 + u9/2)` instead of the mean | merged (#21) |
| B6 | `fix/nusselt-continuity` | `SmoothTransition()` (default, VDI blend) / `SharpTransition()` (original) | merged (#22) |
| 1 | `feature/clapeyron-properties` | `ClapeyronRefrigerant` | merged (#23) |
| 2 | `feature/property-interface` | refrigerant interface + (T, v) state formulation, see below | merged (#27) |
| 3 | `feature/bounded-variables` | transforms as types (`LogBox` for pressure, `Box`, `Unit`); fix the absolute `pr_min` bound (B4) | merged (#28) |
| 4 | `feature/component-core` | lightweight components (bounded unknowns + scaled residuals + ports), ComponentArrays + NonlinearSolve | merged (#30) |
| 5 | `feature/glide-epsntu` | Qiao 2025 counter-flow ε-NTU with glide (eq. 56); pure fluid as the φ → 0 limit; Te := dew temperature | lower priority (natural refrigerants are pure); R125 blocker, see below |
| 6 | `feature/evaporator-operating-envelope` (was `robust-moving-boundary`) | Te lower bound (spurious zero-flow root), validity reports (brine freezing, issue #14), continuation helper; PTC evaluated | merged (#31) |
| 6c | `feature/optimal-flow-components` | port of the optimal-flow study to the component API; compressor/condenser on the property interface; validated against the original (needs PR #32) | PR open, see below |
| 6b | `feature/zone-switching` | zone appearance/collapse (Bell 2015 bounded-Q algorithm, Bonilla switching): needed for the condenser / R744 gas cooler or another expansion-valve control | deferred (cycle closure) |
| 7 | `spike/mtk` | how much ModelingToolkit integration is worthwhile | |
| 8 | `spike/dynamic-vs-steady` | is a dynamic model needed for robustness? coupling to a dynamic borehole model (warm starts, see below) | |
| – | `docs/dev-notes`, `ci/julia-versions` | this file and CLAUDE.md; CI on Julia `1`, `lts`, `min`, `pre` | merged (#24, #25) |

B4 (fixed in step 3): the pressure was mapped onto `(pr_min, pcrit)`, i.e. the reduced pressure
`pr_min` read as Pa (1e-2 Pa) instead of a pressure.

## Refrigerants of interest
- **Priority: natural, low-GWP pure fluids** (European move away from HFCs): R290 (main case),
  R744, R717, R600a, R1270. All have superancillaries (explicit saturation, sub-step 3).
- **Reference:** R134a.
- **Glide examples (step 5):** R407C (R32/R125/R134a, glide ≈ 6.4 K at 4 bar in CoolProp) and
  R410A (R32/R125, near-azeotropic); R454C as an extra low-GWP blend.
- Blocker for step 5: Clapeyron 0.6.30 (latest in the local registry) cannot load **R125**
  ("ResidualHelmholtzLemmon2005 not supported"), so neither R407C nor R410A can be built. Options:
  implement the Lemmon–Jacobsen 2005 term type in Clapeyron (upstream PR); CoolProp HEOS mixture
  through `AbstractState` as reference (no AD); `EmpiricPseudoPure("R410A")` has no glide.
- R744: evaporator subcritical (≈ 30 bar at −5 °C, pr ≈ 0.4), but the high side is usually
  transcritical: `compressor_power`/`condenser_power` assume condensation, a gas cooler needs its
  own branch. Check the Cooper correlation in `href_e` at R744's higher reduced pressure.

## Decisions
- Composition framework: lightweight custom layer, not ModelingToolkit (MTK is evaluated in step 7).
- Laminar-turbulent transition: keep both models; the smooth one is the default.
- Transport properties: CoolProp for pure fluids (reference values), EntropyScaling for mixtures.
- Keep the legacy code on the `legacy` branch.

## Key findings (with numbers)
- **Initial guess vs continuation** (R290 evaporator, 42 cases from 0.2 to 1.5 × design flow,
  0/5/10 °C source):
  - A fixed initial guess converges in only 24–30 of the 42 cases.
  - Continuation from the previous solution converges in all 42, in 3–5 iterations.
  - The corrected sigmoid guess speeds up only cases near the design point (5 instead of 162
    iterations).
- **Solver choice:** `TrustRegion`/`LevenbergMarquardt` stall on the coupled system because the
  residuals have very different scales (W ~ 1e3, dimensionless ε, K). `NewtonRaphson` from a
  nearby point converges in 3–5 iterations, down to −5 °C ground. Residual scaling belongs in
  step 4.
- **Laminar-turbulent transition:** the original sharp model makes the solution jump when the
  borehole Re crosses 2300. At Tb = 5 °C, going from 0.90 to 0.85 × design flow, Te drops by
  2 K and Q2 by 7%. The borehole runs at Re ≈ 2500–3200 at design flow.
- **Properties (Clapeyron, MultiFluid with CoolProp fluid files):**
  - Thermodynamics match CoolProp to ≤ 2e-9; AD works through saturation, dew/bubble and the
    p–s flash.
  - EntropyScaling vapour transport vs CoolProp: λ up to +17% (R290), μ up to 14%.
  - CoolProp cannot compute R32 saturated-vapour λ below ≈ −30 °C.
  - R454C (21.5/78.5 wt%): glide 7.9 K at 4.7 bar; agrees with CoolProp's HEOS mixture
    within 0.3 K.
- **Property speed (µs per call):**

  | Call | Time |
  |---|---|
  | `T_dew`, `h_dew` | 3.5 |
  | `h_pT`, `phase = :unknown` | 50 |
  | `h_pT`, `phase = :vapor` | 2.7 |
  | explicit `VT_enthalpy(v, T)` | 0.5 |
  | p–s flash (cold / with `T0`) | 37 / 18 |
  | R454C dew point (cold / warm `ChemPotDewTemperature`) | 118 / 72 |
  | `vol0` warm start for the volume solve | no gain |

## `feature/property-interface` (done; plan and results)
Idea: carry the thermodynamic state as **(T, v)** inside the model's nonlinear system. Each
property is then an explicit Helmholtz evaluation (`VT_*`, ≈ 0.5 µs, exact AD) with no inner
iterations. The outer Newton iteration acts as the warm start, and nothing is interpolated, so
there is no spline mesh error.

Base the branch on up-to-date `main` (all prerequisite fixes are merged).

1. **Interface (`AbstractRefrigerant`).** One set of functions used by the models:
   - saturation at p: `T_dew`, `T_bubble`, `h_dew`, `ρ_dew`, ...;
   - single phase at (p, T) with a phase hint;
   - transport.

   Implemented by the spline `Refrigerant` (wrapper around the existing interpolants) and by
   `ClapeyronRefrigerant`. Port `evaporator_model!`/`evaporatorborehole_model!` to the
   interface without changing equations or unknowns.

   Check: with the spline backend the residuals are bit-identical to `main`; with Clapeyron
   the solutions differ from the splines only by the spline error.

   **Done** (`src/properties/interface.jl`, `test/property_interface.jl`):
   - `AbstractRefrigerant` with `critical_pressure`, `molar_mass`, `pressure_bounds`, dew-point
     functions, `(p, T)` functions and `*_superheated(r, p, ΔT)` (vapour `ΔT` above the dew
     point; the spline backend interpolates in `(log10 p, ΔT)` directly, so it stays bit-exact).
   - `evaporator_model!`/`evaporatorborehole_model!` share `evaporator_residuals!`; the
     parameter structs are parametric in the refrigerant type (so `EvaporatorParameters(p;
     refrigerant = other)` cannot switch backend: build a new struct).
   - `pressure_bounds(::Refrigerant) = (pr_min, pcrit)` keeps B4 for bit-identity;
     `ClapeyronRefrigerant` uses `(0, pcrit)`.
   - Results (R290, 6 kW design, Ts_h = 278 K): spline residuals `==` the copy of main's; the
     Clapeyron solution matches the spline one to 1.6e-9 in pe, ≤ 1e-8 in the heat flows and ≤ 1e-7 K
     (also coupled to the borehole). Newton (ForwardDiff Jacobian) from the fixed design guess:
     6 iterations with either backend.
   - Cost per call: residual 8 µs (spline) vs 280 µs (Clapeyron + CoolProp transport); Jacobian
     9 µs vs 1.1 ms. Each Clapeyron residual repeats the dew-point solve 8 times and calls
     CoolProp for μ, λ with finite-difference derivatives: this is what sub-steps 2–3 target.
   - Not ported: `compressor_power`/`condenser_power` (still CoolProp by name, `::Refrigerant`).
2. **(T, v) state layer** for `ClapeyronRefrigerant`:
   - `state = (T, v)` with `p_Tv`, `h_Tv`, `s_Tv`, `cp_Tv`, `ρ_Tv`;
   - transport from (ρ, T). To verify: EntropyScaling density-based methods without Unitful;
     CoolProp via `PropsSI(..., "D", ρ, "T", T, ...)`;
   - the residual `p_Tv(T, v)/p - 1`, scaled from the start;
   - the molar volume as a bounded unknown on a log scale (from the co-volume to ∞), in line
     with the pressure mapping.
   **Done** (`src/properties/tv_state.jl`, `test/tv_state.jl`):
   - `p_Tv`, `h_Tv`, `s_Tv`, `cp_Tv`, `ρ_Tv`, `μ_Tv`, `λ_Tv` (v molar [m³/mol], outputs per kg),
     `pressure_residual`, `volume_pT` (for initial guesses), `volume_to_R`/`volume_from_R`
     (`v = vmin (1 + eᵘ)`, `vmin` = co-volume, a new field, independent of T for MultiFluid).
   - They reproduce the (p, T) functions to ≤ 1e-12 (vapour and liquid, R290 and R454C) and AD
     matches central differences.
   - Transport from (ρ, T): EntropyScaling has an unexported, unitless
     `VT_viscosity(model, V, T, z)` / `VT_thermal_conductivity` (the density method is
     Unitful-only). CoolProp: `PropsSI` costs ≈ 100 µs per call for any input pair, so
     `CoolPropTransport` now uses a cached low-level `AbstractState` (`CoolPropState`, with a
     lock), ≈ 2–3 µs, values `==` PropsSI. CoolProp has no analytic transport derivatives
     (`get_dT_drho[viscosity] is invalid`), so Duals still use central differences.
   - Cost per call: `p_Tv` 0.3 µs, `h_Tv` 0.45, `cp_Tv` 0.57, `μ_Tv` 0.6 (EntropyScaling) /
     1.8 (CoolProp); with Duals in (T, v): 0.5, 0.8, –, 1.0 / 11 µs. Compare `h_pT(:vapor)` 2.3,
     `T_dew` 2.7, `h_superheated` 6.4 µs.
   - Side effect on the existing models (CoolProp transport): evaporator residual 280 → 41 µs,
     Jacobian 1.1 ms → 115 µs (EntropyScaling transport: 69 / 163 µs).
   - Benchmark note: hand-made Dual tags (`Dual{Nothing}`) break Clapeyron's nested AD
     ("Cannot determine ordering of Dual tags"); use `ForwardDiff.Tag(f, Float64)`.

3. **Saturation without inner iterations** (pure fluids). Options:
   - Clapeyron superancillaries, if available for `MultiFluid`. To verify: `SuperAncSaturation`
     may only cover cubic/SAFT models;
   - otherwise (T, v_l, v_v) as extra unknowns with equal pressure and Gibbs energy, 2–3
     equations per saturation state.

   Mixture dew/bubble points (n + 2 unknowns each) go to the glide branch (step 5).
   **Done** (`saturation_T`, `p_sat` in `src/properties/tv_state.jl`; tests in `test/tv_state.jl`):
   - `SuperAncSaturation` covers only cubics (and PC-SAFT via EoSSuperancillaries), **but** when a
     CoolProp fluid file has a `SUPERANCILLARY` block Clapeyron loads it as the fluid's ancillary
     equations (`x0_psat`, `x0_sat_pure`). These are Bell's Chebyshev superancillaries: p, v_l,
     v_v explicit in T, phase equilibrium satisfied to ≤ 2e-11 (0.6–0.98 Tc), p within 1e-11 of
     CoolProp, exact AD (Clausius–Clapeyron to 7e-14). `NaN` outside (Ttriple, Tc).
   - Every CoolProp refrigerant tried has them (R290, R32, R1234yf, R134a, R600a, R1270, R744,
     R717, R1234ze(E)/(Z), R1233zd(E), R245fa, R152a, R143a, ...; R125 does not load in
     Clapeyron at all). `ClapeyronRefrigerant` checks them at construction (`superancillary`
     field); without them `saturation_T` falls back to `Clapeyron.saturation_pressure`.
   - Cost: `saturation_T` 0.2 µs (Dual 1.4 µs) vs 3.5 µs for the solver and 4 µs for `T_dew(p)`;
     dew-point state from T (`saturation_T` + `h_Tv`) 0.6 µs vs `h_dew(p)` 2.9 µs.
   - Consequence: for pure fluids the extra unknowns (T, v_l, v_v) with equal-p/equal-g residuals
     are not needed. Use the **saturation temperature** as the unknown instead of the pressure:
     p, v_l, v_v follow explicitly. The superancillaries are in T only, so keeping p as the unknown
     would still need an inner iteration (or a residual `p_sat(T)/p − 1`).

4. **Augmented evaporator**, behind an option. Add the superheat-outlet state (and optionally
   the saturation state) as unknowns with the residuals above. Compare with the nested
   formulation (properties solved inside the residual) on:
   - the 42-case sweep: success rate and iterations, from a fixed guess and with continuation;
   - cost per residual evaluation and per solve;
   - Jacobian conditioning.

   Report the results. Adopt the formulation in `component-core` only if it pays off.

   Benchmark fluids: R290, R744, R717 (and R134a as reference).
   Plan after sub-step 3: unknowns `Te` (sigmoid on (Ttriple, Tcrit), replacing `pme`) and the
   superheat-outlet volume `u_v = volume_to_R(v_sup)` with residual
   `pressure_residual(r, Te + ΔTsup, v_sup, p_sat(Te))`; transport at the mean superheat state from
   (T, v) as well (one more volume unknown at `Te + ΔTsup/2`, or interpolate the volume). Compare
   with the nested formulation (sub-step 1 port on `ClapeyronRefrigerant`).

   **Done** (`src/evaporator_formulations.jl`; tests in `test/property_interface.jl`; benchmark
   `examples/formulations_benchmark.jl`). `evaporator_model!(du, u, p, F)` / `evaporatorborehole_model!`
   take a `PropertyFormulation`: `NestedPressure()` (default, unknown pe, all backends),
   `SaturationTemperature()` (unknown Te, explicit saturation, 2 inner volume solves) and
   `ExplicitTv()` (Te + 2 volume unknowns, no inner iteration). Helpers: `evaporator_unknowns`
   (guess from physical values), `nested_unknowns` (for post-processing). The heat balances are
   shared, so all three give the same root (Te within 5e-11 K, Q within 1.4e-12).

   Sweep (6 kW design, 14 flows 0.2–1.5 × Vs, sources 0/5/10 °C; evaporator and evaporator +
   borehole; NonlinearSolve `NewtonRaphson`/`TrustRegion`, AutoForwardDiff; CoolProp transport):

   | | R290 | R744 | R717 | R134a |
   |---|---|---|---|---|
   | design guess, unshifted (R290 only) | 17/42 | | | |
   | guess = design state shifted to the source T, Newton | 32/42 | 34/42 | 30/42 | 30/42 |
   | continuation from it (up and down from the design flow) | 42/42 | 42/42 | 42/42 | 42/42 |
   | median Newton iterations (continuation) | 4 | 4 | 4 | 4 |
   | median time per solve, Nested / SatT / ExplicitTv [ms] | 3.5 / 1.3 / 1.3 | 7.1 / 2.9 / 3.4 | 4.2 / 1.2 / 1.2 | 3.6 / 1.3 / 1.3 |
   | residual [µs], Nested / SatT / ExplicitTv | 36 / 18 / 11 | 61 / 27 / 19 | 42 / 11 / 9 | 37 / 12 / 9 |
   | Jacobian [µs], Nested / SatT / ExplicitTv | 112 / 119 / 63 | 203 / 100 / 69 | 108 / 35 / 36 | 106 / 43 / 34 |
   | cond(row-scaled J), Nested / SatT / ExplicitTv | 1.3e4 / 3.2e4 / 3.5e4 | 7.3e3 / 7.3e3 / 2.1e4 | 1.3e4 / 2.6e4 / 3.0e4 | 1.2e4 / 2.8e4 / 3.1e4 |

   - **Convergence is the same for all three formulations, case by case** (also with the spline
     backend): failures are set by the model and the initial guess, not by the properties. The
     unshifted design guess puts Te above the secondary fluid for a 0 °C source. From the shifted
     guess only the low flows (≤ 0.4–0.5 × Vs) fail; continuation solves everything.
   - R744 runs with a subcritical stand-in for the high side (condensing at 25 °C for `hei`).
   - **Speed: ≈ 2.5–3× per solve** with Te as the unknown; `ExplicitTv` adds little over
     `SaturationTemperature` (Clapeyron's warm volume solve is cheap), at the price of two
     unknowns and a 2× worse conditioning than nested.
   - **R744** exposed two bugs (fixed in this branch): superancillary detection sampled below the
     triple point (0.6 Tc < Ttriple), and Clapeyron 0.6.30 `saturation_temperature` returns NaN for
     R744 at every pressure (upstream bug, to report). `T_dew(p)` for pure fluids now inverts the
     superancillary instead.
   - Recommendation for `component-core`: **use Te (dew temperature) as the unknown with explicit
     saturation (`SaturationTemperature`)**; keep `ExplicitTv` available but not as default. The
     robustness work belongs to step 6 (initial guesses from the source temperature,
     continuation in flow), plus residual scaling (step 4 of the roadmap).

Out of scope for this branch: residual scaling of the existing equations (step 4 of the
roadmap), glide (step 5) and new components.

## `feature/bounded-variables` (done)
1. **B4** (failing test, then fix): the spline `Refrigerant` stores its sampled pressure range
   (`pmin`, `pmax` fields) and `pressure_bounds` returns it (was `(pr_min, pcrit)` with `pr_min`
   read as Pa). Outside that range the splines extrapolate flat, so the Jacobian is singular there.
2. **Transforms as types** (`src/transforms.jl`): `VariableTransform` with `from_R(t, u)`,
   `to_R(t, x)`, `bounds(t)`; `Free()`, `Unit()` (0, 1), `Box(lo, hi)`, `LogBox(lo, hi)` (sigmoid
   in log x), `Above(lo, scale)` ((lo, ∞), `lo + scale·eᵘ`). The old `sigmoid_*`/`tanhmap_*`
   helpers stay exported.
3. **Models ported**: `pressure_transform(r) = LogBox(pressure_bounds(r)...)` (the pressure is now
   on a log scale), `temperature_transform(r) = Box(Ttriple, Tcrit)`, `volume_transform(r) =
   Above(vmin)`; `ϵe`, `ϵsup`, `xs` use `Unit()`. `evaporator_transforms(F, r; borehole)` lists the
   transforms of all unknowns (a first step towards the component unknowns of step 4), and
   `evaporator_unknowns` maps physical values with it. `ClapeyronRefrigerant` gets a `ptriple`
   field and `pressure_bounds = (ptriple, pcrit)` (LogBox needs `pmin > 0`; lowest component
   value for a mixture).
   **Breaking for scripts:** a pressure unknown built with `sigmoid_ab_to_R(pe, pressure_bounds(r)...)`
   is now wrong; use `to_R(pressure_transform(r), pe)` (examples and tests updated). To make old
   scripts fail instead of silently mis-mapping, `pressure_bounds` is no longer exported (read the
   interval with `bounds(pressure_transform(r))`).

Effect of the log-scale pressure on the sweep (`examples/formulations_benchmark.jl`, `sweep` +
`sweep2`, 4 fluids × evaporator / evaporator + borehole × Newton / TrustRegion, Clapeyron
backend; "cont" = 87 continuation solves):
- Same roots: |ΔTe| ≤ 4e-11 K wherever both converge.
- `NestedPressure` success rate unchanged within ±1 case almost everywhere; better with
  TrustRegion + continuation on the coupled system (R290 84 → 87/87, R717 82 → 87/87) and
  fixed guess (R290 22 → 24/42, median iterations 106 → 66); one regression (R290 evaporator,
  Newton, continuation, 73 → 71/87). Median iterations unchanged.
- `SaturationTemperature` identical case by case; `ExplicitTv` identical up to one iteration count
  (rounding of `vmin + vmin·eᵘ` vs `vmin(1 + eᵘ)`).
- Conclusion: the log map is neutral to slightly positive. Robustness is set by the initial guess and
  the residual scales, as found in step 2, and belongs in steps 4 and 6.

## `feature/component-core` (done)
Scope (decided 2026-11-13): **the evaporator first**, standalone (given secondary-fluid inlet
temperature) and coupled to the borehole. Compressor and condenser stay post-processing; the
design must allow closing the cycle later (e.g. `hei` becoming a port instead of a parameter).
The old API (`evaporator_model!`, `evaporatorborehole_model!`, `*Parameters`) stays unchanged;
tests check that the components give the same roots.

Design:
- **Component** (`AbstractComponent`): declares its internal unknowns (name → transform,
  physical values), its ports (name → transform), its residuals with a nominal scale each, and
  `residuals!(r, c, x, ports)` on **physical** values. The system does all the `from_R` mapping.
- **Ports as shared node variables:** a node joins ports of different components (e.g. the
  evaporator secondary outlet and the borehole inlet) into one unknown, or fixes them to a value
  (e.g. a given source temperature). No connection equations: the evaporator alone stays 8 × 8 and
  evaporator + borehole 9 × 9, as today.
- **System:** unknowns in a `ComponentVector` (`u.evaporator.Qe`, `u.nodes.T_cold`, ...), square
  check, `residual!(du, u, sys)`, `initial_guess(sys; ...)` from physical values and
  `physical(sys, u)` back. Generic in the number type (ForwardDiff).
- **Scaling:** unknowns through the transforms (new `Linear(x0, scale)`, e.g. heat flows in units
  of a nominal duty); residuals divided by their nominal scale (W → duty, K → 1, ε → 1). Newton is
  invariant to this, TrustRegion/LevenbergMarquardt are not: that is the measurement.
- **Formulation:** `SaturationTemperature()` by default for a pure `ClapeyronRefrigerant` with
  superancillaries (step 2 recommendation), `NestedPressure()` otherwise.
- **NonlinearSolve** as a package extension: `NonlinearProblem(sys, u0)` when it is loaded.

Implemented: `src/components/core.jl` (`AbstractComponent`, `ComponentSystem`, `residual!`,
`initial_guess`, `physical`, `port_values`), `src/components/evaporator.jl` (`Evaporator`,
`BoreholeHX`, `default_formulation`, `nominal_duty`, post-processing to `EvaporatorSolution` /
`BoreholeSolution`), `ext/OptimalGshpFlowNonlinearSolveExt.jl` (`NonlinearProblem(sys, u0)`),
`Linear(x0, scale)` in `src/transforms.jl`. The residual functions of the original models now
work on physical values (`refrigerant_state!`, `evaporator_balances!`), with the old functions
as wrappers (bit for bit). Tests: `test/components.jl`.

Checks (R290 design case):
- Evaporator alone: exactly the root of `evaporator_model!` (all three formulations, both
  backends); scaled residuals × scales == original residuals to 1e-12.
- Evaporator + borehole: root of `evaporatorborehole_model!` within 4.4e-6 K / 1.4e-7 relative
  in Qe (density, see below).
- Cost per residual / ForwardDiff Jacobian: the same as the original models (SatT 12 / 42 µs,
  ExplicitTv 10 / 38, Nested 38 / 86, spline 4 / 8).

**Scaling** (`examples/components_benchmark.jl`; default formulation `SaturationTemperature`,
CoolProp transport; 14 flows × 3 source temperatures; "unscaled" = `Qscale = 1`, i.e. the
original unknowns and residuals). Converged / 42 from the guess shifted to the source temperature:

| | R290 | R744 | R717 | R134a |
|---|---|---|---|---|
| Newton, unscaled → scaled (evap / evap+bh) | 32 → 32 / 32 → 32 | 34 → 34 / 34 → 34 | 30 → 30 / 30 → 30 | 30 → 30 / 30 → 30 |
| TrustRegion, unscaled → scaled | 29 → **42** / 29 → **42** | 36 → **42** / 33 → **42** | 26 → **42** / 26 → **41** | 28 → **41** / 29 → **42** |
| LevenbergMarquardt, unscaled → scaled | 10 → **42** / 6 → **42** | 6 → **42** / 3 → **42** | 8 → **42** / 8 → **42** | 8 → **42** / 10 → **41** |

- With continuation in flow everything converges with Newton and TrustRegion (scaled or not);
  LevenbergMarquardt needs the scaling (3–12/42 unscaled, 42/42 scaled).
- Newton is unchanged, as expected (affine invariance); its failures from the shifted guess are
  the low flows.
- Scaled TrustRegion: median iterations as Newton (4–8), maximum ≤ 48 (unscaled up to 196);
  scaled LevenbergMarquardt: median 12–28 iterations (unscaled 24–96).
- All variants reach the same root (|ΔTe| ≤ 6.5e-7 K).
- **Recommendation:** scaled system + `TrustRegion` as the default solve: it converges from a
  guess derived from the source temperature without continuation. Continuation remains the tool
  for zone appearance/collapse (step 6).

Known difference: the borehole component evaluates the fluid density at its own mean temperature
`(Tin + Tout)/2`, as `borehole_model!` does; `evaporatorborehole_model!` uses the evaporator
two-phase mean `(Ts_l + Ts_m)/2`, which a separate component cannot see. Roots differ by the
density change over ≈ 1 K (≈ 1e-4 relative).

## `feature/evaporator-operating-envelope` (done)
Planned as `robust-moving-boundary`; renamed after mapping the failures (2026-11-13, agreed with
the user). Map: R290 (then R744, R717, R134a), flows 0.05–5 × design, sources −5…15 °C, evaporator
alone and with the borehole, scaled `TrustRegion` from the shifted guess.

Findings:
- **No zone collapses in this evaporator model.** With fixed superheat ΔTsup and compressor volume
  flow, Te adjusts: at high flow `Te → Ts_in − ΔTsup` (superheat pinch), at low flow `xs → 0.98`.
  Zone switching (Bell et al. 2015: bounded solve for Q in (0, Q_max) from a pinch analysis;
  Bonilla et al.: switching moving-boundary models) matters only for other components or controls
  (condenser, R744 gas cooler, fixed orifice / flooded outlet): deferred to the cycle closure (6b).
- **Spurious zero-flow root.** As `Te → Ttriple` the vapour density at the compressor inlet, hence
  the refrigerant mass flow and every heat flow, vanish: a true root of the equations (Te = 86 K,
  Q ≈ 0, all brine temperatures = ground). TrustRegion reached it once (R290, borehole, −5 °C,
  0.3 × flow) because Te was only bounded by (Ttriple, Tcrit).
- **Brine freezing (issue #14).** At 0.05–0.2 × flow the brine leaves at −27 °C (MEA-28 % freezes at
  −18.2 °C) and Te reaches −28 °C; the brine splines covered −15…30 °C and extrapolated flat.

Implemented:
- `Evaporator(...; Te_min = default_Te_min(r))`: the first unknown is mapped on `(Te_min, Tcrit)`
  (or the pressure from the dew pressure at `Te_min`); default: dew temperature at pr = 0.01, or
  just above the triple point (R290 −60 °C, R744 −56.5 °C, R717 −31 °C, R134a −44 °C). The
  spurious root is out of reach: the case above now fails honestly from the fixed guess and
  converges with continuation (Te = −15.3 °C).
- `SecondaryFluid`: fields `Tmin`, `Tmax` (table range), `Tfreeze` (`freezing_point`, CoolProp
  `T_freeze`); `SecondaryFluid(name)` tabulates from just above the freezing point to 30 °C.
- **Validity reports (reported only, the solver is not constrained):** `validity(c, x, p)` per
  component, `validity_issues(sys, u)` for a system: brine below freezing or outside its table,
  Te at `Te_min`, two-phase and superheat pinch (0.1 K), a zone almost vanished (`xs` within 1e-3
  of 0 or 1).
- `continuation(build, λs, u0, alg; secant, max_halvings)` (NonlinearSolve extension): solves
  `build(λ)` along `λs` from the previous solution with a secant predictor, halving failed steps.

Results (480 cases = 4 fluids × 2 systems × 12 flows × 5 sources, scaled systems):
- Continuation in flow (TrustRegion, from the design flow up and down): **480/480**, median 4
  iterations per step. Reported issues: superheat pinch at high flow (4–22 cases per fluid and
  system), brine below freezing (3–12) or outside its table at low flow.
- Fixed (shifted) guess, no continuation: LevenbergMarquardt 477/480 (median 13–29 iterations),
  TrustRegion 462/480 (median 8–20; R717 is the hardest, 52–53/60), no wrong roots.
  `PseudoTransient()` does not help: 0/60 with the default initial step (and with α = 1), 56/60
  (R290) and 46/60 (R717) with α = 10, below TrustRegion and LevenbergMarquardt: not adopted.
- Recommendation: continuation for sweeps (what the optimal-flow study does anyway); for single
  points from a guess, LevenbergMarquardt (most robust) or TrustRegion (faster). Then check
  `validity_issues`.

## Warm starts for a dynamic coupling (measured 2026-11-13, for step 8)
Use case: the steady-state evaporator solved at every time step of a dynamic model (e.g. a dynamic
borehole whose wall temperature `Tb` changes slowly). Measured on R290, evaporator + borehole,
design flow, scaled system, after a change ΔTb of the ground temperature:

| ΔTb per step | cold start (TrustRegion, shifted guess) | warm start (Newton, previous `u`) | warm start + tangent predictor |
|---|---|---|---|
| 0.01 K | 6 it, 2.2 ms | 2 it, 0.95 ms | 1 it, 0.64 ms + 0.1 ms |
| 0.1 K | 6 it, 2.2 ms | 2 it, 0.95 ms | 1 it, 0.64 ms + 0.1 ms |
| 1 K | 5 it, 1.9 ms | 3 it, 1.3 ms | 2 it, 0.95 ms + 0.1 ms |

- Tangent predictor: `du/dTb = −J⁻¹ ∂F/∂Tb` at the previous solution (implicit function theorem;
  one ForwardDiff Jacobian and one derivative, 0.1 ms).
- One Newton iteration costs ≈ 60 µs (residual 12 µs + Jacobian 42 µs); the rest of the 0.64 ms
  is NonlinearSolve's per-call setup. Building the `ComponentSystem` costs 24 µs.
- Order of magnitude: hourly steps over 20 years ≈ 175 000 solves: ≈ 6 min cold, ≈ 3 min warm,
  < 30 s with a lean warm-started Newton. The evaporator is unlikely to be the bottleneck of a
  dynamic simulation (g-function convolutions of a dynamic borehole model probably are).

Decisions for step 8:
- No caching layer now: the warm start is "keep the last `u`" (what `continuation` does).
- Make the ground temperature `Tb` a **port** of `BoreholeHX` (connected to a node) instead of a
  field, so that a dynamic borehole component can provide it and a time step changes a node value
  instead of rebuilding the system.
- Only if profiling shows the solver overhead matters: reuse the solver cache (`init`/`reinit!`)
  or a small dedicated Newton loop (≈ 0.1–0.2 ms per step), and/or the tangent predictor.

## `feature/optimal-flow-components` (results 2026-11-13)
Goal (agreed 2026-11-13): port `examples/heatpumpborehole_optimalflowrate.jl` (the study the package
is named after: COP vs brine flow, evaporator + borehole, R290, ground −5 °C, ΔTsup 3.5 K) to the
component API, validate it against the original pipeline, and fix the gaps it exposes.

Content:
- `src/compressor.jl`: `compressor_power` / `condenser_power` accept any `AbstractRefrigerant`
  (`compression(r, pe, Tsuc, pc, ηis)`; spline backend keeps the original CoolProp-by-name calls,
  `ClapeyronRefrigerant` uses `h_pT`/`s_pT`/`h_ps`/`h_bubble`; condensing pressure = bubble pressure,
  subcritical only). Tests: `test/compressor.jl`.
- `examples/optimalflow_components.jl`: `study(r)` (continuation in flow 2 → 0.1 × design,
  post-processing, `validity_issues` per point), `report(rows)`, `original()`, `compare(new, old)`.

Gap found: **compile hangs** (the "~20 min at 100 % CPU" of the handoff was type inference, not the
solver). Fixed in `fix/concrete-solution-types` (PR #32, merged into this branch), found with JET
`@report_opt` and `Core.Box` checks of lowered code: untyped solution/parameter structs and
`SecondaryFluid` splines, abstract `NamedTuple[]` continuation results, non-inferable
`physical`, transport type and default formulation (now purity in the `ClapeyronRefrigerant`
type). Inference with an abstract refrigerant still blows up inside Clapeyron: keep refrigerants
and components concretely typed (beware Julia boxing a captured local keyword function, or a
variable assigned in an inner `@testset`, which is a closure).

Results (R290, MEA-28 %, Tb = −5 °C, Tc = 40 °C, ηis 0.7, ηpump 0.3; 39 flows 2 → 0.1 × design):

| variant | first call (compile) | sweep | max ΔTe vs original | max rel ΔCOP1 vs original |
|---|---|---|---|---|
| components, ClapeyronRefrigerant | 366 s | 1.2 s | 6.1e-5 K | 1.3e-6 |
| components, spline `Refrigerant` | 14 s | 0.3 s | 6.2e-5 K | 1.2e-6 |
| `original()` (old API) | 7.7 s total | | — | — |

- Clapeyron vs spline: |ΔTe| ≤ 7e-7 K, rel ΔCOP1 ≤ 1e-8. The residual difference to `original()`
  is the borehole density at its own mean temperature (B2). `original()` fails at the first point
  (2.0 × design, its hand-written guess) and recovers from 1.95; the component API converges at all.
- **COP optimum**: COP1 = 3.209 at 0.55 × design flow (0.27 l/s); COP2 = 2.302 at 0.65 ×. Both lie
  in the valid region: the brine leaves its table at 0.15 × and drops below its freezing point
  only at 0.10 × design flow.
- Jump between 1.65 and 1.60 × design flow (Te −0.8 K, Q2 −115 W, pump power 352 → 211 W):
  laminar–turbulent transition; the Nusselt blend is smooth, the pressure-drop correlation
  probably not. To check (follow-up).

Open:
1. First-call compile of the Clapeyron study (366 s) vs spline (14 s): profile with JET /
   SnoopCompile (`T_ps` alone ≈ 20 s); candidates: own Newton for `T_ps` of a pure fluid,
   PrecompileTools workload.
2. Pressure-drop continuity at the laminar–turbulent transition (jump above).
3. `test/tv_state.jl` (boxes `p`, `r`) and `test/property_interface.jl` (boxes `u`, 137 s
   runtime): same closure pattern as `operating_envelope.jl`.
4. Possible follow-ups: pump-power helper in the library; `Tb` as a port (step 8 note).
