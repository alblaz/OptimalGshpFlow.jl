# Development roadmap and status

Last updated: 2026-11-13 (feature/bounded-variables).

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
| 4 | `feature/component-core` | lightweight components (bounded unknowns + scaled residuals + ports), ComponentArrays + NonlinearSolve | **next** |
| 5 | `feature/glide-epsntu` | Qiao 2025 counter-flow ε-NTU with glide (eq. 56); pure fluid as the φ → 0 limit; Te := dew temperature | lower priority (natural refrigerants are pure); R125 blocker, see below |
| 6 | `feature/robust-moving-boundary` | zone appearance/collapse, continuation, pseudo-transient continuation | |
| 7 | `spike/mtk` | how much ModelingToolkit integration is worthwhile | |
| 8 | `spike/dynamic-vs-steady` | is a dynamic model needed for robustness? | |
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
