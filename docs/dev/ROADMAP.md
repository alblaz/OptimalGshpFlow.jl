# Development roadmap and status

Last updated: 2026-11-18 (feature/ils-ground).

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
| 4 | `feature/component-core` | lightweight components (bounded unknowns + scaled residuals + ports), ComponentArrays + NonlinearSolve (ComponentArrays dropped in B8) | merged (#30) |
| 5 | `feature/glide-epsntu` | Qiao 2025 counter-flow ε-NTU with glide (eq. 56); pure fluid as the φ → 0 limit; Te := dew temperature | **raised priority** (2026-11-18): the refrigerant pressure drop gives a glide even for pure fluids, see below; R125 blocks only the zeotropic examples |
| 6 | `feature/evaporator-operating-envelope` (was `robust-moving-boundary`) | Te lower bound (spurious zero-flow root), validity reports (brine freezing, issue #14), continuation helper; PTC evaluated | merged (#31) |
| 6c | `feature/optimal-flow-components` | port of the optimal-flow study to the component API; compressor/condenser on the property interface; validated against the original (needs PR #32) | merged (#34) |
| B7 | `fix/pressure-drop-transition` | pipe friction factor jumped ×1.65 at Re = 2300 (pump power 344 → 216 W in the study); smooth blend 2300–4000, `SharpTransition()` keeps the original | merged (#36) |
| B8 | `fix/inferable-component-system` | `ComponentSystem` type not inferred → callers inferred NonlinearSolve/Clapeyron on abstract types; unknowns now a plain `Vector{Float64}` (no ComponentArrays); first call of the Clapeyron study 366 → 52 s | merged (#40) |
| T1 | `test/efficiency` | test suite 386 → 200 s (fresh session): shared ForwardDiff tag and Newton helper in `test/helpers.jl`, one spline `Refrigerant("R290")`, boxed test variables | merged (#41) |
| P1 | `perf/precompile-workload` | PrecompileTools workload: first call of the study 52.5 → 1.0 s, but breaks mixtures (ForwardDiff tag-count collision), see open item 1 | dropped |
| P2 | `perf/t-ps-newton` | `T_ps(…; phase = :vapor)`: Newton in T at constant p; Clapeyron's p-s flash behind an inference barrier; first `compressor_power` 17.9 → 0.2 s | merged (#43) |
| 6d | `feature/optimal-flow-vs-tb` | the optimal-flow study as a sweep over the ground temperature (`study_tb`, 9 → −5 °C) | merged (#45) |
| B9 | `fix/spline-study-compile` | spline `Refrigerant` had untyped fields (all property calls `Any`): first call of the spline study ≈ 450 → 4.9 s | merged (#46) |
| 8a | `feature/tb-port` | ground temperature `Tb` as a port of `BoreholeHX` (fixed node), `set_fixed` for new fixed values without rebuilding (step 8 groundwork) | merged (#44) |
| 6e | `feature/nusselt-gnielinski` | `GnielinskiTransition()`: Gnielinski (VDI 2013) turbulent, C¹ Hermite blend from laminar; sensitivity of the optimal flow to the transition model | merged (#48) |
| B10 | `fix/dittus-boelter-exponent` | Dittus-Boelter `mode = "heating"` used `Pr^0.3` (cooling exponent); now 0.4 | merged (#49) |
| 6f | `feature/gnielinski-default` | `GnielinskiTransition()` as the default transition model (heat transfer; friction unchanged) | merged (#50) |
| 4b | `feature/pump-power` | pump-power helper in the library: `pressure_drop`/`pumping_power` of an `EvaporatorSolution` or `BoreholeSolution` (open item 4) | merged (#51) |
| 5a | `feature/refrigerant-pressure-drop` | refrigerant-side Δp in the BPHE: Amalfi 2016 and Longo-type fits (friction), acceleration, gravity, ports, vapour zone; `surface_tension`, spline `ρ_bubble`. No model change (design: `docs/dev/glide.md`) | PR #55 |
| 6b | `feature/zone-switching` | zone appearance/collapse (Bell 2015 bounded-Q algorithm, Bonilla switching): needed for the condenser / R744 gas cooler or another expansion-valve control | deferred (cycle closure) |
| 7 | `spike/mtk` | how much ModelingToolkit integration is worthwhile | |
| 8 | `spike/dynamic-vs-steady` | is a dynamic model needed for robustness? coupling to a dynamic borehole model (warm starts, see below) | |
| 8b | `feature/ils-ground` | first dynamic ground model: infinite line source with temporal superposition, explicit coupling over a heating season (`examples/season_ils.jl`) | PR #52 |
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
- Glide (step 5) is not only for zeotropic blends (user, 2026-11-18): with the refrigerant-side
  pressure drop in the evaporator, the saturation temperature of a pure fluid falls along the
  two-phase zone, i.e. a glide of the opposite sign to a zeotropic one (temperature decreasing
  along the flow). The ε-NTU with glide must allow that sign; it needs a refrigerant pressure-drop
  model in the evaporator (see the old PR #12 "Pressure drop evap").
- First dynamic ground model (step 8, user, 2026-11-18): infinite line source,
  `ILS(Fo) = ½ E₁(1/(4Fo))`, `Fo = α t / r²`. Not good for multi-year dynamics, which are not the
  focus now.
- Composition framework: lightweight custom layer, not ModelingToolkit (MTK is evaluated in step 7).
- Laminar-turbulent transition: keep all models; `GnielinskiTransition()` (C¹) is the default since
  2026-11-18 (before: `SmoothTransition()`, B6). Results recorded before that date use the old default.
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
- **System:** unknowns in a `ComponentVector` (`u.evaporator.Qe`, `u.nodes.T_cold`, ...; a plain
  `Vector{Float64}` since B8, names through `physical`), square
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
  instead of rebuilding the system. **Done** (8a: `set_fixed`).
- Only if profiling shows the solver overhead matters: reuse the solver cache (`init`/`reinit!`)
  or a small dedicated Newton loop (≈ 0.1–0.2 ms per step), and/or the tangent predictor.

## `feature/optimal-flow-components` (results 2026-11-13)
Goal (agreed 2026-11-13): port `examples/heatpumpborehole_optimalflowrate.jl` (the study the package
is named after: COP vs brine flow, evaporator + borehole, R290, ground −5 °C, ΔTsup 3.5 K) to the
component API, validate it against the original pipeline, and fix the gaps it exposes.

Purpose of the study (user, 2026-11-16): how the optimal brine flow changes over the heating
season. The ground temperature `Tb` is highest at the start of the season and decreases during
extraction; this moves the heat-pump operating point (evaporation temperature, capacity) and the
brine viscosity, hence the flow regime, the heat transfer and the hydraulic performance. So the
study is a sweep over a range of `Tb`, not one value: the −5 °C of the original script is one end
(undisturbed ground is ≈ 8–9 °C in e.g. Stockholm; designs keep the brine ≥ 0 °C over the
lifetime). Next step: optimal flow vs `Tb` (e.g. 9 → −5 °C) with `set_fixed` (8a) and warm starts.

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
- Jump between 1.65 and 1.60 × design flow (Te −0.8 K, Q2 −115 W, pump power 352 → 211 W),
  checked on a 0.01 grid (2026-11-14, `fix/pressure-drop-transition`):
  - pump power: real discontinuity, `fDarcy` switched from 64/Re to the turbulent fit at
    Re = 2300 (×1.65). Fixed: linear blend in Re from 2300 to 4000 (`SmoothTransition()`, the
    `NusseltTransition` types now also select the friction model; `SharpTransition()` keeps
    the original). Pump power now 225.5 → 219.6 → 215.9 W across the crossing; optimum unchanged.
  - Te: continuous but steep (−0.26 K and −0.42 K per 0.01 × design flow at 1.64 → 1.62): the
    linear Nusselt blend goes from 3.66 at Re 2300 to Dittus-Boelter at 1e4 (≈ 100 for this
    brine's Pr), so h rises very fast just above 2300 and has a kink there. It makes a small
    spurious local COP1 maximum at 1.64 × (2.968, vs 2.941 at 1.62 ×). Not a bug; a smoother
    Nu model (e.g. Gnielinski 2013 with a C¹ blend) would remove it if it matters.
  - `f_Martin` (plate HX friction) also switches at Re = 2000, but the plate Re is only 3–56
    in this study.

Open:
1. First-call compile of the Clapeyron study (366 s) vs spline (14 s). **Main cause found and
   fixed (B8, 2026-11-16):** the `ComponentSystem` constructor was not inferable (port-node tuples
   from a generator, offsets in boxed captured variables, `ComponentVector(::NamedTuple)` which
   ComponentArrays cannot infer), so inside `study` the system, `u0` and everything downstream
   (`continuation` → `_solve` → NonlinearSolve, `evaporator_solution` → Clapeyron) were inferred on
   abstract types: ≈ 330 s of compilation for a 2-point study, although inference proper was only
   40 s (SnoopCompile) and each stage compiled at top level took < 25 s. JET `@report_opt` on the
   `study` body pointed at it directly. After the fix: first call 52.5 s, second 0.1 s (was 1.2 s).
   Remaining first-call costs (fresh session, stage by stage): `ClapeyronRefrigerant("R290")`
   22.5 s, `compressor_power` (`T_ps`) 17.9 s, first `solve` 7.8 s. A PrecompileTools workload
   was tried and dropped (P1, below).

   **`T_ps` (P2, 2026-11-16).** SnoopCompile on the first `compressor_power`: 15.8 of its 17.6 s
   of inference were `Clapeyron.PS.temperature`, and its result inferred as `Any`, so `T_ps`,
   `h_pT` and `compressor_power` were abstract too. Now: with `phase = :vapor` and a subcritical
   vapour state (`s ≥ s_dew(p)`), Newton in T at constant p (`ds = cp dT/T`, guess from the dew
   point with constant cp; iterations on primal values, final step with the Duals: exact
   derivatives). Other states (wet, supercritical, `:unknown`) still use the flash, called
   through `Base.inferencebarrier` (compiled only if used) with its result converted and asserted
   to the promoted number type. Note: a branch that is never taken is still inferred and compiled
   when its call is concrete; an inference barrier is the way to keep a heavy fallback out.
   First `compressor_power` after a solve: 17.9 → 0.2 s; `h_ps(:vapor)` 58 → 53 µs; Newton =
   flash within 6e-14 K; CoolProp agreement 1e-9 (R290, R134a, R32, R744).
   Pitfall seen while measuring: a closure capturing a non-`const` global refrigerant makes the
   same blow-up (515 s for `continuation`).

   **Precompile workload (P1, 2026-11-16, dropped).** PrecompileTools workload in the package
   (`ClapeyronRefrigerant("R290")`, saturation and single-phase properties, `T_ps`/`h_ps`,
   compression, the residual of an evaporator + borehole system) and in the NonlinearSolve
   extension (one `NewtonRaphson` and one `TrustRegion` solve, post-processing). Measured (Julia
   1.13, native CPU target):

   | | without | with workload |
   |---|---|---|
   | precompile, package + extension | 4.7 + 4.6 s | 40.2 + 12.8 s |
   | load of the example | 12.0 s | 5.3 s |
   | first call of the study (39 points) | 52.5 s | 1.0 s |
   | test suite | 190 s | ≈ 105 s (partly; tests failed) |

   **Why dropped:** it breaks every mixture. ForwardDiff orders the tags of nested Duals with
   `tagcount`, a `@generated` function that bakes the value of a global counter (`TAGCOUNT`,
   ForwardDiff 1.4.7 `src/config.jl`) into the compiled code. The counter restarts at 0 in every
   session, and the values baked into our package image during precompilation collide with tags
   created at runtime. Clapeyron's mixture critical point (`mixture_critical_constraint`: a Hessian
   nested in a gradient) then throws `DualMismatchError: Cannot determine ordering of Dual tags`, so
   constructing any mixture `ClapeyronRefrigerant` fails (8 tests). Pure fluids passed, no silently
   wrong value was seen. No safe fix on our side: it would mean changing ForwardDiff's internal
   counter, and Clapeyron uses ForwardDiff almost everywhere, so the workload cannot avoid baking
   tags. Possible later: (a) an upstream ForwardDiff fix (tag counts that do not restart per
   session; report it with a minimal reproducer), after which the workload can come back;
   (b) for fast first calls in one's own studies, a local PackageCompiler sysimage (a single
   image also stores the counter, so no collision expected; untested). The experiment was commit
   ae3fe2f on a local branch, not pushed.
2. ~~Pressure-drop continuity at the laminar–turbulent transition~~ (fixed, B7, see above).
3. ~~Test-suite time~~ (T1, 2026-11-16): 386 → 200 s in a fresh session (compile 354 → 185 s).
   Main cause: `ForwardDiff.derivative/gradient/jacobian` tag the Duals with the closure type, so
   every throwaway closure recompiled the model and Clapeyron for a new Dual type (≈ 0.8 s per
   closure for one property, much more for a model Jacobian). `test/helpers.jl` gives one tag
   (`dderivative`, `dgradient`, `djacobian`, `newton_solve`; fixed chunk 8) and one spline
   `SPL_R290` (3.3 s per construction, built 7 times). Boxed variables fixed (`u` in
   property_interface, `p` in tv_state). Per file: property_interface 131 → 37 s, components
   60 → 16 s, tv_state 35 → 9 s. Remaining: `properties_clapeyron.jl` ≈ 105 s, the first use of
   Clapeyron's code paths (same type for all pure fluids, nothing duplicated): a package
   precompile workload would cut it (≈ 102 → 42 s) but was dropped (P1, item 1). Rule for new tests:
   differentiate through the helpers, not `ForwardDiff.*` on closures.
4. ~~Pump-power helper in the library~~ (4b, 2026-11-18): `pressure_drop(es_or_bs, sf)` and
   `pumping_power(es_or_bs, sf, η)` for an `EvaporatorSolution` (plates, Martin, mean brine
   temperature) or a `BoreholeSolution` (U-pipe 2H, its `transition`); the example sums the two,
   same values as before. ~~`Tb` as a port~~ (8a, 2026-11-16):
   `BoreholeHX` has ports `Tin`, `Tout`, `Tb`; `Tb` is usually a fixed node
   (`T_ground = (:borehole => :Tb,)`, `fixed = (T_ground = …,)`), and `set_fixed(sys; T_ground = …)`
   returns a copy of the same type (≈ 1.6 µs, no allocation, nothing recompiled). Time steps with
   ΔTb = -0.1 K, warm start from the previous `u`: 1.3 ms and 3 Newton iterations per step (R290,
   NewtonRaphson). `BoreholeHX(…; Tb)` now throws an `ArgumentError` pointing to the node.

## `feature/optimal-flow-vs-tb` (results 2026-11-16)
`examples/optimalflow_components.jl`: `study_tb(r; Tbs = 9 → −5 °C)` runs the flow sweep of `study`
(2 → 0.1 × design) at each ground temperature. A continuation in `Tb` at the first flow gives the
warm start of each sweep. `optimum(rows, :COP1)` refines the best grid point with a parabola through
its neighbours (no extra solves; valid points only). `report_tb` prints the table. Shared
post-processing (`postprocess`) now also gives the borehole Reynolds number.

Results (spline R290, same settings as the study; 15 Tb × 39 flows, 1.7 s after compilation):

| Tb [°C] | Vf* COP1 | COP1* | Re at opt. | Te [°C] | Vf* COP2 | COP2* |
|---|---|---|---|---|---|---|
| 9 | 0.97 | 4.374 | 2940 | 0.3 | 1.00 | 3.478 |
| 5 | 1.07 | 3.972 | 2670 | −3.2 | 1.10 | 3.082 |
| 1 | 1.24 | 3.586 | 2610 | −6.5 | 1.26 | 2.714 |
| 0 | 0.67 | 3.515 | 1190 | −9.3 | 1.32 | 2.625 |
| −1 | 0.65 | 3.451 | 1130 | −10.1 | 0.72 | 2.543 |
| −5 | 0.57 | 3.209 | 770 | −13.7 | 0.64 | 2.302 |

- Two local COP maxima: one just above the laminar–turbulent transition (Re ≈ 2600–2900, the
  extra heat transfer pays for the pumping) and one in laminar flow (Re < 1200). In warm ground
  the transitional one wins and the optimal flow *rises* as the ground cools (0.97 → 1.24 ×
  design from 9 to 1 °C: the brine gets more viscous, so more flow is needed to stay above the
  transition). Between 1 and 0 °C (COP1) or 0 and −1 °C (COP2) the laminar maximum takes over and
  the optimal flow jumps down to ≈ 0.65–0.7 × design, then falls slowly (0.57 × at −5 °C).
- At −5 °C it reproduces `study` (COP1 = 3.209 at 0.570 ×, grid point 0.55 ×).
- The transitional maximum sits where the Nusselt blend (linear in Re from 2300) is steep, so its
  position depends on the transition model (see the Te kink of B7); a C¹ Nusselt blend (Gnielinski)
  is worth checking before drawing design conclusions.

Regression found here: the first call of the spline `study` took ≈ 445 s in a fresh session
(pre-existing on main); fixed in B9 (#46), see `fix/spline-study-compile` below.
## `fix/spline-study-compile` (2026-11-18)
First call of `study(Refrigerant("R290"))` (spline backend, fresh session, Julia 1.13.3): ≈ 450 s,
same on main 669de5b; the 14 s of the table above was measured earlier (probably on another Julia
version, not checked). Every stage compiled at top level took < 2 s (7 s in total), and
SnoopCompile counted only 20 s of inference, so the time was not in inference proper: a profile
put ≈ 210 of the 390 s of compilation in `store_backedges → _typename_add_backedge`, the
backedges Julia 1.13 registers for calls inferred on abstract types. Cause: `Refrigerant`,
`RefrigerantSaturatedVapor` and `RefrigerantSuperHeatedVapor` had untyped fields, so `T_dew`,
`h_pT`, … on a spline refrigerant inferred as `Any` and the whole study downstream was abstract.
JET `@report_opt` on the study did not finish in 15 min, which is itself a hint of this.
Fix: type parameters for the splines, `Float64` scalars. First call 4.9 s, second 0.11 s, results
unchanged (COP1 3.209 at 0.55 × design flow). Clapeyron study first call unchanged (36 s).
Test: `test/solutiontypes.jl` (concrete fields, `@inferred` `T_dew` / `h_pT`).
Side effect found by the tests: `initial_guess` named each block by searching `sys.blocks` with
`===`. With the refrigerant now stored inline in the (large, immutable) evaporator block, that
search returned `nothing` on Julia 1.13 although `b === b` was true (not pinned down; probably
egal on a copied inline struct). Now the names and blocks are zipped (same order by construction).

## `feature/nusselt-gnielinski` (2026-11-18)
Question from #45: is the jump of the optimal flow (transitional → laminar COP maximum) a product of
the linear Nusselt blend, which has a kink at Re = 2300?

New model `GnielinskiTransition()`: Gnielinski (VDI Heat Atlas, 2013 form, fully developed) for
`Re ≥ 1e4`, joined to the laminar 3.657 by a cubic Hermite blend in Re (value and slope continuous
at 2300 and 1e4, flat on the laminar side). Friction: same blend as `SmoothTransition`. The example
takes `bh = SingleUPipeBorehole(BOREHOLE; transition = …)`; pump power and Re use the borehole of
the solution.

Also found: `nusselt_dittusboelter(…; mode = "heating")` uses `Pr^0.3`. Dittus-Boelter uses 0.4 for a
fluid being heated (the brine in the borehole) and 0.3 for cooling, so the code looks reversed (the
test locks it). For the brine (Pr ≈ 30) the turbulent end of `SmoothTransition` is ≈ 30 % low:
Nu(1e4, Pr 32) = 103 vs 146 with `Pr^0.4`, 151 with Gnielinski. Fixed in B10 (below).

Optimal flow vs Tb (spline R290, 15 Tb × 39 flows, grid 0.05 × design), three transition models;
COP1 local maxima (flow × design, COP1, borehole Re):

| Tb [°C] | linear, DB `Pr^0.3` (current default) | linear, DB `Pr^0.4` | Gnielinski C¹ |
|---|---|---|---|
| 9 | 0.95, 4.374, Re 2940 | 0.95, 4.396 | 1.10, 4.305, Re 3380 |
| 5 | 0.80, 3.872 / **1.05, 3.970** | 0.80, 3.872 / **1.05, 3.995** | 0.80, 3.872 / **1.20, 3.877** |
| 3 | 0.75, 3.723 / **1.15, 3.778** | 0.75, 3.723 / **1.10, 3.802** | **0.75, 3.723** / 1.25, 3.675 |
| 1 | 0.70, 3.582 / **1.25, 3.586** | 0.70, 3.582 / **1.20, 3.612** | **0.70, 3.582** (only one) |
| −5 | **0.55, 3.209** / 1.65, 2.967 (Re 2444) | **0.55, 3.209** / 1.65, 2.989 | **0.55, 3.209** (only one) |

(bold: global maximum)
- **The two-branch structure is robust**: in all three models there is a laminar maximum (Re < 2200)
  and, in warm ground, a transitional one (Re ≈ 2500–3400); the optimal flow jumps from one to the
  other as the ground cools. The laminar branch is identical in all models (same laminar Nu and
  friction).
- **Where the jump happens is not**: COP1 between 1 and 0 °C (linear blends), between 5 and 4 °C
  (Gnielinski C¹); COP2 between 0/−1 °C, −1/−2 °C and 4/3 °C. The C¹ blend makes h rise slowly
  just above 2300, so the transitional optimum moves to higher Re (≈ 3000–3400, 1.1–1.2 × design)
  and is worth less.
- **The cold-ground transitional maximum of the linear blends is an artefact** of the kink at
  Re = 2300 (e.g. 1.65 × at −5 °C, Re 2444, the spurious maximum noted in #36): with the C¹ blend it
  is gone below 3 °C.
- **COP is flat between the branches near the switch**: e.g. 3.877 vs 3.872 at 5 °C (Gnielinski),
  3.586 vs 3.582 at 1 °C (linear): < 0.2 %. In practice, from ≈ 5 °C down the laminar flow
  (≈ 0.7–0.8 × design, lower pump power) costs almost nothing; the exact switch temperature is
  within the uncertainty of the transition heat transfer (± 4 K here).
- Grid 0.05 × design; a finer grid would sharpen the switch temperatures, not change the picture.

## `fix/dittus-boelter-exponent` (B10, 2026-11-18)
`nusselt_dittusboelter(…; mode = "heating")` used `Pr^0.3`; Dittus-Boelter uses 0.4 when the fluid
is heated (the default, and the borehole brine in heat extraction) and 0.3 when it is cooled. Now
0.4 / 0.3. Affects `SmoothTransition` (the default: its blend ends at Dittus-Boelter at Re = 1e4) and
every model above Re = 1e4 except `GnielinskiTransition`. Effect on the optimal-flow study = the
"linear, DB `Pr^0.4`" column above: laminar branch and cold-ground results unchanged; warm ground
COP1 +0.5 % (4.374 → 4.396 at 9 °C), optimal flow −0.02…−0.05 × design; COP1 switch still between
1 and 0 °C, COP2 switch moves from 0/−1 to −1/−2 °C.

## `feature/gnielinski-default` (6f, 2026-11-18)
Default `transition` of `SingleUPipeBorehole`, `nusselt_circularpipe`, `h_circularpipe`, `fDarcy`,
`pressure_drop` and `pumping_power`: `SmoothTransition()` → `GnielinskiTransition()`. The friction
factor is the same for both, so only the borehole heat transfer changes (Re ≥ 2300). The study
results under the new default are the "Gnielinski C¹" column of #48: COP1 optimum 1.10 × design
(4.305) at 9 °C, switch to the laminar branch between 5 and 4 °C, unchanged below ≈ 2 °C
(−5 °C: 0.55 ×, 3.209). The tables of #45 and B7 above were computed with `SmoothTransition()`.

## `feature/ils-ground` (8b, 2026-11-18)
`src/ground.jl`: `ILS(Fo)`, `ILS(t, r, α)` (SpecialFunctions `expint`, new direct dependency, already
in the manifest through Clapeyron), `InfiniteLineSource(; λ, α, rb, Tg)` and
`borehole_wall_temperature(g, tq, q, t)`: temporal superposition of piecewise-constant extraction
steps `q` [W/m] (O(steps) per call, no load aggregation). The wall temperature is the `Tb` port of
`BoreholeHX` (Hellström gives the fluid-to-wall part).

`examples/season_ils.jl`: `season(r; flow, days = 180, dt = 1 h)`. Explicit coupling: `Tb` at the
start of a step from the loads before it, `set_fixed(S; T_ground = Tb)`, Newton from the previous
`u`; the extraction `Qb/H` of the solution is the load of the step. `flow` fixed, or `:optimal`
(COP1 optimum of the flow sweep once a day, grid 0.05 × design).

Results (spline R290, continuous full-load operation, λ = 3, α = 1e-6 m²/s, rb = 0.115 m,
Tg = 8 °C, GnielinskiTransition):

| | fixed 1.0 × design | daily optimal flow |
|---|---|---|
| seasonal COP1 | 3.740 | 3.761 (+0.6 %) |
| heat delivered | 28 164 kWh | 27 960 kWh (−0.7 %) |
| electricity (compressor + pumps) | 7 531 kWh | 7 434 kWh |
| pump energy | 229 kWh | 143 kWh |
| Tb after 180 days | 2.89 °C | 2.94 °C |
| time (first run includes compilation) | 18 s | 32 s |

- Cost (season, compiled 10 s): superposition 0.1 s with `expint` at every step (9.3 M calls),
  0.008 s with the tabulated responses of `UniformLoadHistory` (uniform step: the response depends
  only on n − i; `ILS(k dt)` computed once, each step a dot product, same result to 1e-13).
  20 years hourly (175 200 steps): 13.6 s tabulated (≈ 165 s with the general method, O(n²)).
  Load aggregation (Claesson & Javed 2012, as in pygfunction; Bernier's MLAA) would make it
  ≈ O(n log n) and allow variable steps: postponed to the long-term (g-function) model.
  The rest of the 10 s is post-processing (2.2 ms per step, mostly outside the solve): to fix
  separately.
- Tb falls 2 K on the first day, 4 K in the first month, then slowly (≈ log t), 25–29 W/m.
- The optimal flow starts on the transitional branch (1.1–1.2 × design) and switches to the
  laminar one (0.8 → 0.75 ×) on day 7, at Tb ≈ 4.7 °C, as in the static Tb sweep (5 → 4 °C).
- No validity issues (brine stays above its freezing point).
- Caveats: continuous full-load operation (no load profile or on/off cycling); explicit coupling
  (one-step lag); `rb = 0.115` in the example is used as a radius, but 0.115 m is a typical
  borehole *diameter* (the package default is `rb = 0.115/2`), to check.
