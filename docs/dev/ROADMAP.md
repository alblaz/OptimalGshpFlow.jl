# Development roadmap and status

Last updated: 2026-11-21 (feature/cycle, handoff).

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
| T2 | `test/condenser-newton` | CI runtest +60 s on Julia min after 9b/9c (331 → 390 s): condensation testset ≈ 55 s (first compile of the R454C mixture transport, partly moved from later testsets; the fix is `perf/entropyscaling-barrier`), condenser ≈ 13 s (NonlinearSolve compiled again for the condenser system): condenser test on `newton_solve` (`chunk` keyword), −3.7 s locally | PR #67 |
| P1 | `perf/precompile-workload` | PrecompileTools workload: first call of the study 52.5 → 1.0 s, but breaks mixtures (ForwardDiff tag-count collision), see open item 1 | dropped |
| P2 | `perf/t-ps-newton` | `T_ps(…; phase = :vapor)`: Newton in T at constant p; Clapeyron's p-s flash behind an inference barrier; first `compressor_power` 17.9 → 0.2 s | merged (#43) |
| 6d | `feature/optimal-flow-vs-tb` | the optimal-flow study as a sweep over the ground temperature (`study_tb`, 9 → −5 °C) | merged (#45) |
| B9 | `fix/spline-study-compile` | spline `Refrigerant` had untyped fields (all property calls `Any`): first call of the spline study ≈ 450 → 4.9 s | merged (#46) |
| 8a | `feature/tb-port` | ground temperature `Tb` as a port of `BoreholeHX` (fixed node), `set_fixed` for new fixed values without rebuilding (step 8 groundwork) | merged (#44) |
| 6e | `feature/nusselt-gnielinski` | `GnielinskiTransition()`: Gnielinski (VDI 2013) turbulent, C¹ Hermite blend from laminar; sensitivity of the optimal flow to the transition model | merged (#48) |
| B10 | `fix/dittus-boelter-exponent` | Dittus-Boelter `mode = "heating"` used `Pr^0.3` (cooling exponent); now 0.4 | merged (#49) |
| B11 | `fix/suction-mass-balance` | with `Δp_suction ≠ 0` the compressor took ṁ = ρ(pe − Δp) Vref ≠ the evaporator's ρ(pe) Vref and h_in at (pe − Δp, T_suc): now the evaporator's ṁ and outlet enthalpy, isenthalpic suction (20 kPa: W −0.4 % → +4.7 %) | merged (#59) |
| B12 | `fix/initial-guess-inference` | `initial_guess` looped over the blocks: with four component types (closed cycle, 9d) the block widened to an abstract type and `guess` was inferred on `Any` (Clapeyron): 564 s first call; now recursion over the tuples, 0.3 s | merged (#68) |
| B14 | `fix/component-loops-inference` | same widening as B12 in `validity_issues` (loop over the blocks), `port_values` (runtime position in the block tuple) and `physical` (`_drop_empty` recursion: inference gives up from three entries): with four components (closed cycle) the post-processing of a blend was inferred on `Any` and never finished; recursion, name in the type (`Val`), generated `_drop_empty` | merged (#70) |
| B15 | `perf/entropyscaling-barrier` | EntropyScaling transport of blends inferred `Any`: through the superheat Δp (`Δp_suction`) it reached the compression, and `evaluate` of the open study (R454C) was inferred on `Any` with Clapeyron (no result in 30 min; 0.7 s with the barrier). Type assertion as `T_ph`. First blend solve unchanged (58.5 vs 59.5 s): the 70–85 s first call is not this | merged (#71) |
| B16 | `fix/node-transforms-inference` | `ComponentSystem` collected the node transforms in a `Vector` by port name: with transforms of different types (closed cycle: `LogBox`, `Linear`, `Above`, `Free`) the system built inside a function had an abstract type (`NTuple{11, Any}`), and everything called with it (`set_fixed`, Newton, post-processing) was inferred on it, Clapeyron included: `compare_cycle(r454c)` never started (9d handoff). Now a small lookup constant-propagated from the literal `nodes`, forced with `@constprop :aggressive` on the constructor (by default not propagated for the cycle) | PR #72 |
| 6f | `feature/gnielinski-default` | `GnielinskiTransition()` as the default transition model (heat transfer; friction unchanged) | merged (#50) |
| 4b | `feature/pump-power` | pump-power helper in the library: `pressure_drop`/`pumping_power` of an `EvaporatorSolution` or `BoreholeSolution` (open item 4) | merged (#51) |
| 5a | `feature/refrigerant-pressure-drop` | refrigerant-side Δp in the BPHE: Amalfi 2016 and Longo-type fits (friction), acceleration, gravity, ports, vapour zone; `surface_tension`, spline `ρ_bubble`. No model change (design: `docs/dev/glide.md`) | merged (#55) |
| 5b | `feature/epsntu-pressure-drop` | component `Evaporator`: `ϵe` eliminated (same unknowns with/without Δp), Qiao ε with the refrigerant pressure drop (default `Amalfi2016()`, `nothing` = original equations), `evaporator_pressure_drop`, compressor `Δp_suction`; `examples/evaporator_profiles.jl` (temperature profiles, compressor speed sweep) | merged (#56) |
| 5c | `feature/evaporator-zone-profiles` | `evaporator_zones` (U, capacities, NTU per zone); profiles of both zones from the model's NTUs (far ends reproduced to 1e-13 K) | merged (#57) |
| 5d | `feature/pressure-drop-studies` | effect of the refrigerant Δp on `study_tb` and `season` (COP1 −0.2 %, capacity −0.3 %, optimal flow unchanged; up to −0.5 % / −2.5 % capacity at 1.5× speed, Amalfi / Longo-type); `refrigerant_pressure_drop` kwarg and suction Δp in the examples | merged (#58) |
| 5e | `feature/blend-glide` | zeotropic glide in the two-phase zone (R454C): `T_r,in = T_ph(pe, hei)` (p-h flash), φ = C_s/C_r in Qiao's ε; `effectiveness_counterflow` finite at φ = 1; two-phase pinch against T_r,in | merged (#60) |
| 5f | `feature/blend-pressure-drop` | refrigerant Δp for blends (Amalfi default for all refrigerants): inlet at `T_ph(pe + Δp, hei)`, saturation shift in γ, glide in φ; blend surface tension (mole-fraction average) | merged (#61) |
| 5g | `feature/blend-boiling-htc` | nucleate boiling of blends: Thome & Shakir (1987) factor on Cooper (`mixture_boiling_factor`); R454C: h_ref −17…−23 %, capacity −0.6…−1.7 % | merged (#62) |
| 5h | `feature/blend-optimal-flow` | R454C in the optimal-flow study (same machine size as R290); condenser convention for blends (`:bubble`, `:midpoint` default, `:dew`) | merged (#63) |
| 9 | condenser and closed cycle | condenser on the water side (counter-flow, glide matching, as the evaporator): the fair blend comparison; plan and decisions below (9a–9d) | in progress |
| 9a | `feature/saturated-liquid-properties` | `cp_bubble`, `μ_bubble`, `λ_bubble` (Clapeyron pure/blend, spline via CoolProp) | merged (#64) |
| 9b | `feature/condensation-htc` | Longo 2010: Nusselt film (gravity) / Akers (forced convection), C¹ blend on Re_eq (Longo's criterion); Silver–Bell–Ghaly for blends | merged (#65) |
| 9c | `feature/condenser` | `Condenser` component: desuperheating + condensation, saturated liquid out; water side; `condenser_design`; `examples/condenser.jl` | PR #66 |
| 9d | `feature/cycle` | closed cycle: condensing pressure an unknown, `hei` from the condenser, compressor and expansion valve as links (`heatpump_cycle`, `Compressor`); `examples/cycle.jl` | PR (results below) |
| 10 | `feature/borehole-sizing` | borehole and collector sizing: hydraulic vs thermal optimum over borehole diameter, collector diameter/SDR and depth (larger holes and pipes for deep boreholes: lower linear losses, different Rb and flow regime), with the optimal flow; natural convection in groundwater-filled holes (Gustafsson & Westerlund; Spitler et al. 2016) | planned (user, 2026-11-20) |
| 6b | `feature/zone-switching` | zone appearance/collapse (Bell 2015 bounded-Q algorithm, Bonilla switching): needed for the condenser / R744 gas cooler or another expansion-valve control | deferred (cycle closure) |
| 7 | `spike/mtk` | how much ModelingToolkit integration is worthwhile | |
| 8 | `spike/dynamic-vs-steady` | is a dynamic model needed for robustness? coupling to a dynamic borehole model (warm starts, see below) | |
| 8b | `feature/ils-ground` | first dynamic ground model: infinite line source with temporal superposition, explicit coupling over a heating season (`examples/season_ils.jl`) | merged (#52) |
| – | `docs/dev-notes`, `ci/julia-versions` | this file and CLAUDE.md; CI on Julia `1`, `lts`, `min`, `pre` | merged (#24, #25) |

B4 (fixed in step 3): the pressure was mapped onto `(pr_min, pcrit)`, i.e. the reduced pressure
`pr_min` read as Pa (1e-2 Pa) instead of a pressure.

## Next steps (handoff 2026-11-19)
Done in the last sessions: Tb sweep and season (ILS ground, `UniformLoadHistory`), Gnielinski C¹
default, plate Δp fix, post-processing speed, refrigerant pressure drop (Amalfi default, Longo-type
fits) in the evaporator with Qiao's ε, profiles and speed sweep (`examples/evaporator_profiles.jl`;
design and results in `docs/dev/glide.md`).
Done (5d): the refrigerant pressure drop changes COP1 by only −0.2 % and capacity by −0.3 % at
design speed and leaves the optimal flow unchanged; at high compressor speed capacity drops up to
0.5 % (Amalfi) or 2.5 % (Longo-type) (`feature/pressure-drop-studies` below).
Candidates, roughly in order:
1. (was 2) Longo 2012 (R290 vaporisation, ASME J. Heat Transfer 134, 101801): the user will look for it;
   add its fit to `LONGO_FITS` and redo the Amalfi comparison (our G is only 3–13 kg/(m² s), at the
   low end of Amalfi's database).
2. Blends, next steps (5e–5h done): condenser on the water side (step 9: the COP comparison of a
   blend depends on it, see 5h); first-call compile of a blend solve (≈ 70–85 s: Clapeyron mixture
   code with Duals, EntropyScaling transport not inferable).
3. Validity: the two-phase pinch now compares the brine outlet with T_r,in (5e); a pinch inside a
   zone cannot occur (D monotonic along each zone), so this item is done.
4. Suction line (lumped Δp, L ≈ 0.3–0.5 m for packaged units), low priority, after B11.
5. The `rb` radius/diameter question of the ILS example (8b caveats); post-processing cost already
   fixed in #53 (2.2 ms → 30 µs per point).
6. Step 8 (dynamic vs steady), step 7 (MTK), 6b (zone switching), g-function ground model with
   load aggregation (Claesson & Javed) for multi-year runs.
Test pitfall to remember: in a `@testset`, a closure calling a local keyword function of the
testset (e.g. `map(...) do e; coupled(e); end`) boxes it and makes NonlinearSolve/Clapeyron infer
on abstract types (25 min); write straight-line code.
Example pitfall (5e): rows collected with `push!` into `[]` (a `Vector{Any}`) and a report function
calling model functions on their fields (`reference_error(x.S, x.u)`) infers the model, Clapeyron
included, on `Any`: 320 s of compilation, 32 GB. Build rows with a comprehension (concrete eltype).
Also put an inference barrier (type assertion) on non-inferable Clapeyron calls (`T_ph`: the p-h
flash returns `Any`).

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

## `feature/pressure-drop-studies` (5d, 2026-11-19)
Examples only. `study`, `study_tb`, `postprocess`, `season`, `optimal_point` take
`refrigerant_pressure_drop` (default `default_refrigerant_pressure_drop(r)`, Amalfi 2016 for a pure
fluid; `nothing` = original equations). Shared `evaluate(S, u, r, sf)`: with a Δp model the
superheat-zone Δp now lowers the compressor suction pressure (`Δp_suction`; before, the examples
left it at 0). Rows carry `Δp`, `Δp_sup`, `Tr_in`. `pressure_drop_tb(r)` / `report_pressure_drop_tb`
compare the two cases; `pressure_drop_grid` / `report_pressure_drop_grid`: capacity and COP1 over
Tb × compressor speed × brine flow for none / Amalfi / Longo-type.

Tb sweep (Clapeyron R290, 15 Tb × 39 flows, 7.6 s for both cases after compilation):

| Tb [°C] | Vf* without | COP1* without | Vf* with | COP1* with | Δp two-phase | T_r,in − Te | ΔCOP1* |
|---|---|---|---|---|---|---|---|
| 9 | 1.082 | 4.298 | 1.082 | 4.287 | 3.36 kPa | 0.23 K | −0.25 % |
| 5 | 0.769 | 3.868 | 1.187 | 3.860 | 2.96 kPa | 0.22 K | −0.20 % |
| 0 | 0.662 | 3.512 | 0.664 | 3.505 | 2.61 kPa | 0.22 K | −0.20 % |
| −5 | 0.561 | 3.206 | 0.562 | 3.200 | 2.29 kPa | 0.22 K | −0.18 % |

- The dew temperature at the evaporator outlet drops ≈ 0.09 K. The two-phase inlet is ≈ 0.22 K warmer
  (smaller mean ΔT to the brine), so COP1 falls by 0.16–0.25 % at every Tb and flow.
- The optimal flow does not move (≤ 0.003 × design). At 5 °C the transitional and laminar maxima
  are tied (COP1 3.8677 vs 3.8677 without, 3.8597 vs 3.8591 with): the jump 0.77 → 1.19 is a
  tie-break, not an effect of the pressure drop.
- Suction Δp of the superheat zone: 13 Pa (Tb 0 °C, design flow), no effect on W.
- Capacity (evaporator duty Q2 = ṁ Δh, ṁ = ρ_suction V_ref, user's point 2026-11-19): at design
  flow −0.29/−0.30 % at every Tb, i.e. about 1.5× the COP1 change. ρ_suction follows the dew
  pressure: ≈ 3.3 %/K of Te against ≈ 2.3 %/K for COP1 (W falls too). Over the whole flow grid
  −0.14 % (Vf 2, Tb −5) to −0.37 % (Vf 0.25, Tb −5): the low brine flows lose most.

Capacity vs compressor speed × brine flow (`pressure_drop_grid`, spline R290, 1.6 s; changes
relative to no refrigerant Δp; Longo-type = the upper estimate of Δp):

| Tb | speed | Vf | Q2 (no Δp) | Amalfi Δp | ΔQ2 / ΔCOP1 Amalfi | Longo-type Δp | ΔQ2 / ΔCOP1 Longo |
|---|---|---|---|---|---|---|---|
| −5 | 0.5 | 0.25 | 2015 W | 1.2 kPa | −0.19 / −0.11 % | 2.1 kPa | −0.33 / −0.18 % |
| −5 | 0.5 | 1.0 | 2109 W | 0.8 kPa | −0.01 / −0.01 % | 1.3 kPa | −0.01 / −0.01 % |
| −5 | 1.0 | 0.25 | 3449 W | 2.2 kPa | −0.37 / −0.19 % | 7.4 kPa | −1.23 / −0.62 % |
| −5 | 1.0 | 1.0 | 3849 W | 2.2 kPa | −0.29 / −0.16 % | 7.4 kPa | −0.98 / −0.53 % |
| −5 | 1.5 | 0.25 | 4557 W | 3.1 kPa | −0.52 / −0.24 % | 15.0 kPa | −2.54 / −1.17 % |
| −5 | 1.5 | 1.0 | 5269 W | 3.2 kPa | −0.46 / −0.24 % | 16.0 kPa | −2.30 / −1.16 % |
| 9 | 1.0 | 1.0 | 6006 W | 3.3 kPa | −0.30 / −0.24 % | 10.7 kPa | −0.97 / −0.76 % |
| 9 | 1.5 | 0.25 | 6536 W | 4.2 kPa | −0.51 / −0.30 % | 20.1 kPa | −2.40 / −1.43 % |
| 9 | 1.5 | 2.0 | 8923 W | 4.9 kPa | −0.44 / −0.34 % | 23.4 kPa | −2.14 / −1.62 % |

- Capacity loses about twice as much as COP1. The loss grows with compressor speed (G, Δp) and
  at low brine flow. With Amalfi it stays ≤ 0.5 %. With the Longo-type fit (Δp ∝ G²) it reaches
  2.5 % at 1.5× speed, so the choice of correlation (Longo 2012 R290 data) matters for capacity at
  high speed.
- At half speed and design flow or more, the capacity barely changes (≤ 0.01 %) despite
  0.5–1.5 kPa. Likely cause, not checked: the moving boundary gives the two-phase zone more area
  from an oversized superheat zone.
- Not modelled: the suction line between evaporator and compressor; it would lower ρ_suction
  (capacity) directly. Packaged units ≤ 10 kW have short lines (user, 2026-11-19: take the lower
  end, ≈ 0.3–0.5 m, a few bends). Estimate (R290, −10 °C, ṁ ≈ 0.015 kg/s, 1 m + 4 bends): 5/8" tube
  ≈ 1.5 kPa (0.12 K, capacity ≈ −0.4 %), 1/2" ≈ 4.4 kPa (0.37 K, ≈ −1.2 %); strongly dependent on
  D (friction ∝ D⁻⁵, fittings ∝ D⁻⁴). Same order as the two-phase Δp, not larger. The pressure
  drop inside a hermetic compressor belongs to its measured efficiencies (do not count it twice).
  Plan: B11 first (below), then optionally a lumped `(f L/D + ΣK) ρv²/2` line, off by default.
- B11 (found here): with `Δp_suction ≠ 0` the compressor took ṁ = ρ(pe − Δp, T_suc) Vref while the
  evaporator delivers ρ(pe, T_suc) Vref (mass not conserved), and its inlet enthalpy at (pe − Δp,
  T_suc) instead of the evaporator outlet enthalpy. Fix in `fix/suction-mass-balance`.

Season (spline R290, 180 days hourly, as 8b):

| | fixed 1.0 without | fixed 1.0 with | daily optimal without | daily optimal with |
|---|---|---|---|---|
| seasonal COP1 | 3.733 | 3.726 (−0.17 %) | 3.757 | 3.750 (−0.19 %) |
| heat delivered | 28 164 kWh | 28 104 kWh | 27 938 kWh | 27 881 kWh |
| electricity | 7 545 kWh | 7 542 kWh | 7 436 kWh | 7 435 kWh |
| pump energy | 244 kWh | 244 kWh | 145 kWh | 147 kWh |
| Tb after 180 days | 2.890 °C | 2.904 °C | 2.957 °C | 2.960 °C |

- The gain of the daily optimal flow is the same with and without (+0.65 %); switch to the laminar
  branch on day 5 in both.
- The "without" column differs from the 8b table (3.740, pumps 229 kWh, switch on day 7): changes
  merged since 8b, probably the plate Δp fix (pump energy); not investigated, the 8b numbers are outdated.
- Conclusion: for R290 at these mass fluxes (G ≈ 3–13 kg/(m² s)) the refrigerant pressure drop is a
  second-order effect on the optimal-flow question; it matters more for blends (glide) and higher G.

## `feature/blend-glide` (5e, 2026-11-19)
Zeotropic glide in the two-phase zone of the component `Evaporator` (design: `docs/dev/glide.md`,
item 6 of "Consequences"). R454C (R32/R1234yf 21.5/78.5 %, ≈ 8 K bubble–dew glide) as the case.
- `T_ph(r, p, h)`: Clapeyron p-h flash, with an inference barrier (`::float(promote_type(p, h))`;
  without it everything downstream inferred `Any`). ≈ 0.5–0.7 ms per call (Float64 or Dual).
- Two-phase zone: refrigerant enters at `T_r,in = T_ph(pe, hei)` (below `Te = T_dew(pe)`), linear
  temperature along the zone, pseudo capacity `C_r = ṁ (h_dew − hei)/(Te − T_r,in)`, `φ = C_s/C_r` in
  `effectiveness_counterflow`. Pure fluids: φ = 0 exactly, no flash (dispatch on the type).
  Decision (user, 2026-11-19): exact flash for T_r,in rather than linear bubble–dew interpolation
  (0.25 K too warm at the inlet for R454C, x_in ≈ 0.4).
- `effectiveness_counterflow` was 0/0 = NaN at φ = 1 (balanced counter-flow, crossed by a blend in a
  brine-flow sweep); now `ε = [E + γ G]/(1/NTU + φ E)` with series near `a = (1 − φ) NTU = 0`
  (test commit, then fix).
- `evaporator_solution`: ϵe referred to the coldest refrigerant temperature (≤ 1). `validity`: the
  two-phase pinch compares the brine outlet with `T_r,in` (blend: brine may leave below Te; pure
  fluid with Δp: slightly stricter than before). `evaporator_zones` reports `φ`.
- Cost (R454C residual 0.79 ms before): + one flash ≈ 0.7 ms (×1.9 residual). First call of a blend
  solve ≈ 85 s of compilation (pre-existing mixture cost dominates); the blend test is Float64 only
  (+7 s in the suite).

Results (`examples/blend_glide.jl`, standalone evaporator, brine in 5 °C, ΔTsup 3.5 K, design
Vref of the R290 example, sweep 2 → 0.3 × design brine flow, 0.9 s after compilation):

| Vf | Te | T_r,in | glide | φ | Qe | Ts,out | ΔQe exact profile |
|---|---|---|---|---|---|---|---|
| 2.0 | 1.44 °C | −3.97 °C | 5.40 K | 3.66 | 5909 W | 3.48 °C | −1.06 % |
| 1.0 | 1.04 | −4.35 | 5.39 | 1.85 | 5823 | 2.01 | −1.24 % |
| 0.5 | −0.65 | −5.99 | 5.33 | 0.98 | 5469 | −0.61 | −1.43 % |
| 0.3 | −3.00 | −8.25 | 5.25 | 0.63 | 5008 | −3.57 | −1.56 % |

- The zone glide is ≈ 5.4 K (inlet quality ≈ 0.4 of the 8 K bubble–dew glide). φ passes through 1
  near 0.49 × design flow, smoothly.
- The brine leaves below the dew temperature at design flow and below (Ts,out < Te): only possible
  with the glide; the pinch check uses T_r,in.
- Exact T(h) (flash at every RK4 step, `reference_error`) vs the linear profile: the exact curve is
  concave (T rises faster near the inlet), the refrigerant is warmer in mid-zone, so the linear
  model overestimates the two-phase duty by 1.1–1.6 % (brine end temperature off by ≤ 0.13 K).
  A two-segment linearisation would reduce it, if it matters.
- Superheat pinch at Vf ≥ 1.7: Te + ΔTsup reaches the 5 °C brine inlet (operating point, not a
  model issue).

## `feature/blend-pressure-drop` (5f, 2026-11-19)
- `surface_tension` of a blend: mole-fraction average of the pure components (CoolProp, cached
  states) at `T_bubble(p)`; neither Clapeyron nor CoolProp has a mixture surface tension. R454C at
  4 bar: R32 0.0131, R1234yf 0.0111 N/m, mix 0.0119 (mass-fraction average 0.0116). Amalfi depends
  on σ^0.22, so the mixing rule matters < 1 % in Δp. Amalfi needs no viscosity.
- Two-phase zone with both effects: `p_in = pe + Δp`, `T_r,in = T_ph(p_in, hei)`; saturation shift
  `Te − T_dew(p_in)` in γ, glide `T_dew(p_in) − T_r,in` in φ (`C_r = ṁ (h_dew − hei)/glide`). Pure
  fluids unchanged. `default_refrigerant_pressure_drop` is `Amalfi2016()` for every refrigerant.
- Cost: blend residual 1.2 ms (glide only) → 1.6 ms (with Δp: bubble state, σ, `T_dew(p_in)`).
- Check of the mixing rule (2026-11-20, paper from the user: Khordad 2018, Int. J. Refrig. 86, 73,
  `assets/docs/`): its Table 15 has measured σ of R32 + R1234yf (mole fractions 0.519/0.481,
  293–343 K, Zhao 2014). The mole-fraction average of the CoolProp pure values overestimates them
  by 7–10 % (real blends lie below the linear average: the lower-σ component enriches the
  surface). Effect: Δp ∝ σ^0.22, so ≈ +2 % in Δp and ≈ 0.01 % in capacity: negligible, kept.
  The paper's own method (Ornstein–Zernike/HNC with a Stockmayer potential) is not worth
  implementing: numerically heavy and less accurate for the pure fluids than CoolProp (R1234yf
  +17 %). If σ ever matters: Di Nicola & Pierantozzi 2013 (correlation for refrigerant binaries)
  or the corresponding-states model of Rice & Teja 1982 (≈ 0.3 mN/m mean deviation, their Table 18).

R454C, `examples/blend_glide.jl` (`pressure_drop_effect`), brine in 5 °C, ΔTsup 3.5 K:

| Vf | Δp | γ | T_r,in shift | ΔTe | ΔQe |
|---|---|---|---|---|---|
| 2.0 | 3.5 kPa | 0.025 | +0.21 K | −0.008 K | −0.03 % |
| 1.0 | 5.0 kPa | 0.035 | +0.24 K | −0.068 K | −0.25 % |
| 0.3 | 5.5 kPa | 0.031 | +0.23 K | −0.148 K | −0.56 % |

- Same order as R290 (5d: −0.3 % capacity at design): the pressure drop is a small correction for
  the blend too; the glide (5.2–5.4 K over the zone) dominates.
- The exact-T(h) reference with the pressure linear along the zone gives the same linear-profile
  error as without Δp (+1.1–1.6 % duty).

## `feature/blend-boiling-htc` (5g, 2026-11-20)
The refrigerant side of the two-phase zone is Cooper's nucleate pool-boiling correlation (×1.5),
with the mixture's reduced pressure and molar mass for a blend: the "ideal" coefficient. For a
zeotropic blend the light component evaporates first, the liquid at the bubble interface is
enriched in the heavy one and its local boiling point rises (mass-transfer resistance):
`mixture_boiling_factor` = Thome & Shakir (1987),
`F = [1 + (h_ideal/q) ΔT_bp (1 − exp(−B q/(ρ_l Δh_lv β_L)))]⁻¹`, `ΔT_bp = T_dew − T_bubble`,
`β_L = 0.0003 m/s`, `B = 1`; 1 for pure fluids (dispatch, no cost). Silver–Bell–Ghaly is the
correction for convective (flow-boiling/condensation) correlations, not for Cooper.
- Cost: one bubble and one dew state per residual (blends only).

R454C, `examples/blend_glide.jl` (`glide_sweep`, with Δp), brine in 5 °C, ΔTsup 3.5 K:

| Vf | q | h_ideal | F | h_ref | U | ΔTe vs 5f | ΔQe vs 5f |
|---|---|---|---|---|---|---|---|
| 2.0 | 8.2 kW/m² | 2783 | 0.767 | 2135 | 1536 | −0.17 K | −0.6 % |
| 1.0 | 6.5 | 2372 | 0.793 | 1881 | 1209 | −0.42 K | −1.6 % |
| 0.5 | 5.5 | 2067 | 0.815 | 1685 | 928 | −0.45 K | −1.7 % |
| 0.3 | 4.8 | 1846 | 0.833 | 1537 | 738 | −0.37 K | −1.4 % |

- The largest single correction for the blend so far (Δp: −0.03…−0.56 %). F rises towards 1 at
  low heat flux (`h_ideal/q` falls as `q^−0.33`, the exponential term as `q`).
- Not validated against blend data in a plate exchanger (none at hand); β_L = 0.0003 m/s is
  Thome's generic value. Cooper itself (pure, ×1.5) is the original model's choice; a flow-boiling
  plate correlation (e.g. Amalfi 2016 part 2, Longo) would be a separate step for all fluids.

## `feature/blend-optimal-flow` (5h, 2026-11-20)
`examples/optimalflow_components.jl`: `refrigerant_design(r; condenser)` gives `(hei, Vref, Tc)`.
Pure fluids: the R290 design of the example. Blend: the same machine size (nominal two-phase duty
`DESIGN.Qnom` at the design dew temperature, same exchanger, brine, borehole), liquid at the
condenser bubble point. `condenser_bubble_T(r, condenser)`: the condenser is fixed at 40 °C as
the bubble point (`:bubble`), the mean of bubble and dew (`:midpoint`, default) or the dew point
(`:dew`); R454C condenser glide 6.7 K → bubble 40.0 / 36.6 / 33.05 °C, `Vref` 2.29 / 2.20 / 2.11 l/s
(R290: 2.17). `design` is threaded through `heatpump`, `postprocess`, `study`, `study_tb`;
`evaluate`/`performance` default to `condenser_bubble_T(r)` (consistent with the default design).
`report_refrigerants_tb(pairs...)` prints the optima side by side.

COP1 optimum (valid points), R454C (Clapeyron, Δp, glide, Thome–Shakir) vs R290:

| Tb | R290 Vf* / COP1* | R454C bubble | R454C midpoint | R454C dew |
|---|---|---|---|---|
| 9 °C | 1.082 / 4.287 | 1.087 / 3.679 | 1.083 / 3.975 (−7.3 %) | 1.079 / 4.331 |
| 5 °C | 1.187 / 3.860 | 0.677 / 3.364 | 0.666 / 3.616 (−6.3 %) | 0.655 / 3.915 |
| 0 °C | 0.664 / 3.505 | 0.588 / 3.092 | 0.578 / 3.305 (−5.7 %) | 0.569 / 3.555 |
| −5 °C | 0.562 / 3.200 | 0.504 / 2.852 | 0.496 / 3.033 (−5.2 %) | 0.488 / 3.244 |

- The COP comparison is set by the condenser convention: R454C −14…−11 % (bubble), −7…−5 %
  (midpoint), +1…+1.4 % (dew). The evaporator side (5e–5g) moves Te by < 0.1 K between them. Only a
  condenser on the water side (step 9) decides it.
- Independent of the convention: on the laminar branch the optimal brine flow of R454C is
  8–13 % lower (0.58 vs 0.66 × design at 0 °C); the laminar maximum takes over at 5 °C (R290: 4 °C;
  6 °C with `:dew`); the dew temperature is 0.6–1.1 K higher (the glide follows the brine); Q2 within
  ±2 % (same nominal duty).
- Cost: R454C `study_tb` 116 s first call (70 s compilation), R290 11 s.

## Step 9: condenser and closed cycle (plan, 2026-11-20)
Decisions (user, 2026-11-20):
- Refrigerant-side condensation: Longo 2010 (Int. J. Refrig. 33, 944; `assets/docs/`), R290/R600a/
  R1270 in a BPHE (65°, b = 2 mm, the plate of `LONGO_FITS`; the condensation Δp fit ζ = 1900 is
  already there). Two regimes, transition at Re_eq ≈ 1600–1700 (G ≈ 15–18 kg/(m² s)): Nusselt film
  condensation × Φ (gravity controlled; explicit in the heat flux, h ∝ q^−1/3; 14.8 % mean
  deviation) and Akers × Φ (forced convection, local in x; 21 %). Joined by Longo's criterion (Re_eq
  at the mean quality, C¹ blend 1500–1800; a smooth maximum was tried and rejected, see 9b). Our G
  (3–13) is in the gravity regime. Blends: Silver–Bell–Ghaly
  (1/h = 1/h_film + Z_g/h_g, Z_g = x cp_g dT/dh).
- Zones: desuperheating + condensation, saturated liquid at the outlet (as the current cycle);
  subcooling later with zone switching (6b).
- Closed cycle: the condenser is a component, the condensing pressure an unknown set by the water
  side, `hei` from the condenser outlet. Expansion valve: isenthalpic, superheat set by `ΔTsup`
  (as now). Compressor: fixed displacement and isentropic efficiency (as now); a better compressor
  model is of interest, later: volumetric efficiency (mass flow below ρ_suc × displacement), and
  for speed-modulating compressors the inverter efficiency (user, 2026-11-20).
- Water side: floor heating 30 → 35 °C at design (for now); Wanniarachchi. Geometry: the plate
  type of the evaporator, area from a condenser analogue of `EvaporatorDesign`.
Tools checked: saturated-liquid ρ, h available; μ_L, λ_L, cp_L work (CoolProp for pure fluids,
exact; EntropyScaling for blends) and enter the interface in 9a. Desuperheating: the vapour
correlation of the evaporator superheat zone; the spline covers 40 K of superheat (discharge
25–35 K above condensation; configurable).
Found in 9a: EntropyScaling transport (blends: `μ_pT`, `μ_dew`, `μ_bubble`, …) infers `Any`
(`_transport` passes EntropyScaling's result through without a barrier). An assertion as in
`T_ph` would fix it, and may cut the 70–85 s first-call compile of a blend solve: separate branch
(`perf/entropyscaling-barrier`), to measure.

## `feature/condensation-htc` (9b, 2026-11-20)
`src/condensation.jl`: `nusselt_film_htc` (written with the heat flux: `h = [0.943⁴ λ³ρ²gΔh_lv/(μ L q)]^⅓`),
`akers_htc` (mean over the quality in closed form: Re_eq linear in x), `reynolds_equivalent`,
`condensation_htc(Longo2010(), r, phx, m, p, x_in, x_out, q; L)` → `(; h, gravity, forced, Re_eq)`.
h and q refer to the actual area `A` (as the evaporator); Longo's projected-area values are Φ × h,
Φ × q (so Φ drops out of his formulas here). `d_h = 2b` as Longo. Blends: Silver–Bell–Ghaly with
the vapour alone (Wanniarachchi) at the mean quality.

Check against Longo's data (his BPHE, x 1 → 0, h × Φ; reference = my reading of his Figs. 3–5):

| fluid | Tsat | G = 6 | G = 10 | G = 15 | G = 20 | G = 28 |
|---|---|---|---|---|---|---|
| R290 | 25 °C | +6 % | −11 % | +42 % | +32 % | +13 % |
| R290 | 40 °C | 0 % | −16 % | +30 % | +21 % | +4 % |
| R600a | 25 °C | −18 % | −31 % | +23 % | +33 % | +29 % |
| R1270 | 25 °C | +15 % | −3 % | +35 % | +33 % | +19 % |

- Gravity regime (G ≤ 10, Nusselt): within Longo's 14.8 % mean deviation for R290 and R1270, R600a
  −20…−30 % (the data are flat in G, Nusselt falls with q ∝ G). Forced (G ≥ 20, Akers): +4…+33 %
  (Longo: 21 %, Akers over-predicts). Just above the transition (G 15–20) Akers jumps while the
  data rise gradually: +16…+42 %.
- First version: the smooth maximum (p-norm, n = 4) of the two: +22…+41 % at G = 6–15 (adds them
  where they coincide; takes Akers below the transition). Replaced by Longo's Re_eq criterion.
- Our condenser (the evaporator's mass flow, G ≈ 3–13) is in the Nusselt regime. If the transition
  matters (high compressor speed): Longo's J_H(Re_eq) curve (his Fig. 6, read off the figure).

## `feature/condenser` (9c, 2026-11-20)
Done:
- `src/components/condenser.jl`: `Condenser(; refrigerant, secondaryfluid, plateheatexchanger, Vs,
  Qscale, condensation = Longo2010())`. Ports `Tw_in`, `Tw_out`, `m`, `h_in` (refrigerant mass flow
  and inlet enthalpy: from the compressor in 9d; fixed nodes standalone). Unknowns `pc` (LogBox),
  `Qd`, `Qc`, `xd`, `Tw_m`. Six residuals: refrigerant and water balance and counter-flow ε-NTU per
  zone (`effectiveness_counterflow`; desuperheat with φ = C_w/C_vapour, condensation with φ for a
  blend's glide). U: Longo condensation (over the zone length, q = Qc/A_zone floored at 1e-3 Qscale)
  ‖ water (Wanniarachchi); desuperheat: vapour ‖ water (Wanniarachchi). `vapour_T_ph` (Newton on
  h_pT, no flash), `condenser_zones`, `validity` (pinches at both ends and at the boundary),
  `condenser_design(Q, Tc, Tw_in, Tw_out, water; U)`. Zones inferred concretely (spline, Clapeyron).
- `test/condenser.jl` (12 tests, spline R290, < 1 s after compilation): converges at the design
  point, energy balances, Td − Tw_m = (Td − Tw_in) e^−NTU, validity. Full suite passes (247 s, fresh session).
- R290 design point (discharge from Te 0 °C, ΔTsup 3.5 K, ηis 0.7; Q1 7641 W; A 0.530 m², water
  0.368 l/s, 30 → 35 °C): Td 41.38 °C (sized for 40 at U 2000; actual U 1724), Tw_out 34.94 °C,
  T_dis 58.4 °C, xd 0.184, Qd 836 W (11 %); condensation h_ref 2044, U 1724; desuperheat h_ref 669
  (vapour side limits), U 632.
- Note: G = 21 kg/(m² s) in this condenser (smaller area than the evaporator: fewer channels), i.e.
  in Longo's transition/forced region where Akers over-predicts (~ +13 % at G = 21 vs his curve).
  The J_H(Re_eq) fit of his Fig. 6 would fix it if needed.
- `examples/condenser.jl` (written, first run was compiling at handoff): `report_condenser([R290,
  R454C])` (R454C in the R290-sized condenser at its own design discharge: where bubble, dew and
  mean temperature fall against the 30 → 35 °C water, vs the 5h conventions) and `water_sweep`.
  Results (first run 57 s, sweep 6 s): R290-sized condenser (A 0.528 m², water 0.366 l/s):
  | | pc | T_bubble | T_dew | mean | Tdis | Q | h_ref cond. | U cond. | G |
  |---|---|---|---|---|---|---|---|---|---|
  | R290 | 14.14 bar | 41.40 | 41.40 | 41.40 °C | 58.5 °C | 7521 W | 2041 | 1721 | 20.9 |
  | R454C | 17.07 bar | 39.79 | 46.45 | 43.12 °C | 65.1 °C | 7621 W | 1556 | 1363 | 41.9 |
  The water side puts R454C at bubble ≈ 40 °C: the `:bubble` convention of 5h (the least
  favourable, COP −11…−14 % vs R290), not the midpoint: a 6.7 K glide against a 5 K water rise
  does not match, and its condensation coefficient is lower (twice the mass flux, lower latent heat,
  Silver–Bell–Ghaly). `water_sweep` (Tw_in 25 → 40 °C): T_bubble − Tw_out ≈ 4.6–4.9 K throughout,
  xd 0.155–0.167. To be confirmed in the closed cycle (9d).

Next: 9d (closed cycle: evaporator +
compressor link + condenser + borehole; `hei` from the condenser outlet, pc from the water side;
node scaling for `m`, `h_in` ports, now `Free()`); then `perf/entropyscaling-barrier`.

## `fix/borehole-groundwater` (B13, 2026-11-20)
The study borehole (`examples/optimalflow_components.jl`) had three geometry errors: `rb = 0.115`
(the 115 mm diameter as the radius: a 230 mm borehole), `rpi = 0.04 (1 − 1/17)/2` (SDR 17 wall halved:
ID 37.6 instead of 35.3 mm) and no `rpo` (struct default 0.0223, a 44.6 mm OD). Found from the
effective borehole resistance, 0.193 m K/W (user: far too large). Now `borehole(; Db = 0.115,
Dpo = 0.040, SDR = 17, λfill = 0.57, H = 200., gap = 0.0025)`: the Swedish reference, 115 mm
groundwater-filled (user, 2026-11-20; stagnant water, conduction only), PE 40 × 2.4, legs 2.5 mm
from the wall (±35 mm). `ORIGINAL_BOREHOLE` keeps the old geometry for `original()`/`compare`.
Naming (user): `Tb` borehole wall, `Tg` undisturbed ground (as `InfiniteLineSource`); the node
`T_ground` of examples, tests and docstrings is now `Tb`. `season_ils` takes `rb` from `BOREHOLE`:
the ILS radius item of 8b is fixed with it.

`study_tb(R290)`, optimum of COP1 (old geometry → reference borehole):

| Tb [°C] | Vf* | Vs* [l/s] | COP1* | Te [°C] |
|---|---|---|---|---|
| 9 | 1.082 → 1.008 | 0.531 → 0.495 | 4.287 → 4.221 | 0.32 → −0.28 |
| 5 | 1.187 → 0.698 (laminar) | 0.582 → 0.343 | 3.860 → 3.823 | −3.39 → −5.41 |
| 0 | 0.664 → 0.600 | 0.326 → 0.294 | 3.505 → 3.477 | −9.37 → −9.69 |
| −5 | 0.562 → 0.508 | 0.276 → 0.249 | 3.200 → 3.178 | −13.79 → −14.10 |

- COP1 −0.7…−1.5 %, optimal flow −7…−10 %; the laminar optimum takes over at 5 °C (was 4 °C).
- The corrections cancel only partly: the smaller borehole lowers the filling resistance, but
  stagnant water (0.57) conducts less than the grout (2), and the smaller ID raises the pump power
  (Re +7 % at the same flow helps the heat transfer). Natural convection in the groundwater would
  lower Rb: a candidate for step 10 (borehole sizing).

## `feature/cycle` (9d, 2026-11-20)
Closed cycle as one `ComponentSystem` (design agreed with the user, 2026-11-20): links are ports
with residuals in the component that computes them; the isenthalpic valve is a shared node.
- `src/components/cycle.jl`: `CyclePorts(r)` (transforms of the refrigerant-loop ports, the same in
  every component of a refrigerant: pressure `LogBox` over `pressure_bounds`, enthalpy
  `Linear(h_bubble, Δh_lv)` at 0 °C, mass flow `Above(0, 1)` (log scale); residual scales);
  `Compressor(; refrigerant, ηis, ξ, Δp_suction)` (ports `p_suc`, `T_suc`, `m`, `pc`, `h_dis`; one
  residual), `compression_state` (property interface, AD-able; spline: `vapour_h_ps`, Newton on
  `s_pT`), `compressor_power(c::Compressor, p)`.
- `Evaporator(…; hei = nothing, hei_nominal)`: ports `hei` (in), `p_suc`, `T_suc`, `m` (out, three
  residuals). Dispatch on the `cycleports` field: a method on `Evaporator{…, Nothing}` was not taken
  as more specific than `ports(::Evaporator)` (constrained type variable `H <: Union{Nothing, T}`).
- `Condenser`: `pc` is now a port (was an internal unknown), new port `h_out` (`h_bubble(pc)`).
- `src/components/heatpump.jl`: `heatpump_cycle(; evaporator, compressor, condenser, borehole,
  Tb, Tw_in)` (nodes `T_cold`, `T_warm`, `Tb`, `p_suc`, `T_suc`, `m`, `pc`, `h_dis`,
  `h_liq`, `Tw_in`, `Tw_out`; 19 unknowns), `heatpump_guess(sys; Te, Tc, T_cold, T_warm, Tw_out)`.
- `test/cycle.jl` (16 tests, spline R290, `newton_solve`): energy balance to 1e-10, W and Q1 equal the
  open-cycle post-processing (`compressor_condenser_power`, CoolProp compression) within 1e-7.
- R290 spline, old study borehole (before B13), Tb 8 °C, water 30 °C in: Te −1.37 °C, condensation 41.1 °C, Tdis 59.7 °C, water out
  34.80 °C, Q2 5703 W, W 1753 W, Q1 7334 W, COP1 (no pumps) 4.18. 6 Newton steps from the guess;
  residual 19 µs, second solve 5 ms; `residual!` inferred, JET clean except a pre-existing runtime
  dispatch in `uniform_Tb_coeffs` (keyword call of `h_circularpipe` with an abstract `transition`):
  separate fix.
- `heatpump_guess` first call took 564 s (Clapeyron R290): `initial_guess` inferred `guess` on
  `Any` with four component types; fixed in B12 (#68), now 0.3 s; first closed-cycle solve 19 s.

### 9d handoff (2026-11-21)
State: `feature/cycle` pushed (not a PR). It contains the 9d work and merges of the open PRs #70 (B14)
and #71 (B15); after they are merged into main, merge main into it.
- Closed cycle works for R290 and R454C (`cycle_newton` in `examples/cycle.jl`: plain Newton, fixed
  ForwardDiff tag; R454C 5 iterations, 77 s first call incl. compilation). `cycle_point` R454C 0.07 s.
- R290 results (reference borehole, water 30 °C in, design brine flow), closed vs open (Tc = 40 °C):
  | Tb | Te | T_cond | Tdis | Tw_out | Q1 | W | COP1 closed | open |
  |---|---|---|---|---|---|---|---|---|
  | 8 | −1.25 | 41.20 | 59.8 | 34.83 | 7356 | 1757 | 4.006 | 4.110 (+2.6 %) |
  | 0 | −9.05 | 39.49 | 61.5 | 33.93 | 5976 | 1645 | 3.452 | 3.418 (−1.0 %) |
  | −5 | −13.20 | 39.75 | 64.0 | 33.48 | 5304 | 1601 | 3.093 | 3.079 (−0.5 %) |
  To check: T_cond rises from 0 to −5 °C while Q1 falls (Longo gravity/forced transition at lower G?).
- Resolved (2026-11-21): `compare_cycle(r454c)` hung in inference for two reasons. (1) B16 (#72):
  the closed-cycle system type was abstract when built inside a function (node transforms in a
  `Vector`), so `set_fixed`, `cycle_newton` and `cycle_point` were inferred on it; fixed in
  `ComponentSystem` (`@constprop :aggressive`), locked by `@inferred heatpump_cycle` in
  `test/cycle.jl`. (2) `open_point` passed `sol.u` of `solve(NonlinearProblem…)` (not inferred) to
  `evaluate`/`performance`: inference of `open_point(R454C)` > 13 min, 4.8 s with a
  `convert(Vector{Float64}, sol.u)::Vector{Float64}` barrier (as the extension's `_solve`).
  `refrigerant_design` was fine. Not rechecked: whether the 30-min `solve` of the blend closed cycle
  (reason for `cycle_newton`) was B16 too.
- Results (`compare_cycle`, R290 31 s, R454C 135 s incl. compilation; open cycle COP1 and Δ vs closed):
  | | Tb | Te | T_bub | T_dew | Tdis | Tw_out | Q1 | W | COP1 closed | bubble | midpoint | dew |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | R290 | 8 | −1.25 | 41.20 | 41.20 | 59.8 | 34.83 | 7356 | 1757 | 4.006 | | 4.110 (+2.6 %) | |
  | R290 | 0 | −9.05 | 39.49 | 39.49 | 61.5 | 33.93 | 5976 | 1645 | 3.452 | | 3.418 (−1.0 %) | |
  | R290 | −5 | −13.20 | 39.75 | 39.75 | 64.0 | 33.48 | 5304 | 1601 | 3.093 | | 3.079 (−0.5 %) | |
  | R454C | 8 | −0.40 | 39.84 | 46.50 | 67.3 | 35.04 | 7671 | 2064 | 3.579 | 3.547 (−0.9 %) | 3.823 (+6.8 %) | 4.151 (+16.0 %) |
  | R454C | 0 | −8.18 | 37.77 | 44.52 | 68.6 | 34.02 | 6114 | 1861 | 3.142 | 3.011 (−4.2 %) | 3.209 (+2.1 %) | 3.439 (+9.4 %) |
  | R454C | −5 | −12.42 | 36.76 | 43.55 | 69.6 | 33.53 | 5371 | 1742 | 2.895 | 2.741 (−5.3 %) | 2.902 (+0.3 %) | 3.088 (+6.7 %) |
  Open-cycle Tc: bubble 40.00, midpoint 36.60, dew 33.05 °C. Closed cycle R454C vs R290: COP1
  −10.7 / −9.0 / −6.4 %, Q1 +4.3 / +2.3 / +1.3 % (Tb 8 / 0 / −5 °C).
  - With the water side, R454C condenses with bubble 37–40 °C, dew 44–46 °C. No fixed convention
    fits the whole range: bubble is closest at 8 °C (−0.9 %) and drifts to −5.3 % at −5 °C, midpoint
    goes from +6.8 % to +0.3 %, dew overestimates throughout (+7 to +16 %). The open-study
    blend comparison (5h, midpoint) favoured R454C by up to 7 % at warm ground; in the closed cycle
    R454C is 6–11 % below R290 in COP1.
  - Still to check: T_cond (R290) rises from 0 to −5 °C while Q1 falls (also the R454C bubble from
    0 to −5 °C falls only 1 K): Longo gravity/forced transition at lower G?
- Lessons (inference with 4 components / blends): B12, B14, B15, B16 were all "abstract type reaches
  Clapeyron mixture code" (minutes to hours of compilation, looks like a hang). Before a long run,
  check the entry points with `Base.return_types`/`@inferred`; log progress to a file
  (`compare_cycle(r; log = io)`) and follow it with a monitor (never delete a file Julia has open).
  Stopping an MCP call does not stop Julia: `julia_restart` the session (or the process keeps running).
  Open the log in its own call: `open(log) do io … end` around the runs compiles the whole closure
  first, and the log stays empty during compilation (looks like a hang).
