# Development roadmap and status

Last updated: 2026-11-12.

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
| 2 | `feature/property-interface` | refrigerant interface + (T, v) state formulation, see below | **in progress** (sub-step 1 done) |
| 3 | `feature/bounded-variables` | transforms as types (`LogBox` for pressure, `Box`, `Unit`); fix the absolute `pr_min` bound (B4) | |
| 4 | `feature/component-core` | lightweight components (bounded unknowns + scaled residuals + ports), ComponentArrays + NonlinearSolve | |
| 5 | `feature/glide-epsntu` | Qiao 2025 counter-flow ε-NTU with glide (eq. 56); pure fluid as the φ → 0 limit | |
| 6 | `feature/robust-moving-boundary` | zone appearance/collapse, continuation, pseudo-transient continuation | |
| 7 | `spike/mtk` | how much ModelingToolkit integration is worthwhile | |
| 8 | `spike/dynamic-vs-steady` | is a dynamic model needed for robustness? | |
| – | `docs/dev-notes`, `ci/julia-versions` | this file and CLAUDE.md; CI on Julia `1`, `lts`, `min`, `pre` | merged (#24, #25) |

Known bug, not fixed yet (B4): `models.jl` maps the pressure with
`sigmoid_R_to_ab(u, refrigerant.pr_min, refrigerant.pcrit)`, so the lower bound is
`pr_min` read as Pa (1e-4 or 1e-2 Pa) instead of `pr_min·pcrit`. This belongs in step 3.

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

## Plan: `feature/property-interface` (next)
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
3. **Saturation without inner iterations** (pure fluids). Options:
   - Clapeyron superancillaries, if available for `MultiFluid`. To verify: `SuperAncSaturation`
     may only cover cubic/SAFT models;
   - otherwise (T, v_l, v_v) as extra unknowns with equal pressure and Gibbs energy, 2–3
     equations per saturation state.

   Mixture dew/bubble points (n + 2 unknowns each) go to the glide branch (step 5).
4. **Augmented evaporator**, behind an option. Add the superheat-outlet state (and optionally
   the saturation state) as unknowns with the residuals above. Compare with the nested
   formulation (properties solved inside the residual) on:
   - the 42-case sweep: success rate and iterations, from a fixed guess and with continuation;
   - cost per residual evaluation and per solve;
   - Jacobian conditioning.

   Report the results. Adopt the formulation in `component-core` only if it pays off.

Out of scope for this branch: residual scaling of the existing equations (step 4 of the
roadmap), glide (step 5) and new components.
