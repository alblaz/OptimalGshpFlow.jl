# Step 5: ε-NTU with glide and refrigerant pressure drop (design note)

Status: decisions Q1–Q4 taken (user, 2026-11-18, see the end); waiting for the correlation papers.

## Why (user, 2026-11-18)
The refrigerant-side pressure drop in the evaporator lowers the saturation temperature along the
flow, so even a pure fluid (R290) evaporates with a temperature change ("glide" from pressure
drop). Blends add the thermodynamic glide. The current model has a constant evaporation
temperature `Te` in the two-phase zone and no refrigerant pressure drop at all.

## Reference: Qiao 2025 (Appl. Therm. Eng. 279, 127896), `assets/docs/epsilon-ntu glide.pdf`
Two parameters, both assumed linear along the area:
- **Glide** (blends): pseudo two-phase heat capacity `c_p,TP = (h_dew − h_bub)/(T_dew − T_bub)`
  (eq. 3), `C_r = ṁ c_p,TP`, capacity ratio `φ = C_f / C_r` (`φ = 0` for a pure fluid).
- **Pressure drop**: `ΔT_sat = T_sat(p_out) − T_sat(p_in)`, the saturation-temperature shift over
  the zone (eq. 6), `γ = ΔT_sat / (T_r,in − T_f,in)`.

Effectiveness `ε = (T_f,out − T_f,in)/(T_r,in − T_f,in)` (eq. 46), counter-flow:
- general (eq. 56):
  `ε = [1 + γ/((1−φ)NTU) − (γ + γ/((1−φ)NTU) + 1) e^{−(1−φ)NTU}] / (1 − φ e^{−(1−φ)NTU})`
- pure fluid with pressure drop, `φ → 0` (eq. 58): `ε = (1 − e^{−NTU})(1 + γ/NTU) − γ e^{−NTU}`
- no pressure drop, `γ → 0`: the single-phase counter-flow formula with `C_r` (blends), and the
  classical `1 − e^{−NTU}` (current model) when also `φ → 0` (eq. 59).

In an evaporator `ΔT_sat < 0` and `T_r,in − T_f,in < 0`, so `γ > 0`: the sign is built in. The
refrigerant enters the two-phase zone *warmer* than it leaves (higher pressure), facing the
outgoing brine.

## Consequences for the model
1. **Two-phase zone residuals** (`evaporator_balances!`, `du[3:4]`): `Q_e = ε C_s (T_s,m − T_r,in)`
   with `T_r,in = T_sat(p_in)` and `ε` from eq. 56/58, instead of `ε C_s (T_s,m − T_e)` with
   `ε = 1 − e^{−NTU}`.
2. **Bound of `ϵe`**: with `γ > 0` the effectiveness can exceed 1 (e.g. `γ = 1`, `NTU = 2` →
   `ε ≈ 1.16`); the unknown is mapped onto `(0, 1)` now (`Unit()`). Options: drop `ϵe` as an
   unknown (it is explicit in NTU and γ), or a wider `Box`. Recommendation: eliminate it (one
   unknown and one residual fewer, no bound to choose).
3. **Reference pressure**: keep `pe` (the first unknown) as the pressure at the end of the
   two-phase zone (dew point), so `Te = T_dew(pe)` and the superheat `ΔTsup` keep their meaning;
   `p_in = pe + Δp_TP`. The superheat zone has its own (single-phase) pressure drop, so the
   compressor suction pressure is `pe − Δp_sup`; `compression` must use it (it uses `pe` now).
4. **Inlet state**: the inlet enthalpy `hei` (isenthalpic expansion) is fixed; the inlet quality
   now depends on `p_in`.
5. **Implicit coupling**: `Δp` depends on the refrigerant flow, quality range and pressure, i.e. on
   the unknowns; evaluate it inside the residual (cheap, AD through the correlation) rather than
   lagging it.
6. **Blends**: `φ` from the zone's own range, `C_r = ṁ (h_dew − h_in)/(T_dew − T_in)` at constant
   pressure (Qiao eq. 3 with the inlet quality instead of the bubble point). Natural refrigerants
   (pure) first; R454C as the blend case (R407C/R410A still blocked by R125 in Clapeyron).

## Refrigerant pressure-drop correlation (needed first)
Brazed plate heat exchanger, flow boiling. Candidates:
- **Amalfi, Vakili-Farahani & Thome (2016)**, Int. J. Refrig. 61, 185–203 (part 2): fitted to a
  large multi-fluid database of flow boiling in plate exchangers; frictional gradient from a
  two-phase multiplier on the all-liquid value, plus the momentum (acceleration) term.
- **Longo et al.**: measured R290 (and other natural refrigerants) boiling in BPHEs; frictional
  pressure drop linear in the kinetic energy per unit volume, `Δp_f ∝ (G²/2ρ_m)`; simple,
  fluid-specific.
- Single phase (superheat zone): the Martin correlation already used for the brine side
  (`f_Martin`), with the vapour properties.

Plus manifold/port losses (≈ 1.5 velocity heads, as in Shah & Focke) and the acceleration term
`G² Δ(x²/ρ_v + (1−x)²/ρ_l…)` (homogeneous or separated-flow void fraction).

## Plan (separate branches, as usual)
1. `feature/refrigerant-pressure-drop`: two-phase (chosen correlation) and single-phase
   refrigerant pressure drop in the plate exchanger, as functions of the zone state; tests against
   published values. No model change.
2. `feature/glide-epsntu`: Qiao eq. 56 in the two-phase zone (φ, γ), `ϵe` eliminated, compressor
   at the suction pressure. Tests: the limits (γ → 0, φ → 0, both), and a **reference
   discretised evaporator** (integrate the local balances along the area with the linear
   pressure profile) that the closed form must match.
3. Study: effect on Te, capacity, COP and on the optimal brine flow (the refrigerant Δp depends
   on the refrigerant flow, i.e. on the capacity, so it couples weakly to the brine flow).

## Open questions
- **Q1** Pressure-drop correlation: Amalfi et al. 2016 (general, recommended) or Longo (R290
  specific), or one you already use? Include port/manifold losses?
- **Q2** Eliminate `ϵe` (recommended) or keep it with a wider bound?
- **Q3** Also the superheat-zone pressure drop and the lower suction pressure in the compressor
  (recommended: yes, it is part of the same effect on COP)?
- **Q4** Plate geometry for the refrigerant side: the `PlateHeatExchanger` has the brine-side
  data (plates, gap `b`, enlargement `ϕ`, chevron `β`); same channels for the refrigerant
  (symmetric BPHE), or separate data?

## Decisions (user, 2026-11-18)
- **Q1**: implement the general correlation (Amalfi et al. 2016) **and** the R290-specific one
  (Longo) and compare them on our operating range; if they agree, use the general one.
- **Q2**: eliminate `ϵe`. The set of unknowns is then **the same with and without** the
  refrigerant pressure drop (γ = 0 gives `ε = 1 − e^{−NTU}` explicitly): the pressure drop is a
  model option of the `Evaporator` (e.g. `refrigerant_pressure_drop = nothing | Amalfi() |
  Longo()`), so the two cases can be compared at the same points, with the same guesses, or by a
  continuation that scales the pressure drop from 0 to 1. A test locks that the option off
  reproduces the current model.
- **Q3**: yes, superheat-zone pressure drop and suction pressure in the compressor.
- **Q4**: same channels for both fluids (brazed plates stacked: symmetric geometry), i.e. the
  `PlateHeatExchanger` data serve both sides.

## Sources found so far
- Desideri et al. 2017 (Int. J. Heat Mass Transfer 113, 6–21; open access, DTU Orbit) gives the
  non-frictional terms used in BPHE data reduction (from Claesson): momentum
  `Δp_mom = G² (1/ρ_tp,out − 1/ρ_tp,in)` with the homogeneous density
  `ρ_tp = (x/ρ_v + (1−x)/ρ_l)⁻¹` (eqs. 16–17); gravity `ρ_m g L_p` (18–19); manifolds and ports
  `Δp_mp = 1.5 G_port²/(2ρ_m)` (20); and the Longo-type frictional form
  `Δp_f = C G²/(2ρ_m)` (linear in the kinetic energy per unit volume, eq. 38–39; their C = 138
  is fitted to ORC fluids at 100–130 °C, **not** Longo's R290 value). It does not give Amalfi's
  equation. Their finding: at high saturation temperature all general correlations (Amalfi
  included) overpredict strongly; Amalfi's refrigerant data cover T_sat −3…35 °C, our evaporation
  is ≈ −6…−14 °C, slightly below: one more reason for the comparison with Longo.
- Still needed (exact equations): Amalfi, Vakili-Farahani & Thome 2016, Int. J. Refrig. 61,
  185–203 (part 2, frictional pressure gradient); Longo et al., R290 (HC290) vaporisation in BPHE
  (frictional constant of the kinetic-energy correlation, and its geometry/range).

## Open
- **Q5** Original API (`evaporator_model!`, `evaporatorborehole_model!`, `original()` in the
  example, the residual comparison tests): they share `evaporator_balances!` with the component.
  Recommendation: leave them as they are (8 unknowns, no pressure drop) as the frozen reference
  of the original formulation, and give the component `Evaporator` its own residual with the
  explicit ε and the pressure-drop option; the "option off" test then compares the component with
  that reference.
