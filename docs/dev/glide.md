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
- **Amalfi, Vakili-Farahani & Thome 2016** (Int. J. Refrig. 61; parts 1 and 2 in `assets/docs`):
  - frictional (part 1 eq. 18): `Δp_f = 2 f_tp L_p G²/(d_h ρ_m)`, `ρ_m` homogeneous at the mean
    quality `x_m` (part 1 eq. 17: `1/ρ_m = x_m/ρ_v + (1−x_m)/ρ_l`), `d_h = 2b/ϕ` (part 1 eq. 6,
    `b` the gap = 2 × corrugation amplitude: our `PlateHeatExchanger.b`);
  - `f_tp = C · 15.698 · We_m^−0.475 · Bd^0.255 · ρ*^−0.571` (part 2 eq. 7),
    `C = 2.125 β*^9.993 + 0.955` (eq. 8), `β* = β/70°`, β the chevron angle from the flow
    (vertical) axis, as Martin's angle in `f_Martin` (our `PlateHeatExchanger.β`);
    `We_m = G² d_h/(ρ_m σ)`, `Bd = (ρ_l − ρ_v) g d_h²/σ`, `ρ* = ρ_l/ρ_v` (part 2 eq. 5);
  - fit: 1513 points, MAE 21.5 %, 90.9 % within ±50 %; refrigerants T_sat −3…35 °C;
  - other terms (part 1 eqs. 8–11): `Δp = Δp_g + Δp_acc + Δp_f + Δp_p`, `Δp_g = ρ_m g L_p`,
    `Δp_acc = G² Δx (1/ρ_v − 1/ρ_l)`, ports and manifolds (Shah & Focke)
    `Δp_p = 0.75 N_pass [(G_p²/2ρ)_in + (G_p²/2ρ)_out]`, `G_p` the port mass flux.
  - Needs the surface tension σ (not in the refrigerant interface yet: CoolProp `I` through the
    cached state, which also handles Duals by finite differences).
- **Longo et al. 2016** (Appl. Therm. Eng. 108, 1401–1408): same BPHE family as their R290 work
  (L = 278 mm, W = 72 mm, chevron 65°, corrugation depth b = 2 mm, pitch 8 mm, 10 plates,
  `d_h = 2b` without enlargement factor); `Δp_f [kPa] = 1.667 KE/V [J m⁻³]`, `KE/V = G²/(2ρ_m)`
  (eqs. 22–23), MAPD 10.8 %, **fitted to HFC404A at 20 °C**; HC290 and HC1270 lie above HFC404A
  at the same mass flux (Fig. 11) but their own fit is in **Longo 2012** (ASME J. Heat Transfer
  134, 101801, "Hydrocarbon refrigerant vaporisation inside a brazed plate heat exchanger"),
  not in the folder. Longo et al. 2023 (local HTC of R290/R1270) refers to Longo 2012 for the
  pressure drop. To transfer a KE/V fit to another plate: `Δp_f = ζ KE/V` with `ζ = 4 f_eq L/d_h`,
  i.e. an equivalent constant Fanning factor `f_eq = ζ_ref d_h,ref/(4 L_ref)` (HFC404A:
  ζ = 1667, f_eq ≈ 6.0): an assumption to state.
- **Q1 comparison status**: Amalfi can be implemented now; the R290-specific comparison needs
  Longo 2012 (or the Longo et al. 2015 computational procedure, IJHMT, if it gives Δp).

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

## Q1 comparison: Amalfi vs Longo-type (2026-11-18, `feature/refrigerant-pressure-drop`)
More Longo fits on the same BPHE (L 278 mm, W 72 mm, 65°, b 2 mm), `Δp_f [kPa] = ζ KE/V [J m⁻³]`:
ζ = 1.49 (HFC vaporisation, Longo & Gasparella 2007), 1.667 (R404A vaporisation, Longo et al.
2016), 1.90 (R600a/R290/R1270 **condensation**, Longo 2010). R290 vaporisation (Longo 2012) is not
available: the band 1.49–1.90 is used as the Longo-type estimate.

On Longo's plate (ϕ = 1.14 from the corrugation, x 0.2 → 1), frictional Δp [kPa]:

| T_sat | G [kg/(m² s)] | Amalfi | ζ = 1.49 … 1.90 |
|---|---|---|---|
| 10 °C | 10 | 3.4 | 3.3 … 4.2 |
| 10 °C | 30 | 10.9 | 29.7 … 37.9 |
| 20 °C | 10 | 3.4 | 2.5 … 3.2 |
| 20 °C | 20 | 7.0 | 10.1 … 12.9 |

Longo 2016 Fig. 11 (R290 data, 20 °C): ≈ 2.4 kPa at G = 10, ≈ 8.5 kPa at G = 20 (read from the
plot), between the two. Amalfi grows ≈ G^1.05 (the Weber exponent), the Longo form as G².

Our evaporator (example, design flow, Tb = 5 °C): Te = −4.6 °C, ṁ = 19.2 g/s,
**G ≈ 10 kg/(m² s)**, x 0.31 → 1 over 87 % of the plate length:

| | two-phase Δp | ΔT_sat |
|---|---|---|
| Amalfi | 2.9 kPa | −0.22 K |
| Longo-type ζ 1.49 … 1.90 (scaled with L/2b) | 8.4 … 10.7 kPa | −0.64 … −0.82 K |

Acceleration 7 Pa, gravity 58 Pa, superheat zone (vapour, Martin) 15 Pa, ports 13 Pa: negligible.
Most of the factor ≈ 3 comes from transferring Longo's fit to our plate: it ignores the chevron
angle (65° there, 30° here); Amalfi's β-dependence alone is a factor 2.07 between the two angles
(C = 1.98 vs 0.955), which would bring the Longo-type estimate to ≈ 0.3–0.4 K.
Recommendation: Amalfi as the default (general, contains the chevron angle), the Longo-type fits
as a sensitivity bound. Either way the effect is a few tenths of a kelvin at design flow, i.e.
γ = ΔT_sat/(T_r,in − T_f,in) ≈ 0.05–0.2: small but not negligible against a 3–5 K approach.

## Step 5b: the model (2026-11-19, `feature/epsntu-pressure-drop`)
- Component `Evaporator` only (the original API `evaporator_model!` is the frozen reference, Q5):
  unknowns `(pe|Te, Qe, Qsup, ϵsup, xs, Ts_m, …)`: `ϵe` eliminated, 7 residuals.
- Two-phase zone: `Qe = ε C_s (T_s,m − T_r,in)`, `ε = effectiveness_counterflow(NTU, γ)` (Qiao
  eq. 56, φ = 0), `T_r,in = T_dew(pe + Δp)`, `Δp = twophase_pressure_drop(model, …; L = xs L)`
  (friction + acceleration + gravity, upward flow), inlet quality at `pe`.
  `refrigerant_pressure_drop = nothing`: `T_r,in = Te`, γ = 0: exactly the original equations
  (test: component residuals = original residuals with `ϵe = 1 − e^{−NTU}`).
- Default `default_refrigerant_pressure_drop(r)`: `Amalfi2016()` for pure fluids, `nothing` for
  mixtures (type-determined, inferred); a mixture with a pressure-drop model is refused.
- Post-processing: `evaporator_pressure_drop(e, x, p)` → Δp, `T_r,in`, Te, γ, NTU, superheat-zone
  Δp; `compressor_power(…; Δp_suction)` (also condenser and `compressor_condenser_power`).
- Validation: Qiao's closed form against an RK4 integration of the zone (1e-8); limits.
- Design point (R290, Tb = 5 °C, design flow; Clapeyron): Δp = 3.2 kPa, `T_r,in − Te` = 0.21 K,
  γ = 0.037, NTU = 0.82; superheat-zone Δp 29 Pa.
- Pitfall found (tests): a closure (`map(...) do`) calling a local keyword function of the
  testset (`coupled`) made the testset compile for 25 min (abstract NonlinearSolve/Clapeyron
  inference); straight-line code: 39 s.
- Example `examples/evaporator_profiles.jl` (spline R290, Tb = 5 °C, design flow; figure in
  `examples/figures/evaporator_profiles.png`, not tracked):

  | case | two-phase Δp | T_r,in | Te | γ | Qe |
  |---|---|---|---|---|---|
  | no refrigerant Δp | 0 | −4.60 °C | −4.60 °C | 0 | 5062 W |
  | Amalfi 2016 (default) | 2.9 kPa | −4.47 °C | −4.69 °C | 0.046 | 5047 W |
  | Longo-type ζ = 1.667 | 9.5 kPa | −4.16 °C | −4.89 °C | 0.160 | 5012 W |

  The refrigerant enters the two-phase zone warmer and leaves it colder; the brine profiles hardly
  change. Te (and so the compressor suction) falls by 0.09 K (Amalfi) to 0.29 K (Longo-type), the
  two-phase duty by 0.3–1 %.
- Variable speed (`speed_sweep` in the same example; figure `examples/figures/evaporator_speed.png`):
  compressor speed (displaced volume flow `Vref`) 0.3 → 1.5 × design at design brine flow,
  Tb = 5 °C; suction pressure lowered by the superheat-zone Δp.

  | speed | G [kg/(m² s)] | Amalfi Δp / T_r,in − Te | Longo-type Δp / T_r,in − Te | COP1 change Amalfi / Longo-type |
  |---|---|---|---|---|
  | 0.3 | 3.4 | 0.4 kPa / 0.03 K | 0.4 kPa / 0.03 K | ≈ 0 / ≈ 0 |
  | 0.5 | 5.5 | 1.3 kPa / 0.09 K | 2.1 kPa / 0.15 K | −0.05 % / −0.07 % |
  | 1.0 | 9.8 | 2.9 kPa / 0.22 K | 9.5 kPa / 0.73 K | −0.2 % / −0.7 % |
  | 1.5 | 13.4 | 4.2 kPa / 0.35 K | 20.1 kPa / 1.65 K | −0.3 % / −1.4 % |

  The effect grows with speed: Amalfi ≈ linearly in G (`f_tp ∝ We^−0.475`), the Longo form as
  G², so the two models differ most at high speed (factor 4.8 at 1.5 ×). At part load (≤ 0.5 ×)
  the effect is negligible in both. The mass fluxes of this evaporator (3–13 kg/(m² s)) are at the
  low end of Amalfi's database (5.5–610): the Longo 2012 R290 data would settle which one holds.
  Without refrigerant Δp the COP1 maximum is at ≈ 0.5 × speed (4.097, part-load benefit).
- Profiles from the model (5c): `evaporator_zones(e, x, p)` gives U, capacities and NTU of each
  zone; the superheat profile is integrated with `NTU_s`, `NTU_r` of the model (before: an
  exponential fitted to the end temperatures). In both zones the far-end temperatures are
  reproduced to 6e-14 K (test in `test/components.jl`).
- Cost of the refrigerant pressure drop (compiled, per call; design point, Tb = 5 °C):

  | | residual | warm Newton | cold TrustRegion |
  |---|---|---|---|
  | spline, no Δp / Amalfi | 1.2 / 4.5 µs | 0.04 / 0.21 ms | 0.15 / 0.82 ms (7 / 8 it) |
  | Clapeyron, no Δp / Amalfi | 10.8 / 26.0 µs | 0.56 / 0.89 ms | 2.15 / 3.42 ms (7 / 8 it) |

  The extra cost is property calls in the residual: ρ_bubble, h_bubble and σ (three separate
  CoolProp updates; with Duals, three each for the finite differences) and `T_dew(pe + Δp)`
  (an inversion for Clapeyron). Possible savings if it matters: one CoolProp update for the
  three saturated-liquid properties (`_coolprop_values`), and `T_dew(pe + Δp) ≈ Te + (dT/dp) Δp`.
