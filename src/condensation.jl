#=
Refrigerant-side heat transfer during condensation in a brazed plate heat exchanger (step 9b).

Longo 2010 (Int. J. Refrig. 33, 944, `assets/docs/`): R290, R600a and R1270 condensing in a BPHE
(65° chevron, b = 2 mm). Two regimes, with the transition at an equivalent Reynolds number of
about 1600–1700 (G ≈ 15–18 kg/(m² s)):
- gravity controlled: Nusselt (1916) film condensation on a vertical surface (14.8 % mean
  deviation from the data);
- forced convection: Akers et al. (1959), local in the vapour quality (21 %).
Longo refers h and q to the projected plate area (and multiplies by the enlargement factor Φ);
here both are referred to the actual area `A` of the `PlateHeatExchanger`, as in the evaporator,
so Φ drops out: h_projected = Φ h, q_projected = Φ q.
=#

"""
    CondensationModel

Refrigerant-side condensation heat-transfer model: [`Longo2010`](@ref).
"""
abstract type CondensationModel end

"""
    Longo2010(; Re_low = 1500.0, Re_high = 1800.0)

Condensation in a brazed plate heat exchanger after Longo 2010: Nusselt film condensation
([`nusselt_film_htc`](@ref), gravity controlled) below the transition and Akers forced-convection
condensation ([`akers_htc`](@ref)) above it, with Longo's criterion: the equivalent Reynolds number
at the mean quality, transition at 1600–1700 for all his fluids. The two are joined by a C¹ blend
(cubic Hermite weight) between `Re_low` and `Re_high`, as the solver needs a differentiable model.
(A smooth maximum of the two is not used: below the transition it takes Akers, which over-predicts
there, and adds the two where they coincide, up to +40 % against Longo's R290 data.)
"""
struct Longo2010 <: CondensationModel
    Re_low::Float64
    Re_high::Float64
end
Longo2010(; Re_low = 1500.0, Re_high = 1800.0) = Longo2010(float(Re_low), float(Re_high))

"""
    nusselt_film_htc(λl, ρl, μl, Δh_lv, L, q)

Nusselt (1916) laminar film condensation on a vertical surface of height `L` [m], mean
heat-transfer coefficient [W/(m² K)] written with the heat flux `q` [W/m²] instead of the wall
superheat ΔT = q/h: `h = 0.943 [λ³ρ²gΔh_lv/(μ L ΔT)]^¼` becomes `h = [0.943⁴ λ³ρ²gΔh_lv/(μ L q)]^⅓`
(liquid conductivity `λl`, density `ρl`, viscosity `μl`; latent heat `Δh_lv`). Falls with the
heat flux (thicker film).
"""
nusselt_film_htc(λl, ρl, μl, Δh_lv, L, q) = (0.943^4 * λl^3 * ρl^2 * GRAVITY * Δh_lv / (μl * L * q))^(1/3)

"""
    akers_htc(G, x_in, x_out, ρl, ρv, μl, λl, cpl, dh)

Akers et al. (1959) forced-convection condensation, `h = 5.03 (λl/dh) Re_eq^⅓ Pr_l^⅓` with
`Re_eq = G [(1 − x) + x (ρl/ρv)^½] dh/μl`, averaged over the vapour quality from `x_in` to
`x_out` (uniform heat flux: x linear in the area; Re_eq linear in x, closed form).
"""
function akers_htc(G, x_in, x_out, ρl, ρv, μl, λl, cpl, dh)
    a = G * dh / μl                           # Re_eq = a + b x
    b = a * (sqrt(ρl / ρv) - 1)
    Δx = x_in - x_out
    Re13 = abs(b * Δx) < 1e-8 * a ? cbrt(a + b * (x_in + x_out) / 2) :
           3 * ((a + b * x_in)^(4/3) - (a + b * x_out)^(4/3)) / (4b * Δx)
    return 5.03 * λl / dh * Re13 * cbrt(cpl * μl / λl)
end

"Equivalent Reynolds number of Akers/Longo at the vapour quality `x`."
reynolds_equivalent(G, x, ρl, ρv, μl, dh) = G * ((1 - x) + x * sqrt(ρl / ρv)) * dh / μl

"""
    condensation_htc(model::Longo2010, r, phx, m, p, x_in, x_out, q; L = phx.L) -> NamedTuple

Mean refrigerant-side heat-transfer coefficient [W/(m² K), actual area] of the refrigerant mass
flow `m` condensing at the pressure `p` from the quality `x_in` to `x_out` with the heat flux `q`
[W/m², actual area] over the flow length `L` of the zone, in the plate exchanger `phx` (`d_h = 2b`,
as Longo). Returns `(; h, gravity, forced, Re_eq)` with `Re_eq` at the mean quality.

Zeotropic blend: Silver–Bell–Ghaly, `1/h = 1/h_film + Z_g/h_g` with `Z_g = x cp_g dT/dh` (the
sensible cooling of the vapour as a fraction of the total heat, `dT/dh = (T_dew − T_bubble)/Δh_lv`)
and `h_g` the vapour flowing alone (Wanniarachchi), at the mean quality.
"""
function condensation_htc(model::Longo2010, r::AbstractRefrigerant, phx::PlateHeatExchanger, m, p, x_in, x_out, q;
                          L = phx.L)
    ρl, ρv = ρ_bubble(r, p), ρ_dew(r, p)
    μl, λl, cpl = μ_bubble(r, p), λ_bubble(r, p), cp_bubble(r, p)
    Δh_lv = h_dew(r, p) - h_bubble(r, p)
    G = channel_mass_flux(phx, m)
    dh = 2phx.b
    hg = nusselt_film_htc(λl, ρl, μl, Δh_lv, L, q)
    hf = akers_htc(G, x_in, x_out, ρl, ρv, μl, λl, cpl, dh)
    xm = (x_in + x_out) / 2
    Re = reynolds_equivalent(G, xm, ρl, ρv, μl, dh)
    t = clamp((Re - model.Re_low) / (model.Re_high - model.Re_low), zero(Re), one(Re))
    w = t^2 * (3 - 2t)                       # C¹ weight of the forced-convection regime
    h = _sbg_correction(r, phx, m, p, xm, (1 - w) * hg + w * hf, ρv, Δh_lv)
    return (h = h, gravity = hg, forced = hf, Re_eq = Re)
end

# Silver–Bell–Ghaly for blends (pure fluids: no correction, no cost)
_sbg_correction(::AbstractRefrigerant, phx, m, p, xm, h, ρv, Δh_lv) = h
function _sbg_correction(r::ClapeyronRefrigerant, phx, m, p, xm, h, ρv, Δh_lv)
    ispure(r) && return h
    cpv = cp_dew(r, p)
    Zg = xm * cpv * (T_dew(r, p) - T_bubble(r, p)) / Δh_lv
    hv = h_platehex_wanniarachchi(xm * m / ρv, phx, (ρv, μ_dew(r, p), λ_dew(r, p), cpv))
    return inv(inv(h) + Zg / hv)
end
