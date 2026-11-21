#=
Condenser as a component (step 9c): plate condenser, refrigerant against water in counter-flow.
Two zones along the refrigerant: desuperheating (vapour from the compressor discharge to the dew
point) and condensation (dew point to saturated liquid at the outlet; no subcooling zone, as the
current cycle). ε-NTU in each zone, as the evaporator.
=#

"""
    Condenser(; refrigerant, secondaryfluid, plateheatexchanger, Vs, Qscale,
              condensation = LongoJH())

Two-zone (desuperheating + condensation) moving-boundary plate condenser; the water (secondary
fluid, volume flow `Vs`) flows against the refrigerant. The refrigerant leaves as saturated
liquid.

- Ports: `Tw_in`, `Tw_out` (water), `m` (refrigerant mass flow [kg/s]), `h_in` (refrigerant
  enthalpy at the inlet, the compressor discharge [J/kg]), `pc` (condensing pressure) and `h_out`
  (outlet enthalpy, saturated liquid: to the expansion valve). In a closed cycle `m` and `h_in`
  come from the compressor; standalone, fix them and leave `pc`, `h_out` free. Transforms of the
  refrigerant ports from [`CyclePorts`](@ref).
- Unknowns: `Qd`, `Qc` (duties of the desuperheating and condensation zones,
  `Linear(0, Qscale)`), `xd` (area fraction of the desuperheating zone, `Unit()`), `Tw_m` (water
  at the zone boundary).
- Residuals (W, scaled by `Qscale`): in each zone the refrigerant balance, the water balance and
  the counter-flow ε-NTU relation ([`effectiveness_counterflow`](@ref); desuperheating with the
  capacity ratio of the vapour, condensation with `φ` for the glide of a blend); and
  `h_out − h_bubble(pc)` (scaled by the latent heat).
- U values: condensation [`condensation_htc`](@ref) (`condensation`, default [`LongoJH`](@ref),
  over the zone length) in series with the water (Wanniarachchi); desuperheating: vapour and water
  both with Wanniarachchi. No refrigerant pressure drop.

Guess: `guess(cond; Qd, Qc, Tw_m, xd = 0.2)`; the nodes `pc` (e.g. the dew pressure at the
expected condensing temperature) and `h_out` (`h_bubble(pc)`) are guessed with the other nodes.
Post-processing: [`condenser_zones`](@ref).
"""
struct Condenser{R <: AbstractRefrigerant, P <: PlateHeatExchanger, T <: Real, S <: SecondaryFluid,
                 M <: CondensationModel, X <: CyclePorts} <: AbstractComponent
    refrigerant::R
    secondaryfluid::S
    plateheatexchanger::P
    Vs::T
    Qscale::T
    condensation::M
    cycleports::X
end

function Condenser(; refrigerant, secondaryfluid, plateheatexchanger, Vs, Qscale, condensation::CondensationModel = LongoJH())
    Vs, Qscale = promote(float(Vs), float(Qscale))
    return Condenser(refrigerant, secondaryfluid, plateheatexchanger, Vs, Qscale, condensation, CyclePorts(refrigerant))
end

function unknowns(c::Condenser)
    Q = Linear(0, c.Qscale)
    return (Qd = Q, Qc = Q, xd = Unit(), Tw_m = Free())
end
function ports(c::Condenser)
    cp = c.cycleports
    return (Tw_in = Free(), Tw_out = Free(), m = cp.m, h_in = cp.h, pc = cp.p, h_out = cp.h)
end
residual_scales(c::Condenser) = (ntuple(_ -> c.Qscale, 6)..., c.cycleports.hscale)

guess(c::Condenser; Qd, Qc, Tw_m, xd = 0.2) = (Qd = Qd, Qc = Qc, xd = xd, Tw_m = Tw_m)

"""
    vapour_T_ph(r, p, h) -> T

Temperature of superheated vapour at pressure `p` and enthalpy `h` (≥ the dew enthalpy): Newton in
T on `h_pT(…; phase = :vapor)` from the dew point (no flash; four steps, quadratic convergence).
"""
function vapour_T_ph(r::AbstractRefrigerant, p, h)
    Td = T_dew(r, p)
    T = Td + (h - h_dew(r, p)) / cp_dew(r, p)
    for _ in 1:4
        T -= (h_pT(r, p, T; phase = :vapor) - h) / cp_pT(r, p, T; phase = :vapor)
        T = max(T, Td)
    end
    return T
end

# zone quantities shared by the residuals and the post-processing
function _condenser_zones(c::Condenser, x, p)
    r, sf, phx = c.refrigerant, c.secondaryfluid, c.plateheatexchanger
    pc, m = p.pc, p.m
    Td, Tb = T_dew(r, pc), T_bubble(r, pc)
    hd, hb = h_dew(r, pc), h_bubble(r, pc)
    Tdis = vapour_T_ph(r, pc, p.h_in)
    wprops(T) = (sf.ρ(T), sf.μ(T), sf.λ(T), sf.cp(T))
    # condensation zone: water from Tw_in to Tw_m
    Twc = (p.Tw_in + x.Tw_m) / 2
    Cwc = sf.ρ(Twc) * sf.cp(Twc) * c.Vs
    hwc = h_platehex_wanniarachchi(c.Vs, phx, wprops(Twc))
    Ac = phx.A * (1 - x.xd)
    q = max(abs(x.Qc), 1e-3 * c.Qscale) / Ac
    hrc = condensation_htc(c.condensation, r, phx, m, pc, one(q), zero(q), q; L = (1 - x.xd) * phx.L).h
    Uc = inv(inv(hwc) + inv(hrc))
    NTUc = Uc * Ac / Cwc
    φc = Cwc * (Td - Tb) / (m * (hd - hb))                  # glide of a blend (0 for a pure fluid)
    εc = effectiveness_counterflow(NTUc, zero(NTUc), φc)
    # desuperheating zone: water from Tw_m to Tw_out, vapour from Tdis to Td
    Twd = (x.Tw_m + p.Tw_out) / 2
    Cwd = sf.ρ(Twd) * sf.cp(Twd) * c.Vs
    hwd = h_platehex_wanniarachchi(c.Vs, phx, wprops(Twd))
    Tv = (Tdis + Td) / 2
    ρv = ρ_pT(r, pc, Tv; phase = :vapor)
    hrd = h_platehex_wanniarachchi(m / ρv, phx, (ρv, μ_pT(r, pc, Tv; phase = :vapor), λ_pT(r, pc, Tv; phase = :vapor),
                                                cp_pT(r, pc, Tv; phase = :vapor)))
    Ud = inv(inv(hwd) + inv(hrd))
    NTUd = Ud * phx.A * x.xd / Cwd
    Crd = m * (p.h_in - hd) / max(Tdis - Td, 1e-6)        # mean vapour capacity
    εd = effectiveness_counterflow(NTUd, zero(NTUd), Cwd / Crd)
    return (Td = Td, Tb = Tb, hd = hd, hb = hb, Tdis = Tdis,
            condensation = (U = Uc, h_ref = hrc, h_water = hwc, Cw = Cwc, NTU = NTUc, φ = φc, ε = εc),
            desuperheat = (U = Ud, h_ref = hrd, h_water = hwd, Cw = Cwd, Cref = Crd, NTU = NTUd, ε = εd))
end

function residuals!(r, c::Condenser, x, p)
    z = _condenser_zones(c, x, p)
    zc, zd = z.condensation, z.desuperheat
    r[1] = p.m * (p.h_in - z.hd) - x.Qd                    # refrigerant, desuperheating
    r[2] = p.m * (z.hd - z.hb) - x.Qc                      # refrigerant, condensation
    r[3] = zd.Cw * (p.Tw_out - x.Tw_m) - x.Qd              # water, desuperheating zone
    r[4] = zc.Cw * (x.Tw_m - p.Tw_in) - x.Qc               # water, condensation zone
    r[5] = zc.ε * zc.Cw * (z.Td - p.Tw_in) - x.Qc          # ε-NTU, condensation (refrigerant enters at the dew point)
    r[6] = zd.ε * zd.Cw * (z.Tdis - x.Tw_m) - x.Qd         # ε-NTU, desuperheating
    r[7] = p.h_out - z.hb                                   # saturated liquid at the outlet
    return nothing
end

"""
    condenser_zones(c::Condenser, x, p) -> NamedTuple

Zone quantities of a solution (`x` physical unknowns, `p` port values): dew and bubble
temperatures `Td`, `Tb` and enthalpies at the condensing pressure, discharge temperature `Tdis`,
and per zone (`condensation`, `desuperheat`) the U value, the refrigerant- and water-side
coefficients, capacities, NTU and effectiveness. The refrigerant leaves at `hb` (saturated liquid).
"""
condenser_zones(c::Condenser, x, p) = _condenser_zones(c, x, p)

"""
    condenser_profiles(c::Condenser, x, p; n = 100) -> NamedTuple

Temperature profiles along the condenser from a solution (`x` physical unknowns, `p` port values):
area fraction `a` from the refrigerant inlet (0) to the outlet (1), refrigerant `Tr` and water `Tw`,
`n` points per zone, the zone boundary `xd`, and `end_error`: the largest difference [K] between
the far ends of the integrated profiles and the solution (`Td`, `Tw_m`, `Tw_in`).
- Desuperheating (0 … xd): counter-flow single-phase streams with the zone's U value and capacities
  (`condenser_zones`): with D = Tr − Tw along the refrigerant flow, dTr/dt = −NTU_r D and
  dTw/dt = −NTU_w D (the water leaves at t = 0), from Tdis and Tw_out.
- Condensation (xd … 1): refrigerant linear in the heat from Td to Tb (the glide of a blend), water
  from the balance; the heat fraction z obeys dz/dt = NTU (Cw D₀/Qc + (1 − φ) z), from Td and Tw_m.
"""
function condenser_profiles(c::Condenser, x, p; n::Integer = 100)
    z = _condenser_zones(c, x, p)
    s = range(0, 1, n)
    # desuperheating zone
    d = z.desuperheat
    NTUw = d.NTU
    NTUr = d.NTU * d.Cw / d.Cref
    k = NTUw - NTUr
    D0 = z.Tdis - p.Tw_out
    I(t) = abs(k) < 1e-12 ? D0 * t : D0 * expm1(k * t) / k                  # ∫₀ᵗ D
    Tr1 = z.Tdis .- NTUr .* I.(s)
    Tw1 = p.Tw_out .- NTUw .* I.(s)
    # condensation zone
    cz = z.condensation
    ΔTg = z.Td - z.Tb
    b = cz.NTU * (1 - cz.φ)
    a0 = cz.NTU * cz.Cw * (z.Td - x.Tw_m) / x.Qc
    Z = [abs(b) < 1e-12 ? a0 * t : a0 * expm1(b * t) / b for t in s]
    Tr2 = z.Td .- ΔTg .* Z
    Tw2 = x.Tw_m .- (x.Qc / cz.Cw) .* Z
    end_error = maximum(abs, (Tr1[end] - z.Td, Tw1[end] - x.Tw_m, Tw2[end] - p.Tw_in, (Z[end] - 1) * (z.Td - x.Tw_m)))
    return (a = [x.xd .* s; x.xd .+ (1 - x.xd) .* s], Tr = [Tr1; Tr2], Tw = [Tw1; Tw2], xd = x.xd, end_error = end_error)
end

"""
    validity(c::Condenser, x, p; ΔTpinch = 0.1, xd_tol = 1e-3)

Issues of a condenser solution: water outside its property table; pinch (within `ΔTpinch`) at the
refrigerant inlet (`Tdis − Tw_out`), at the zone boundary (`Td − Tw_m`) and at the refrigerant
outlet (`Tb − Tw_in`); a zone almost vanished (`xd`).
"""
function validity(c::Condenser, x, p; ΔTpinch = 0.1, xd_tol = 1e-3)
    issues = NamedTuple[]
    _secondary_issues!(issues, c.secondaryfluid, (Tw_in = p.Tw_in, Tw_m = x.Tw_m, Tw_out = p.Tw_out))
    z = _condenser_zones(c, x, p)
    z.Tdis - p.Tw_out < ΔTpinch && push!(issues, _issue(:Tdis, z.Tdis, p.Tw_out + ΔTpinch, "condenser pinch at the refrigerant inlet: Tdis - Tw_out below ΔTpinch"))
    z.Td - x.Tw_m < ΔTpinch && push!(issues, _issue(:Td, z.Td, x.Tw_m + ΔTpinch, "condenser pinch at the zone boundary: Td - Tw_m below ΔTpinch"))
    z.Tb - p.Tw_in < ΔTpinch && push!(issues, _issue(:Tb, z.Tb, p.Tw_in + ΔTpinch, "condenser pinch at the refrigerant outlet: Tb - Tw_in below ΔTpinch"))
    x.xd < xd_tol && push!(issues, _issue(:xd, x.xd, xd_tol, "desuperheating zone almost vanished"))
    x.xd > 1 - xd_tol && push!(issues, _issue(:xd, x.xd, 1 - xd_tol, "condensation zone almost vanished"))
    return issues
end

"""
    condenser_design(Q, Tc, Tw_in, Tw_out, water::SecondaryFluid; U = 2000.) -> (; A, Vs)

Basic sizing, as [`EvaporatorDesign`](@ref): heat-transfer area `A = Q/(U ΔT_lm)` for condensation at
`Tc` against water from `Tw_in` to `Tw_out` (desuperheating neglected), and the water volume flow
`Vs = Q/(ρ cp (Tw_out − Tw_in))`.
"""
function condenser_design(Q, Tc, Tw_in, Tw_out, water::SecondaryFluid; U = 2000.)
    ΔTlm = (Tw_out - Tw_in) / log((Tc - Tw_in) / (Tc - Tw_out))
    Tm = (Tw_in + Tw_out) / 2
    return (A = Q / (U * ΔTlm), Vs = Q / (water.ρ(Tm) * water.cp(Tm) * (Tw_out - Tw_in)))
end
