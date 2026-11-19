#= Temperature profiles along the plate evaporator, with and without the refrigerant pressure drop
(ROADMAP step 5, docs/dev/glide.md).

The refrigerant pressure drop lowers its saturation temperature along the two-phase zone: the
refrigerant enters warmer (T_r,in = T_dew(pe + Δp)) than it leaves the zone (Te = T_dew(pe)).
Profiles from the solution of the evaporator + borehole system (`optimalflow_components.jl`):
- two-phase zone (area fraction 0 … xs): refrigerant linear from T_r,in to Te (pressure drop linear
  along the flow, as in Qiao 2025); brine from the local balance
  dD/ds = NTU D − ΔT_sat, D = T_brine − T_r, s = a/xs (closed form, as the ε of Qiao eq. 56/58);
- superheat zone (xs … 1): counter-flow single-phase profiles, exponential in the area, matched to
  the four end temperatures of the zone.

Needs NonlinearSolve and CairoMakie in the environment. Usage:
    include("examples/evaporator_profiles.jl")
    cases = profile_cases(Refrigerant("R290"); Tb = 273.15 + 5, Vf = 1.0)
    fig = plot_profiles(cases); save("examples/figures/evaporator_profiles.png", fig)
=#
include(joinpath(@__DIR__, "optimalflow_components.jl"))
using CairoMakie

"""
    solve_point(r, Vf, Tb; refrigerant_pressure_drop) -> (S, u)

Evaporator + borehole at flow factor `Vf` and ground temperature `Tb` (continuation in the flow
from 2 × design).
"""
function solve_point(r, Vf, Tb; refrigerant_pressure_drop = default_refrigerant_pressure_drop(r))
    build(λ) = heatpump(r, λ; Tb = Tb, refrigerant_pressure_drop = refrigerant_pressure_drop)
    flows = Vf < 2 ? (2.0:-0.05:Vf) : [Vf]
    res = continuation(build, flows, first_guess(build(first(flows)), Tb), TrustRegion(); abstol = 1e-9)
    last(res).ok || error("no solution at $Vf × design flow, Tb = $(Tb - 273.15) °C")
    return build(Vf), last(res).u
end

"""
    profiles(S, u; n = 200) -> NamedTuple

Area fraction `a` (refrigerant inlet at 0) and the refrigerant `Tr` and brine `Ts` temperatures
along the evaporator, plus the summary `dp` of [`evaporator_pressure_drop`](@ref).
"""
function profiles(S, u; n = 200)
    e = S.components.evaporator
    x = physical(S, u).evaporator
    p = port_values(S, u, :evaporator)
    dp = evaporator_pressure_drop(e, x, p)
    xs, Te, Tr_in, N = x.xs, dp.Te, dp.Tr_in, dp.NTU
    ΔTsat = Te - Tr_in
    # two-phase zone
    s = range(0, 1, n)
    Tr2 = Tr_in .+ ΔTsat .* s
    D0 = p.Ts_out - Tr_in
    D = (D0 - ΔTsat/N) .* exp.(N .* s) .+ ΔTsat/N
    Ts2 = Tr2 .+ D
    # superheat zone: D = Ts − Tr exponential in the area, matched to the end temperatures
    Tsup = Te + e.ΔTsup
    Da, Db = x.Ts_m - Te, p.Ts_in - Tsup
    k = log(Db / Da)
    shape(t) = abs(k) < 1e-9 ? t : (exp(k * t) - 1) / (exp(k) - 1)
    Tr3 = Te .+ e.ΔTsup .* shape.(s)
    Ts3 = x.Ts_m .+ (p.Ts_in - x.Ts_m) .* shape.(s)
    a = [xs .* s; xs .+ (1 - xs) .* s]
    return (a = a, Tr = [Tr2; Tr3], Ts = [Ts2; Ts3], xs = xs, dp = dp, Qe = x.Qe, Qsup = x.Qsup)
end

"""
    profile_cases(r; Tb = 278.15, Vf = 1.0) -> Vector

Profiles without refrigerant pressure drop, with Amalfi et al. 2016 (default) and with the
Longo-type fit (R404A, transferred to this plate: an upper estimate, see docs/dev/glide.md).
"""
function profile_cases(r; Tb = 273.15 + 5, Vf = 1.0)
    cases = (("no refrigerant Δp", nothing), ("Amalfi 2016", Amalfi2016()),
             ("Longo-type, ζ = 1.667", LONGO_FITS.r404a_vaporisation))
    return [(label = l, profiles(solve_point(r, Vf, Tb; refrigerant_pressure_drop = m)...)...) for (l, m) in cases]
end

function plot_profiles(cases)
    C = 273.15
    fig = Figure(size = (1000, 470))
    ax1 = Axis(fig[1, 1]; xlabel = "area fraction (refrigerant inlet → outlet)", ylabel = "temperature [°C]",
               title = "Evaporator temperature profiles")
    ax2 = Axis(fig[1, 2]; xlabel = "area fraction", ylabel = "refrigerant temperature [°C]",
               title = "Two-phase zone: refrigerant")
    styles = (:dash, :solid, :dot)
    for (c, ls) in zip(cases, styles)
        lines!(ax1, c.a, c.Ts .- C; color = :steelblue, linestyle = ls)
        lines!(ax1, c.a, c.Tr .- C; color = :firebrick, linestyle = ls, label = c.label)
        i = findall(≤(c.xs), c.a)
        lines!(ax2, c.a[i], c.Tr[i] .- C; color = :firebrick, linestyle = ls,
               label = string(c.label, ": Δp = ", round(c.dp.Δp_twophase/1e3, digits = 1), " kPa"))
        vlines!(ax1, [c.xs]; color = (:gray, 0.4), linestyle = ls)
    end
    text!(ax1, 0.02, 0.97; text = "brine (blue), refrigerant (red)", space = :relative, align = (:left, :top))
    Legend(fig[2, 1], ax1; orientation = :horizontal, framevisible = false)
    axislegend(ax2; position = :rt, framevisible = false)
    return fig
end

function report_profiles(cases)
    for c in cases
        @printf("%-24s Δp = %6.0f Pa  T_r,in = %6.2f °C  Te = %6.2f °C  γ = %.3f  Qe = %5.0f W  xs = %.3f\n", c.label,
                c.dp.Δp_twophase, c.dp.Tr_in - 273.15, c.dp.Te - 273.15, c.dp.γ, c.Qe, c.xs)
    end
end
