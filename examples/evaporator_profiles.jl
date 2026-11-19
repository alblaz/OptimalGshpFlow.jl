#= Temperature profiles along the plate evaporator, with and without the refrigerant pressure drop
(ROADMAP step 5, docs/dev/glide.md).

The refrigerant pressure drop lowers its saturation temperature along the two-phase zone: the
refrigerant enters warmer (T_r,in = T_dew(pe + Δp)) than it leaves the zone (Te = T_dew(pe)).
Profiles from the solution of the evaporator + borehole system (`optimalflow_components.jl`):
- two-phase zone (area fraction 0 … xs): refrigerant linear from T_r,in to Te (pressure drop linear
  along the flow, as in Qiao 2025); brine from the local balance
  dD/ds = NTU D − ΔT_sat, D = T_brine − T_r, s = a/xs (closed form, as the ε of Qiao eq. 56/58);
- superheat zone (xs … 1): counter-flow single-phase streams with the model's U value and
  capacities (`evaporator_zones`): dT_r/dt = NTU_r D, dT_s/dt = NTU_s D, so D = D(0) e^{(NTU_s − NTU_r) t},
  t = (a − xs)/(1 − xs), from the boundary values Te and Ts_m.
In both zones the profile is integrated from one end of the zone; the far end (Ts_m, and Te + ΔTsup,
Ts_in) is not imposed and is reproduced by the solution of the model (`end_error`).

`speed_sweep` varies the compressor speed (displaced volume flow, i.e. the refrigerant flow, of a
variable-speed heat pump) at a fixed brine flow: pressure drop, refrigerant temperatures and COP.

Needs NonlinearSolve and CairoMakie in the environment. Usage:
    include("examples/evaporator_profiles.jl")
    r = Refrigerant("R290")
    cases = profile_cases(r; Tb = 273.15 + 5, Vf = 1.0)
    fig = plot_profiles(cases); save("examples/figures/evaporator_profiles.png", fig)
    sw = speed_sweep(r); report_speed(sw)
    save("examples/figures/evaporator_speed.png", plot_speed(sw))
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
along the evaporator, the summary `dp` of [`evaporator_pressure_drop`](@ref), and `end_error`: the
largest difference [K] between the far-end temperatures of the integrated profiles and the solution.
"""
function profiles(S, u; n = 200)
    e = S.components.evaporator
    x = physical(S, u).evaporator
    p = port_values(S, u, :evaporator)
    dp = evaporator_pressure_drop(e, x, p)
    zn = evaporator_zones(e, x, p)
    xs, Te, Tr_in, N = x.xs, zn.twophase.Te, zn.twophase.Tr_in, zn.twophase.NTU
    ΔTsat = Te - Tr_in
    # two-phase zone
    s = range(0, 1, n)
    Tr2 = Tr_in .+ ΔTsat .* s
    D0 = p.Ts_out - Tr_in
    D = (D0 - ΔTsat/N) .* exp.(N .* s) .+ ΔTsat/N
    Ts2 = Tr2 .+ D
    # superheat zone: model U value and capacities, integrated from the boundary (Te, Ts_m)
    (; NTU_s, NTU_r) = zn.superheat
    k = NTU_s - NTU_r
    D0 = x.Ts_m - Te
    I(t) = abs(k) < 1e-12 ? D0 * t : D0 * (exp(k * t) - 1) / k       # ∫₀ᵗ D
    Tr3 = Te .+ NTU_r .* I.(s)
    Ts3 = x.Ts_m .+ NTU_s .* I.(s)
    a = [xs .* s; xs .+ (1 - xs) .* s]
    end_error = maximum(abs, (Ts2[end] - x.Ts_m, Tr3[end] - (Te + e.ΔTsup), Ts3[end] - p.Ts_in))
    return (a = a, Tr = [Tr2; Tr3], Ts = [Ts2; Ts3], xs = xs, dp = dp, zones = zn, end_error = end_error,
            Qe = x.Qe, Qsup = x.Qsup)
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
        @printf("%-24s Δp = %6.0f Pa  T_r,in = %6.2f °C  Te = %6.2f °C  γ = %.3f  Qe = %5.0f W  xs = %.3f  end error %.1e K\n",
                c.label, c.dp.Δp_twophase, c.dp.Tr_in - 273.15, c.dp.Te - 273.15, c.dp.γ, c.Qe, c.xs, c.end_error)
    end
end

const PRESSURE_DROP_CASES = (("no refrigerant Δp", nothing), ("Amalfi 2016", Amalfi2016()),
                             ("Longo-type, ζ = 1.667", LONGO_FITS.r404a_vaporisation))

# one operating point of the speed sweep
function speed_point(r, S, u, sf = BRINE_PROPS)
    e = S.components.evaporator
    ph = physical(S, u)
    pv = port_values(S, u, :evaporator)
    dp = evaporator_pressure_drop(e, ph.evaporator, pv)
    es = evaporator_solution(e, ph.evaporator, pv)
    bs = borehole_solution(S.components.borehole, port_values(S, u, :borehole))
    m = e.Vref * ρ_superheated(e.refrigerant, OptimalGshpFlow._pe(e, ph.evaporator), e.ΔTsup)
    perf = performance(es, bs, r, sf; Δp_suction = dp.Δp_superheat)
    return (G = channel_mass_flux(e.plateheatexchanger, m), m = m, Δp = dp.Δp_twophase, Tr_in = dp.Tr_in, Te = dp.Te,
            γ = dp.γ, Δp_sup = dp.Δp_superheat, perf...)
end

"""
    speed_sweep(r; speeds = 0.3:0.1:1.5, Tb = 278.15, Vf = 1.0) -> Vector

For each refrigerant pressure-drop case (none, Amalfi, Longo-type), continuation in the compressor
speed (from 1 down to `minimum(speeds)` and up to `maximum(speeds)`) at the brine flow `Vf` ×
design and ground temperature `Tb`. One entry per case: `(label, rows)`, rows sorted by speed.
"""
function speed_sweep(r; speeds = 0.3:0.1:1.5, Tb = 273.15 + 5, Vf = 1.0)
    out = []
    for (label, model) in PRESSURE_DROP_CASES
        build(s) = heatpump(r, Vf; Tb = Tb, speed = s, refrigerant_pressure_drop = model)
        S1, u1 = solve_point(r, Vf, Tb; refrigerant_pressure_drop = model)      # design speed
        down = sort(filter(≤(1.0), collect(speeds)); rev = true)
        up = sort(filter(≥(1.0), collect(speeds)))
        rows = []
        for branch in (down, up)
            isempty(branch) && continue
            res = continuation(build, branch, u1, TrustRegion(); abstol = 1e-9)
            for x in res
                x.ok || (@warn "$label: no solution at speed $(x.λ)"; continue)
                push!(rows, (speed = x.λ, speed_point(r, build(x.λ), x.u)...))
            end
        end
        rows = unique(x -> x.speed, sort(rows; by = x -> x.speed))
        push!(out, (label = label, rows = rows))
    end
    return out
end

function report_speed(sw)
    for c in sw
        println("== ", c.label)
        @printf("%6s %8s %8s %8s %8s %8s %7s %7s %7s\n", "speed", "G", "Δp[kPa]", "ΔTsat[K]", "Te[°C]", "Q2[W]", "W[W]",
                "COP1", "γ")
        for x in c.rows
            @printf("%6.2f %8.1f %8.2f %8.3f %8.2f %8.0f %7.0f %7.3f %7.3f\n", x.speed, x.G, x.Δp/1e3, x.Tr_in - x.Te,
                    x.Te - 273.15, x.Q2, x.W, x.COP1, x.γ)
        end
    end
end

function plot_speed(sw)
    fig = Figure(size = (1100, 380))
    ax1 = Axis(fig[1, 1]; xlabel = "compressor speed (× design)", ylabel = "two-phase Δp [kPa]",
               title = "Refrigerant pressure drop")
    ax2 = Axis(fig[1, 2]; xlabel = "compressor speed (× design)", ylabel = "T_r,in − Te [K]",
               title = "Saturation temperature drop")
    ax3 = Axis(fig[1, 3]; xlabel = "compressor speed (× design)", ylabel = "COP1 change [%]",
               title = "COP1 vs no refrigerant Δp")
    ref = Dict(x.speed => x.COP1 for x in first(sw).rows)
    for (c, ls) in zip(sw, (:dash, :solid, :dot))
        s = [x.speed for x in c.rows]
        lines!(ax1, s, [x.Δp/1e3 for x in c.rows]; color = :firebrick, linestyle = ls, label = c.label)
        lines!(ax2, s, [x.Tr_in - x.Te for x in c.rows]; color = :firebrick, linestyle = ls)
        lines!(ax3, s, [100 * (x.COP1 / ref[x.speed] - 1) for x in c.rows]; color = :firebrick, linestyle = ls)
    end
    Legend(fig[2, 1:3], ax1; orientation = :horizontal, framevisible = false)
    return fig
end
