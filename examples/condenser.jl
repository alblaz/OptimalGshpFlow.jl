#= Plate condenser against floor-heating water (step 9c): R290 and R454C.

The condenser of the optimal-flow study (`optimalflow_components.jl`): sized for the R290 design
duty with `condenser_design` (U = 2000, condensation at 40 °C, water 30 → 35 °C), same plate type as
the evaporator. The refrigerant inlet is the compressor discharge at the design evaporator state
(dew temperature `DESIGN.Te`, superheat `ΔTSUP`, `refrigerant_design` displacement, ηis). The
condensing pressure is the solution (not fixed at 40 °C as in the study): for a blend it shows where
bubble and dew temperatures actually fall against the water, i.e. which condenser convention of
step 5h (`:bubble`, `:midpoint`, `:dew`) the water side supports. Not coupled to the evaporator yet
(step 9d).

Usage:
    include("examples/condenser.jl")
    r454c = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785], name = "R454C")
    report_condenser([ClapeyronRefrigerant("R290"), r454c])
    sw = water_sweep(r454c); report_water_sweep(sw)
=#
include(joinpath(@__DIR__, "optimalflow_components.jl"))

const WATER = SecondaryFluid(range(273.15 + 5, 273.15 + 70, length = 40), "Water")
const TW_IN, TW_OUT = 273.15 + 30, 273.15 + 35

# compressor discharge of the design evaporator state, for the refrigerant `r` (its design)
function design_discharge(r; design = refrigerant_design(r))
    pe = OptimalGshpFlow._p_dew(r, DESIGN.Te)
    pc = OptimalGshpFlow._p_bubble(r, design.Tc)
    c = OptimalGshpFlow.compression(r, pe, DESIGN.Te + ΔTSUP, 0.0, pc, ηIS)
    m = ρ_pT(r, pe, DESIGN.Te + ΔTSUP; phase = :vapor) * design.Vref
    return (m = m, h_in = c.h_out, Q1 = m * (c.h_out - h_bubble(r, pc)))
end

# the condenser of the R290 design (area and water flow), for any refrigerant
const CONDENSER_DESIGN = let d = design_discharge(ClapeyronRefrigerant(REFRIGERANT))
    merge(condenser_design(d.Q1, TCOND, TW_IN, TW_OUT, WATER), (Q1 = d.Q1,))
end
const CPHX = PlateHeatExchanger(CONDENSER_DESIGN.A, 0.5, 0.15; b = 0.002, ϕ = 1.189, β = 30.)

function condenser_system(r, m, h_in; Tw_in = TW_IN, Vs = CONDENSER_DESIGN.Vs)
    cond = Condenser(; refrigerant = r, secondaryfluid = WATER, plateheatexchanger = CPHX, Vs = Vs,
                     Qscale = CONDENSER_DESIGN.Q1)
    return ComponentSystem((condenser = cond,); nodes = (Tw_in = (:condenser => :Tw_in,), Tw_out = (:condenser => :Tw_out,),
            m = (:condenser => :m,), h_in = (:condenser => :h_in,)), fixed = (Tw_in = Tw_in, m = m, h_in = h_in))
end

condenser_guess(S, Q) = initial_guess(S; condenser = (Tc = TCOND, Qd = 0.15Q, Qc = 0.85Q, Tw_m = TW_OUT - 0.5),
                                      nodes = (Tw_out = TW_OUT,))

# one condenser solution: temperatures and zone quantities
function condenser_point(r, m, h_in; Tw_in = TW_IN, u0 = nothing)
    S = condenser_system(r, m, h_in; Tw_in = Tw_in)
    sol = solve(NonlinearProblem(S, u0 === nothing ? condenser_guess(S, m * (h_in - 2.6e5)) : u0), TrustRegion(); abstol = 1e-9)
    SciMLBase.successful_retcode(sol) || return nothing
    x = physical(S, sol.u).condenser
    pv = port_values(S, sol.u, :condenser)
    z = condenser_zones(S.components.condenser, x, pv)
    return (u = sol.u, pc = x.pc, Td = z.Td, Tb = z.Tb, Tdis = z.Tdis, Tw_in = pv.Tw_in, Tw_out = pv.Tw_out, Tw_m = x.Tw_m,
            xd = x.xd, Qd = x.Qd, Qc = x.Qc, U_c = z.condensation.U, h_c = z.condensation.h_ref, φ = z.condensation.φ,
            G = OptimalGshpFlow.channel_mass_flux(CPHX, m), issues = String[i.message for i in validity_issues(S, sol.u)])
end

"""
    report_condenser(refrigerants)

Each refrigerant at its design discharge in the R290-sized condenser: condensing pressure, dew,
bubble and mean temperature, the bubble temperatures of the three conventions of
[`condenser_bubble_T`](@ref) (5h), discharge temperature, zone split and coefficients.
"""
function report_condenser(refrigerants)
    @printf("condenser: A = %.3f m², water %.3f l/s, 30 → 35 °C at the R290 design duty %.0f W\n",
            CONDENSER_DESIGN.A, 1e3CONDENSER_DESIGN.Vs, CONDENSER_DESIGN.Q1)
    for r in refrigerants
        d = design_discharge(r)
        c = condenser_point(r, d.m, d.h_in)
        c === nothing && (println(r.name, ": no solution"); continue)
        @printf("%-6s pc %.3f bar: T_dew %.2f, T_bubble %.2f, mean %.2f °C (glide %.2f K); Tdis %.1f °C; water out %.2f °C\n",
                r.name, c.pc/1e5, c.Td - 273.15, c.Tb - 273.15, (c.Td + c.Tb)/2 - 273.15, c.Td - c.Tb, c.Tdis - 273.15, c.Tw_out - 273.15)
        @printf("       Q %.0f W (desuperheat %.0f, xd %.3f); condensation h_ref %.0f, U %.0f W/(m² K), φ %.2f, G %.1f kg/(m² s)  %s\n",
                c.Qd + c.Qc, c.Qd, c.xd, c.h_c, c.U_c, c.φ, c.G, join(c.issues, "; "))
        if !OptimalGshpFlow.ispure(r)
            @printf("       convention bubble temperatures at 40 °C: bubble %.2f, midpoint %.2f, dew %.2f °C\n",
                    condenser_bubble_T(r, :bubble) - 273.15, condenser_bubble_T(r, :midpoint) - 273.15, condenser_bubble_T(r, :dew) - 273.15)
        end
    end
end

"""
    water_sweep(r; Tw_ins = 273.15 .+ (25:1:40)) -> Vector

The design discharge of `r` condensed against water entering at each `Tw_in` (same water flow),
continuation from the coldest.
"""
function water_sweep(r; Tw_ins = 273.15 .+ (25:1:40))
    d = design_discharge(r)
    rows = typeof(condenser_point(r, d.m, d.h_in; Tw_in = first(Tw_ins)))[]
    u = nothing
    for Tw in Tw_ins
        c = condenser_point(r, d.m, d.h_in; Tw_in = Tw, u0 = u)
        c === nothing && (@warn "no solution at Tw_in = $(Tw - 273.15) °C"; continue)
        push!(rows, c); u = c.u
    end
    return rows
end

function report_water_sweep(rows)
    @printf("%7s %8s %8s %8s %8s %8s %7s\n", "Tw_in", "Tw_out", "T_bub", "T_dew", "Tdis", "pc[bar]", "xd")
    for c in rows
        @printf("%7.1f %8.2f %8.2f %8.2f %8.1f %8.3f %7.3f\n", c.Tw_in - 273.15, c.Tw_out - 273.15, c.Tb - 273.15,
                c.Td - 273.15, c.Tdis - 273.15, c.pc/1e5, c.xd)
    end
end
