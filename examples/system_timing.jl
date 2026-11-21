#= Cost of the three systems of the study, same refrigerant and operating point:

- `:evaporator`: standalone evaporator, brine inlet fixed (7 unknowns);
- `:evaporator_borehole`: evaporator + borehole, `heatpump` of the optimal-flow study;
- `:cycle`: closed cycle (evaporator + compressor + condenser + borehole, `heatpump_cycle`).

Each is solved with NonlinearSolve (`TrustRegion`, as the study, from the example's guess): first
solve (in this order, so it includes the compilation not shared with the previous systems),
residual, warm solve (same guess), iterations; the closed cycle also with `NewtonRaphson` (a plain
Newton loop in `examples/cycle.jl` was the workaround for B16, now `cycle_solve`). Brine inlet of the
standalone evaporator ≈ the warm brine of the evaporator + borehole solution (`TB − 3.5`). Run each refrigerant in a fresh session for comparable
first-call times.

Usage:
    include("examples/system_timing.jl")
    rows = system_timing(ClapeyronRefrigerant("R290")); report_timing(rows)
=#
include(joinpath(@__DIR__, "cycle.jl"))

# standalone evaporator at the brine inlet temperature `Ts_in`
function evaporator_system(r; Vf = 1.0, Ts_in = TB - 3.5, design = refrigerant_design(r), sf = BRINE_PROPS)
    e = Evaporator(refrigerant = r, secondaryfluid = sf, plateheatexchanger = PHX, hei = design.hei,
                   Vs = Vf * DESIGN.Vs, Vref = design.Vref, ΔTsup = ΔTSUP)
    return ComponentSystem((evaporator = e,); nodes = (T_out = (:evaporator => :Ts_out,), T_in = (:evaporator => :Ts_in,)),
                           fixed = (T_in = Ts_in,))
end
evaporator_guess(S, Ts_in) = initial_guess(S; evaporator = (Te = Ts_in - 8, Qe = 0.9DESIGN.Qnom, Qsup = 0.1DESIGN.Qnom,
                                                            Ts_m = Ts_in - 0.5), nodes = (T_out = Ts_in - 3,))

# first solve, residual, warm solve; results through barriers (`solve` is inferred as `Any`)
function time_system(S, u0; alg = TrustRegion(), abstol = 1e-9)
    t_first = @elapsed sol = solve(NonlinearProblem(S, u0), alg; abstol = abstol)
    u = convert(Vector{Float64}, sol.u)::Vector{Float64}
    du = similar(u)
    OptimalGshpFlow.residual!(du, u, S)
    t_res = minimum(@elapsed(OptimalGshpFlow.residual!(du, u, S)) for _ in 1:100)
    t_solve = minimum(@elapsed(solve(NonlinearProblem(S, u0), alg; abstol = abstol)) for _ in 1:3)
    return (n = length(S), ok = SciMLBase.successful_retcode(sol)::Bool, iterations = Int(sol.stats.nsteps)::Int,
            t_first = t_first, t_res = t_res, t_solve = t_solve, u = u)
end

function system_timing(r; Tb = TB, log = devnull)
    _log(log, r.name, ": evaporator")
    Se = evaporator_system(r; Ts_in = Tb - 3.5)
    e = time_system(Se, evaporator_guess(Se, Tb - 3.5))
    _log(log, r.name, ": evaporator + borehole")
    Sb = heatpump(r, 1.0; Tb = Tb)
    b = time_system(Sb, first_guess(Sb, Tb))
    _log(log, r.name, ": closed cycle")
    Sc = cycle_system(r; Tb = Tb)
    c = time_system(Sc, cycle_guess(Sc, Tb))
    _log(log, r.name, ": closed cycle, NewtonRaphson")
    cn = time_system(Sc, cycle_guess(Sc, Tb); alg = NewtonRaphson())
    _log(log, r.name, ": done")
    return [(name = r.name, system = :evaporator, e...), (name = r.name, system = :evaporator_borehole, b...),
            (name = r.name, system = :cycle, c...), (name = r.name, system = :cycle_newton, cn...)]
end

function report_timing(rows)
    @printf("%-6s %-20s %3s %3s %4s %10s %11s %11s\n", "", "system", "n", "ok", "its", "first (s)", "residual µs", "solve ms")
    for x in rows
        @printf("%-6s %-20s %3d %3s %4d %10.1f %11.1f %11.2f\n", x.name, x.system, x.n, x.ok ? "yes" : "no", x.iterations,
                x.t_first, 1e6x.t_res, 1e3x.t_solve)
    end
end
