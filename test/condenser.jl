# Condenser component (step 9c): desuperheating + condensation against water, counter-flow
using NonlinearSolve

@testset "condenser" begin
    r = SPL_R290
    water = SecondaryFluid(range(273.15 + 5, 273.15 + 70, length = 40), "Water")
    # vapour temperature from (p, h): inverse of h_pT on the vapour side
    pc0 = PropsSI("P", "T", 313.15, "Q", 0, "R290")
    Tv = 330.0
    @test vapour_T_ph(r, pc0, h_pT(r, pc0, Tv)) ≈ Tv atol = 1e-8

    # design point of the examples: R290 discharge from Te = 0 °C (ΔTsup 3.5 K, ηis 0.7) at 40 °C,
    # floor heating 30 → 35 °C
    m, h_in, Q1 = 0.0221, 6.532e5, 7641.0
    d = condenser_design(Q1, 313.15, 303.15, 308.15, water)
    ΔTlm = 5 / log(10 / 5)
    @test d.A ≈ Q1 / (2000 * ΔTlm) && d.Vs ≈ Q1 / (water.ρ(305.65) * water.cp(305.65) * 5)
    phx = PlateHeatExchanger(d.A, 0.5, 0.15; b = 0.002, ϕ = 1.189, β = 30.)
    cond = Condenser(; refrigerant = r, secondaryfluid = water, plateheatexchanger = phx, Vs = d.Vs, Qscale = Q1)
    @test keys(unknowns(cond)) == (:pc, :Qd, :Qc, :xd, :Tw_m) && keys(ports(cond)) == (:Tw_in, :Tw_out, :m, :h_in)
    S = ComponentSystem((condenser = cond,); nodes = (Tw_in = (:condenser => :Tw_in,), Tw_out = (:condenser => :Tw_out,),
            m = (:condenser => :m,), h_in = (:condenser => :h_in,)), fixed = (Tw_in = 303.15, m = m, h_in = h_in))
    u0 = initial_guess(S; condenser = (Tc = 313.15, Qd = 0.15Q1, Qc = 0.85Q1, Tw_m = 307.), nodes = (Tw_out = 308.15,))
    sol = solve(NonlinearProblem(S, u0), NewtonRaphson(); abstol = 1e-9)
    @test SciMLBase.successful_retcode(sol)
    x = physical(S, sol.u).condenser
    pv = port_values(S, sol.u, :condenser)
    z = condenser_zones(cond, x, pv)
    # energy: refrigerant from the discharge to saturated liquid, water 30 → ≈ 35 °C
    @test x.Qd + x.Qc ≈ m * (h_in - h_bubble(r, x.pc)) rtol = 1e-9
    @test x.Qd ≈ m * (h_in - h_dew(r, x.pc)) rtol = 1e-9
    @test 308.0 < pv.Tw_out < 308.3
    @test 313.15 < z.Td < 316.15                                     # sized for 40 °C at U = 2000
    @test z.Tdis > z.Td > x.Tw_m > pv.Tw_in && 0 < x.xd < 0.5
    # pure fluid: condensation zone with ε = 1 − e^{−NTU}: Td − Tw_m = (Td − Tw_in) e^{−NTU}
    @test z.Td - x.Tw_m ≈ (z.Td - pv.Tw_in) * exp(-z.condensation.NTU) rtol = 1e-8
    @test isempty(validity_issues(S, sol.u))
    # pinch at the refrigerant outlet when the water comes in at the bubble temperature
    @test any(i -> occursin("refrigerant outlet", i.message), validity(cond, x, merge(pv, (Tw_in = z.Tb,))))
end
