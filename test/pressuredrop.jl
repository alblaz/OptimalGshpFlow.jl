# Pipe pressure drop: no jump at the laminar-turbulent transition (it showed up as a 40 % drop of
# the pump power in the optimal-flow study when the borehole Re crossed 2300).
@testset "pressure drop: laminar-turbulent continuity" begin
    δ = 1e-6
    for Rec in (2300., 4000.)
        @test OptimalGshpFlow.fDarcy(Rec - δ) ≈ OptimalGshpFlow.fDarcy(Rec + δ) rtol = 1e-6
    end
    # unchanged regimes: laminar Hagen-Poiseuille and the turbulent smooth-pipe fit
    @test OptimalGshpFlow.fDarcy(1500.) == 64 / 1500
    @test OptimalGshpFlow.fDarcy(2e4) ≈ (1 / (0.838 * OptimalGshpFlow.lambertw(0.629 * 2e4)))^2

    # pressure drop continuous in the flow rate across Re = 2300
    sf = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), "INCOMP::MEA-28%")
    H, D, T = 400., 0.04 * (1 - 1/17), 268.15
    V(Re) = Re * sf.μ(T) / sf.ρ(T) * π * D / 4
    @test pressure_drop(H, V(2300. - δ), D, sf, T) ≈ pressure_drop(H, V(2300. + δ), D, sf, T) rtol = 1e-6
    @test pumping_power(0.3, H, V(2300. - δ), D, sf, T) ≈ pumping_power(0.3, H, V(2300. + δ), D, sf, T) rtol = 1e-6

    # ForwardDiff Duals pass through (number-type generic)
    @test OptimalGshpFlow.ForwardDiff.derivative(OptimalGshpFlow.fDarcy, 3000.) ≈
          (OptimalGshpFlow.fDarcy(4000.) - 64 / 2300) / 1700
end

# The original sharp switch is still available, and the choice reaches pressure_drop / pumping_power.
@testset "pressure drop: transition models" begin
    fD = OptimalGshpFlow.fDarcy
    fturb(Re) = (1 / (0.838 * OptimalGshpFlow.lambertw(0.629 * Re)))^2
    for Re in (1500., 2299., 2300., 3000., 3999., 4000., 2e4)
        @test fD(Re, SharpTransition()) == (Re < 2300 ? 64 / Re : fturb(Re))
        @test fD(Re) == fD(Re, GnielinskiTransition()) == fD(Re, SmoothTransition())
        Re ≥ 4000 && @test fD(Re) == fturb(Re)
    end
    @test fD(2300. + 1e-6, SharpTransition()) > 1.6 * fD(2300. - 1e-6, SharpTransition())

    sf = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), "INCOMP::MEA-28%")
    H, D, T = 400., 0.04 * (1 - 1/17), 268.15
    V3000 = 3000. * sf.μ(T) / sf.ρ(T) * π * D / 4                 # Re = 3000
    @test pressure_drop(H, V3000, D, sf, T; transition = SharpTransition()) > pressure_drop(H, V3000, D, sf, T)
    @test pumping_power(0.3, H, V3000, D, sf, T; transition = SharpTransition()) ≈
          pressure_drop(H, V3000, D, sf, T; transition = SharpTransition()) * V3000 / 0.3
end

# Pressure drop and pumping power of solved components (secondary-fluid loop)
@testset "pressure drop and pumping power of solutions" begin
    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), "INCOMP::MEA-28%")
    phx = PlateHeatExchanger(1.5, 0.5, 0.15)
    bh  = SingleUPipeBorehole(H = 200., D = 4., rpi = 0.0188)
    V   = 4e-4
    hx(Q, Thh, Thl) = HeatExchangerSolution(Q, Thh, Thl, 270., 268., 1e3, 1e3, V, 1e-4, 1e3, 1e3, 300., 0.6, 1.0)
    es  = OptimalGshpFlow.EvaporatorSolution(phx, 0.9, hx(5e3, 275., 272.), hx(300., 276., 275.))
    bs  = BoreholeSolution(bh, 3., 272., 276., 281.15, V, 1.5e3, 6e3, 0.1)
    Δp_e = pressure_drop_evap(V, phx, sf, (276. + 272.)/2)
    Δp_b = pressure_drop(400., V, 0.0376, sf, 274.; transition = bh.transition)
    @test pressure_drop(es, sf) == Δp_e
    @test pressure_drop(bs, sf) ≈ Δp_b
    @test pumping_power(es, sf, 0.3) ≈ Δp_e * V / 0.3
    @test pumping_power(bs, sf, 0.3) ≈ Δp_b * V / 0.3
    @test (@inferred pumping_power(es, sf, 0.3)) isa Float64
    @test (@inferred pumping_power(bs, sf, 0.3)) isa Float64
    # the borehole's transition model reaches the friction factor
    bs_sharp = BoreholeSolution(SingleUPipeBorehole(H = 200., D = 4., rpi = 0.0188, transition = SharpTransition()),
                                3., 272., 276., 281.15, V, 1.5e3, 6e3, 0.1)
    @test pressure_drop(bs_sharp, sf) ≈ pressure_drop(400., V, 0.0376, sf, 274.; transition = SharpTransition())
end

# Plate evaporator, brine side: the channel cross-section is width × gap (W b per channel), as in
# the heat transfer correlation (V/(H W)); it used the plate length L instead of W.
@testset "plate pressure drop: channel velocity" begin
    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), "INCOMP::MEA-28%")
    phx = PlateHeatExchanger(1.5, 0.5, 0.15)                # L = 0.5 m, W = 0.15 m
    V, T = 4e-4, 274.
    Nc = (phx.np - 1)/2
    dh = 2phx.b/phx.ϕ
    v = V/(phx.W * phx.b * Nc)
    @test v ≈ V/(phx.H * phx.W)                             # same velocity as h_platehex_wanniarachchi
    Re = sf.ρ(T) * v * dh / sf.μ(T)
    Δp_ref = OptimalGshpFlow.f_Martin(Re, phx.β) * phx.L/dh * v^2 * sf.ρ(T)/2
    @test pressure_drop_evap(V, phx, sf, T) ≈ Δp_ref
end
