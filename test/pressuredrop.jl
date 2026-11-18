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
