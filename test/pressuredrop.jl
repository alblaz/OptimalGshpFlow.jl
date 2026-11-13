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
end
