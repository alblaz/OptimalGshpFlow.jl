# SingleUPipeBorehole must be constructible with the default keyword values and
# give a borehole inlet coefficient k_i in (0, 1).
@testset "SingleUPipeBorehole" begin
    bh = SingleUPipeBorehole(H = 200., D = 4., λg = 3., pipe_position = ((0.04, 0.0), (-0.04, 0.0)))
    @test bh isa SingleUPipeBorehole
    @test OptimalGshpFlow.get_H(bh) == 200.

    Trange = range(-15. + 273.15, 30 + 273.15, length = 40)
    sf = SecondaryFluid(Trange, "INCOMP::MEA-28%")
    k_i = OptimalGshpFlow.uniform_Tb_coeffs(bh, 3.5, 0.5, 278.15, sf)
    @test 0 < k_i < 1
end
