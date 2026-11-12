@testset "nusselt_circularpipe" begin
    Nu = OptimalGshpFlow.nusselt_circularpipe
    δ = 1e-6

    for Pr in (0.7, 7., 30.), mode in ("heating", "cooling")
        # no jump at the laminar/transition and transition/turbulent boundaries
        for Rec in (2300., 1e4)
            @test Nu(Rec - δ, Pr; mode = mode) ≈ Nu(Rec + δ, Pr; mode = mode) rtol = 1e-6
        end

        # Nu does not decrease with Re
        Res = range(500., 3e4, 500)
        Nus = [Nu(Re, Pr; mode = mode) for Re in Res]
        @test all(diff(Nus) .>= 0)

        # unchanged regimes: laminar (fully developed) and turbulent (Dittus-Boelter)
        @test Nu(1500., Pr; mode = mode) == 3.657
        n = mode == "heating" ? 0.3 : 0.4
        @test Nu(2e4, Pr; mode = mode) ≈ 0.023 * 2e4^0.8 * Pr^n
    end
end
