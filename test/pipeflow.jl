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
        n = mode == "heating" ? 0.4 : 0.3          # Dittus-Boelter: 0.4 for a heated fluid
        @test Nu(2e4, Pr; mode = mode) ≈ 0.023 * 2e4^0.8 * Pr^n
    end
end

# Both transition models are available: the smooth blend (default) and the sharp
# switch laminar -> Gnielinski (2300) -> Dittus-Boelter (1e4) used originally.
@testset "laminar-turbulent transition models" begin
    Nu = OptimalGshpFlow.nusselt_circularpipe
    function nu_sharp_reference(Re, Pr, mode)
        Re < 2300 && return 3.657
        if Re < 1e4
            f = (1.58 * log(Re) - 3.28)^(-2)
            return (0.5 * f * (Re - 1000.) * Pr) / (1 + 12.7 * sqrt(0.5 * f) * (Pr^(2/3) - 1))
        end
        return 0.023 * Re^0.8 * Pr^(mode == "heating" ? 0.4 : 0.3)
    end

    for Pr in (0.7, 7., 30.), mode in ("heating", "cooling"), Re in (1500., 2299., 2300., 2500., 5000., 9999., 1e4, 2e4)
        @test Nu(Re, Pr; mode = mode, transition = SharpTransition()) == nu_sharp_reference(Re, Pr, mode)
        @test Nu(Re, Pr; mode = mode) == Nu(Re, Pr; mode = mode, transition = SmoothTransition())
    end
    # the sharp model really has a jump at Re = 2300
    @test Nu(2300. + 1e-6, 7.; transition = SharpTransition()) > 2 * Nu(2300. - 1e-6, 7.; transition = SharpTransition())

    # the choice reaches the borehole model
    sf = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), "INCOMP::MEA-28%")
    bh_smooth = SingleUPipeBorehole(H = 200., D = 4.)
    bh_sharp  = SingleUPipeBorehole(H = 200., D = 4., transition = SharpTransition())
    @test bh_smooth.transition isa SmoothTransition
    @test bh_sharp.transition isa SharpTransition
    rp, T = bh_smooth.rpi, 278.15
    Re(mf) = 2 * mf / (π * rp * sf.μ(T))
    mf_trans = 3000. * π * rp * sf.μ(T) / 2          # Re = 3000: transition range
    mf_lam   = 1500. * π * rp * sf.μ(T) / 2          # Re = 1500: laminar
    @test Re(mf_trans) ≈ 3000.
    @test h_circularpipe(mf_trans, T, rp, sf; transition = SharpTransition()) > h_circularpipe(mf_trans, T, rp, sf)
    k_smooth(mf) = OptimalGshpFlow.uniform_Tb_coeffs(bh_smooth, 3., mf, T, sf)
    k_sharp(mf)  = OptimalGshpFlow.uniform_Tb_coeffs(bh_sharp, 3., mf, T, sf)
    @test k_sharp(mf_trans) != k_smooth(mf_trans)
    @test k_sharp(mf_lam) == k_smooth(mf_lam)
end

# Gnielinski (VDI 2013) turbulent branch joined to the laminar value by a C¹ blend
@testset "GnielinskiTransition" begin
    Nu = OptimalGshpFlow.nusselt_circularpipe
    G = GnielinskiTransition()
    gn(Re, Pr) = (ξ = (1.8 * log10(Re) - 1.5)^(-2); (ξ/8) * Re * Pr / (1 + 12.7 * sqrt(ξ/8) * (Pr^(2/3) - 1)))
    for Pr in (0.7, 7., 30.)
        f = Re -> Nu(Re, Pr; transition = G)
        @test f(1500.) == 3.657
        @test f(2e4) ≈ gn(2e4, Pr)
        @test f(1e4) ≈ gn(1e4, Pr)
        # analytic slope of the turbulent correlation
        @test OptimalGshpFlow.dnusselt_gnielinski2013(1e4, Pr) ≈ dderivative(Re -> gn(Re, Pr), 1e4)
        # C¹: value and slope continuous at both ends of the transition range
        δ = 1e-6
        for Rec in (2300., 1e4)
            @test f(Rec - δ) ≈ f(Rec + δ) rtol = 1e-6
            @test dderivative(f, Rec - δ) ≈ dderivative(f, Rec + δ) rtol = 1e-5 atol = 1e-9
        end
        Res = range(500., 3e4, 500)
        @test all(diff(f.(Res)) .>= 0)
    end
    @test SingleUPipeBorehole(H = 200., D = 4., transition = G).transition isa GnielinskiTransition
end
