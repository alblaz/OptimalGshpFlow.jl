
Trange = range(-15. + 273.15, 30 + 273.15, length = 40)
ethanol =  "INCOMP::MEA-28%"
secondaryfluid = OptimalGshpFlow.secondary_fluid(Trange, ethanol)
Ttest = [270., 273, 288., 290]

@testset "secondaryfluid" begin
    for T in Ttest 
        @test secondaryfluid.cp.(T) .≈ PropsSI("C","T",T,"P",101325,ethanol) 
        @test secondaryfluid.λ.(T)  .≈ PropsSI("CONDUCTIVITY","T",T,"P",101325,ethanol) 
        @test secondaryfluid.μ.(T)  .≈ PropsSI("V","T",T,"P",101325,ethanol)  atol = 1e-8
        @test secondaryfluid.ρ.(T)  .≈ PropsSI("D","T",T,"P",101325,ethanol) 
    end
    
    #test the extrapolation
    @test secondaryfluid.cp.(maximum(Trange)) .≈  secondaryfluid.cp.(maximum(Trange) + 5.)
    @test secondaryfluid.λ.(maximum(Trange))  .≈  secondaryfluid.λ.(maximum(Trange) + 5.)
    @test secondaryfluid.μ.(maximum(Trange))  .≈  secondaryfluid.μ.(maximum(Trange) + 5.)
    @test secondaryfluid.ρ.(maximum(Trange))  .≈  secondaryfluid.ρ.(maximum(Trange) + 5.)

    @test secondaryfluid.cp.(minimum(Trange)) .≈  secondaryfluid.cp.(minimum(Trange) + -10.)
    @test secondaryfluid.λ.(minimum(Trange))  .≈  secondaryfluid.λ.(minimum(Trange) + -10.)
    @test secondaryfluid.μ.(minimum(Trange))  .≈  secondaryfluid.μ.(minimum(Trange) + -10.)
    @test secondaryfluid.ρ.(minimum(Trange))  .≈  secondaryfluid.ρ.(minimum(Trange) + -10.)
    
end



refrigerant_tp = OptimalGshpFlow.refrigerant_twophaseregion("R134a", 1. );
p = [1e5, 3e5, 1e6];   # inside the sampled range (pr_min = 1e-2 ⇒ p ≥ 4.06e4 Pa for R134a)

@testset "refrigerant_twophase" begin
    for p in p 
        @test refrigerant_tp.h(log10(p))  .≈  PropsSI("H", "P",  p, "Q", 1, "R134a")  rtol = 1e-6
        @test refrigerant_tp.T(log10(p))  .≈  PropsSI("T", "P",  p, "Q", 1, "R134a")  rtol = 1e-6
        @test refrigerant_tp.μ(log10(p))  .≈  PropsSI("V", "P",  p, "Q", 1, "R134a")  rtol = 1e-6
        @test refrigerant_tp.λ(log10(p))  .≈  PropsSI("CONDUCTIVITY", "P",  p, "Q", 1, "R134a")  rtol = 1e-6
        @test refrigerant_tp.cp(log10(p)) .≈  PropsSI("C", "P",  p, "Q", 1, "R134a")  rtol = 1e-6
        @test refrigerant_tp.ρ(log10(p))  .≈  PropsSI("D", "P",  p, "Q", 1, "R134a")  rtol = 1e-6
    end
end


refrigerant_sh = OptimalGshpFlow.refrigerant_superheated("R134a") 
p, ΔTsup = 3e5, 3.
Tsat = PropsSI("T", "P",  p, "Q", 1, "R134a")

@testset "refrigerant_superheated" begin
    for p in p 
        @test refrigerant_sh.h(log10(p), ΔTsup)  .≈  PropsSI("H", "P",  p, "T", Tsat + ΔTsup, "R134a")  rtol = 1e-6
        @test refrigerant_sh.μ(log10(p), ΔTsup)  .≈  PropsSI("V", "P",  p, "T", Tsat + ΔTsup, "R134a")  rtol = 1e-6
        @test refrigerant_sh.λ(log10(p), ΔTsup)  .≈  PropsSI("CONDUCTIVITY", "P",  p, "T", Tsat + ΔTsup, "R134a")  rtol = 1e-6
        @test refrigerant_sh.cp(log10(p), ΔTsup) .≈  PropsSI("C", "P",  p, "T", Tsat + ΔTsup, "R134a")  rtol = 1e-6
        @test refrigerant_sh.ρ(log10(p), ΔTsup)  .≈  PropsSI("D", "P",  p, "T", Tsat + ΔTsup, "R134a")  rtol = 1e-6
    end
end

# Refrigerant construction must work for common heat-pump refrigerants and the
# spline must reproduce CoolProp over the operating range (evaporator and condenser).
@testset "Refrigerant construction" begin
    for f in ("R290", "R134a", "R1234yf", "R600a", "R744", "R410A")
        @test (Refrigerant(f); true)
    end
    # CoolProp cannot evaluate the saturated-vapour conductivity of R32 below ≈ -30 °C
    # (conformal state solver fails), so the sampling range has to start higher.
    @test (Refrigerant("R32"; pr_min = 0.05); true)

    # keyword arguments are forwarded to the interpolants
    r = Refrigerant("R134a"; pr_min = 0.02)
    @test r.pr_min == 0.02
end

@testset "Refrigerant spline accuracy" begin
    f = "R134a"
    r = Refrigerant(f)
    pc = r.pcrit
    for pr in (0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7)       # Tsat ≈ -33 °C … +80 °C
        p  = pr * pc
        lp = log10(p)
        Ts = PropsSI("T", "P", p, "Q", 1, f)
        @test r.saturated_vapor.h(lp)  ≈ PropsSI("H", "P", p, "Q", 1, f) rtol = 1e-6
        @test r.saturated_vapor.cp(lp) ≈ PropsSI("C", "P", p, "Q", 1, f) rtol = 1e-2
        for ΔT in (2.5, 5., 10.)
            @test r.superheated_vapor.h(lp, ΔT)  ≈ PropsSI("H", "P", p, "T", Ts + ΔT, f) rtol = 1e-5
            @test r.superheated_vapor.cp(lp, ΔT) ≈ PropsSI("C", "P", p, "T", Ts + ΔT, f) rtol = 1e-2
            @test r.superheated_vapor.λ(lp, ΔT)  ≈ PropsSI("CONDUCTIVITY", "P", p, "T", Ts + ΔT, f) rtol = 1e-2
        end
    end
end

# B4: the evaporation-pressure unknown is mapped onto `pressure_bounds(r)`. For the spline
# backend this must be the sampled pressure range: the original code used the field `pr_min`
# (a reduced pressure, 1e-2) as a bound in Pa, so the unknown could reach pressures far below
# the splines, where they extrapolate flat and the Jacobian is singular.
@testset "pressure bounds of the spline Refrigerant (B4)" begin
    for (name, pr_min) in (("R134a", 1e-2), ("R290", 2e-2))
        r = Refrigerant(name; pr_min = pr_min)
        _, ps = OptimalGshpFlow.log10p_grid(name; pr_min = pr_min)
        pmin, pmax = pressure_bounds(r)
        @test pmin ≈ first(ps) rtol = 1e-12
        @test pmax ≈ last(ps) rtol = 1e-12
        @test pmin ≥ pr_min * r.pcrit
    end
end
