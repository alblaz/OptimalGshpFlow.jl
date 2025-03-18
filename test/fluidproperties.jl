
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
p = [1e3,1e4,1e5];

@testset "refrigerant_twophase" begin
    for p in p 
        @test refrigerant_tp.h(log10(p))  .≈  PropsSI("H", "P",  p, "Q", 1, "R134a")
        @test refrigerant_tp.T(log10(p))  .≈  PropsSI("T", "P",  p, "Q", 1, "R134a")
        @test refrigerant_tp.μ(log10(p))  .≈  PropsSI("V", "P",  p, "Q", 1, "R134a")
        @test refrigerant_tp.λ(log10(p))  .≈  PropsSI("CONDUCTIVITY", "P",  p, "Q", 1, "R134a")
        @test refrigerant_tp.cp(log10(p)) .≈  PropsSI("C", "P",  p, "Q", 1, "R134a")
        @test refrigerant_tp.ρ(log10(p))  .≈  PropsSI("D", "P",  p, "Q", 1, "R134a")  atol = 1e-6
    end
end


refrigerant_sh = OptimalGshpFlow.refrigerant_superheated("R134a") 
p, ΔTsup = 10000., 3.
Tsat = PropsSI("T", "P",  p, "Q", 1, "R134a")

@testset "refrigerant_superheated" begin
    for p in p 
        @test refrigerant_sh.h(log10(p), ΔTsup)  .≈  PropsSI("H", "P",  p, "T", Tsat + ΔTsup, "R134a")
        @test refrigerant_sh.μ(log10(p), ΔTsup)  .≈  PropsSI("V", "P",  p, "T", Tsat + ΔTsup, "R134a")
        @test refrigerant_sh.λ(log10(p), ΔTsup)  .≈  PropsSI("CONDUCTIVITY", "P",  p, "T", Tsat + ΔTsup, "R134a")
        @test refrigerant_sh.cp(log10(p), ΔTsup) .≈  PropsSI("C", "P",  p, "T", Tsat + ΔTsup, "R134a")
        @test refrigerant_sh.ρ(log10(p), ΔTsup)  .≈  PropsSI("D", "P",  p, "T", Tsat + ΔTsup, "R134a")  atol = 1e-7
    end
end