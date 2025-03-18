
@testset "mapping.jl" begin

    #test linear mapping
    x,a,b = 0.1, -1., 3.
    @test linearmap_01_to_ab(x,a,b) >= a 
    @test linearmap_01_to_ab(x,a,b) <= b 

    y,a,b = 2.3, 1., 3.
    @test linearmap_ab_to_01(y,a,b) <= 1.
    @test linearmap_ab_to_01(y,a,b) >= 0.
       
    x,a,b = 0.1, -1., 3.
    @test linearmap_ab_to_01(linearmap_01_to_ab(x,a,b),a,b) == x
    @test linearmap_01_to_ab(linearmap_ab_to_01(x,a,b),a,b) ≈ x

    #test tanhmap_R_to_01 and tanhmap_R_to_01
    @test tanhmap_R_to_01(0.) == 0.5
    @test tanhmap_R_to_01(-3.) >= 0.
    @test tanhmap_R_to_01(5.) <= 1.

    @test tanhmap_01_to_R(0.5) == 0.
    @test tanhmap_01_to_R(0.1) <= 0
    @test tanhmap_01_to_R(0.6) >= 0.

    # test tanhmap_R_to_ab and tanhmap_ab_to_R
    @test tanhmap_R_to_ab(tanhmap_ab_to_R(x,a,b), a,b) ≈ x rtol = 1e-12
    @test tanhmap_ab_to_R(tanhmap_R_to_ab(x,a,b), a,b) ≈ x rtol = 1e-12
end