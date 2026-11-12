
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

@testset "sigmoid inverse" begin
    # the inverse of the sigmoid composed with the sigmoid must be the identity
    xs = range(0.01, 0.99, 50)       # points in (0,1)
    ys = range(-10., 10., 50)        # points in (-∞,+∞)

    @test sigmoid_01_to_R(0.5) == 0.
    @test sigmoid_R_to_01(0.) == 0.5

    for x in xs
        @test sigmoid_R_to_01(sigmoid_01_to_R(x)) ≈ x rtol = 1e-12
    end
    for y in ys
        @test sigmoid_01_to_R(sigmoid_R_to_01(y)) ≈ y rtol = 1e-10 atol = 1e-12
    end

    # bounded interval (a,b), e.g. a pressure range in Pa
    a, b = 1e3, 4.25e6
    for z in range(a + 1., b - 1., 50)
        @test sigmoid_R_to_ab(sigmoid_ab_to_R(z, a, b), a, b) ≈ z rtol = 1e-10
    end
    for y in ys
        @test sigmoid_ab_to_R(sigmoid_R_to_ab(y, a, b), a, b) ≈ y rtol = 1e-8 atol = 1e-8
    end

    # same round trips for the tanh maps (regression guard)
    for x in xs
        @test tanhmap_R_to_01(tanhmap_01_to_R(x)) ≈ x rtol = 1e-12
    end
    for z in range(a + 1., b - 1., 50)
        @test tanhmap_R_to_ab(tanhmap_ab_to_R(z, a, b), a, b) ≈ z rtol = 1e-10
    end
end