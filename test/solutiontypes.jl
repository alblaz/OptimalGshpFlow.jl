# Solution types must be concretely typed: with untyped fields, post-processing that reads them
# (e.g. Clapeyron saturation/flash calls in a compressor model) is inferred with `Any` arguments,
# which can take practically forever to compile.
using ForwardDiff

@testset "concrete solution types" begin
    allconcrete(x) = isconcretetype(typeof(x)) && all(isconcretetype, fieldtypes(typeof(x)))
    sfname = "INCOMP::MEA-28%"
    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), sfname)
    ev  = EvaporatorDesign(6000., 273., 273 + 40, "R290", sfname; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15)
    r   = SPL_R290
    p   = EvaporatorParameters(refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx, hei = ev.hei,
                               Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., Ts_h = 278.)
    u   = [to_R(pressure_transform(r), PropsSI("P", "T", 270., "Q", 1., "R290")), 5000., 300., 0., 0., 0., 274., 276.]

    es = evaporator_model(u, p)
    @test allconcrete(es) && allconcrete(es.twophase_region) && allconcrete(es.superheat_region)
    @test (@inferred (s -> s.twophase_region.Tcl)(es)) isa Float64
    @test es.twophase_region.Cc == Inf                       # Cref = Inf in the two-phase zone
    @test es.ΔT_l == es.twophase_region.Thl - es.twophase_region.Tcl

    # ForwardDiff Duals pass through (mixed with Float64 fields such as Vs, Inf)
    esd = evaporator_model(ForwardDiff.Dual.(u, 1.), p)
    @test allconcrete(esd) && esd.twophase_region.Tcl isa ForwardDiff.Dual
    @test ForwardDiff.value(esd.twophase_region.Q) ≈ es.twophase_region.Q

    bh = SingleUPipeBorehole(H = 200., D = 4.)
    bs = BoreholeSolution(bh, 3., 270., 272., 268.15, ev.Vs, sf)
    @test allconcrete(bs) && bs.Qb ≈ 2bs.Cs
    @test (@inferred (s -> s.Tout)(bs)) isa Float64
    @test allconcrete(EvaporatorBoreholeSolution(es, bs))
end

# Parameter structs and SecondaryFluid must be concretely typed as well: otherwise the solution
# returned by `evaporator_model` / `evaporator_solution` is inferred abstractly inside a function
# (e.g. a parameter study), and post-processing it with ClapeyronRefrigerant does not finish
# compiling.
@testset "concrete parameter types" begin
    allconcrete(x) = isconcretetype(typeof(x)) && all(isconcretetype, fieldtypes(typeof(x)))
    sfname = "INCOMP::MEA-28%"
    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), sfname)
    ev  = EvaporatorDesign(6000., 273., 273 + 40, "R290", sfname; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15)
    bh  = SingleUPipeBorehole(H = 200., D = 4.)
    r   = ClapeyronRefrigerant("R290")
    @test allconcrete(sf)
    @test (@inferred (s -> s.ρ(270.))(sf)) isa Float64
    # spline refrigerant: untyped fields made every property call `Any`; a study calling it inside
    # a function then took ≈ 450 s to compile (backedges for abstract calls, Julia 1.13)
    @test allconcrete(SPL_R290) && allconcrete(SPL_R290.saturated_vapor) && allconcrete(SPL_R290.superheated_vapor)
    @test (@inferred T_dew(SPL_R290, 4e5)) isa Float64
    @test (@inferred h_pT(SPL_R290, 4e5, 270.; phase = :vapor)) isa Float64

    p = EvaporatorParameters(refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx, hei = ev.hei,
                             Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., Ts_h = 278.)
    @test allconcrete(p)
    @test p.Ts_h === 278.0                                    # scalars promoted to a common type
    @test allconcrete(EvaporatorParameters(p, Vs = 1))        # Parameters.jl reconstruction
    @test allconcrete(BoreholeParameters(secondaryfluid = sf, Vs = ev.Vs, Tin = 270, borehole = bh, λground = 3., Tb = 280.))
    pb = EvaporatorBoreholeParameters(refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx, hei = ev.hei,
                                      Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., borehole = bh, λground = 3., Tb = 280.)
    @test allconcrete(pb)

    # inferred return types of the post-processing used by parameter studies
    @test isconcretetype(only(Base.return_types(evaporator_model, (Vector{Float64}, typeof(p)))))
    e = Evaporator(refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx, hei = ev.hei, Vs = ev.Vs,
                   Vref = ev.Vref, ΔTsup = 5.)
    b = BoreholeHX(borehole = bh, secondaryfluid = sf, Vs = ev.Vs, λground = 3.)
    @test allconcrete(e) && allconcrete(b)
    x  = (Te = 268., Qe = 5e3, Qsup = 300., ϵe = 0.5, ϵsup = 0.5, xs = 0.5, Ts_m = 276.)
    pv = (Ts_in = 278., Ts_out = 274.)
    @test isconcretetype(only(Base.return_types(evaporator_solution, (typeof(e), typeof(x), typeof(pv)))))
    @test isconcretetype(only(Base.return_types(borehole_solution, (typeof(b), typeof((Tin = 270., Tout = 272., Tb = 280.))))))

    # physical values of a system solution (input of the post-processing above)
    sys = ComponentSystem((evaporator = e, borehole = b);
        nodes = (T_cold = (:evaporator => :Ts_out, :borehole => :Tin), T_warm = (:evaporator => :Ts_in, :borehole => :Tout),
                 T_ground = (:borehole => :Tb,)), fixed = (T_ground = 280.,))
    u = initial_guess(sys; evaporator = (Te = 268., Qe = 5e3, Qsup = 300., Ts_m = 276.), nodes = (T_cold = 274., T_warm = 279.))
    ph = @inferred physical(sys, u)
    @test keys(ph) == (:evaporator, :nodes)                  # no entry for the borehole (no unknowns)
    @test ph.evaporator.Te ≈ 268. && ph.nodes.T_warm ≈ 279.
end

# The refrigerant itself must be inferred: code that builds it inside a compiled block (a
# function, or a top-level block with loops such as a @testset) otherwise propagates a partially
# abstract ClapeyronRefrigerant into Clapeyron/EntropyScaling, and inference did not finish
# (test/operating_envelope.jl run on its own).
@testset "inferred ClapeyronRefrigerant construction" begin
    @test isconcretetype(only(Base.return_types(ClapeyronRefrigerant, (String,))))
    r = @inferred ClapeyronRefrigerant("R290")
    @test r.transport isa CoolPropTransport
    # EntropyScaling transport: the backend is inferred, the RefpropRES parameter type is not
    # (EntropyScaling.RefpropRES itself is not inferable); known limitation (mixtures)
    @test only(Base.return_types(n -> ClapeyronRefrigerant(n; transport = EntropyScalingTransport()), (String,))) <:
          ClapeyronRefrigerant{<:Any, EntropyScalingTransport}
    # vector of components, default transport chosen at run time: a small Union at worst
    T = only(Base.return_types(c -> ClapeyronRefrigerant(c; mass_fractions = [0.5, 0.5]), (Vector{String},)))
    @test T isa Union || isconcretetype(T)
end

# Purity is part of the refrigerant type, so that the default formulation (SaturationTemperature
# for a pure fluid, NestedPressure for a mixture) and hence the Evaporator, the system and its
# unknowns are inferred; with a run-time choice they were a Union of two systems with different
# unknowns (Te vs pe), which made a parameter study / test block practically impossible to compile.
@testset "purity in the refrigerant type" begin
    r = ClapeyronRefrigerant("R290")
    mx = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785])
    @test OptimalGshpFlow.ispure(r) && !OptimalGshpFlow.ispure(mx)
    @test ClapeyronRefrigerant(["R290"]) isa typeof(r)                 # a 1-element vector is pure too
    @test only(Base.return_types(default_formulation, (typeof(r),))) === SaturationTemperature
    @test only(Base.return_types(default_formulation, (typeof(mx),))) === NestedPressure
    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), "INCOMP::MEA-28%")
    ev  = EvaporatorDesign(6000., 273., 273 + 40, "R290", "INCOMP::MEA-28%"; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15)
    mk(r) = Evaporator(refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx, hei = ev.hei, Vs = ev.Vs,
                       Vref = ev.Vref, ΔTsup = 5.)
    @test isconcretetype(only(Base.return_types(mk, (typeof(r),))))
    @test isconcretetype(only(Base.return_types(n -> mk(ClapeyronRefrigerant(n)), (String,))))
end
