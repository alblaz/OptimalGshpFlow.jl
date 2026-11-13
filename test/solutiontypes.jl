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
    r   = Refrigerant("R290")
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
