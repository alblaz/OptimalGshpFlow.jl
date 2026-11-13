# The borehole equation inside evaporatorborehole_model! (du[9]) must evaluate the
# borehole coefficient at the mean secondary-fluid temperature (Tin + Tout)/2.
@testset "evaporator-borehole coupling" begin
    sfname = "INCOMP::MEA-28%"
    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), sfname)
    ref = Refrigerant("R134a")
    ev  = EvaporatorDesign(6000., 273., 273 + 40, "R134a", sfname; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15)
    bh  = SingleUPipeBorehole(H = 200., D = 4., λg = 3., pipe_position = ((0.04, 0.), (-0.04, 0.)))
    λground, Tb = 3.5, 283.15
    p = EvaporatorBoreholeParameters(refrigerant = ref, secondaryfluid = sf, plateheatexchanger = phx,
            hei = ev.hei, Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., borehole = bh, λground = λground, Tb = Tb)

    pe = PropsSI("P", "T", 273.15, "Q", 1., "R134a")
    pm = to_R(pressure_transform(ref), pe)
    Tin, Tm, Tout = 274.15, 276.15, 278.15                  # Ts_l, Ts_m, Ts_h
    u  = [pm, 5000., 500., 0., 0., 0., Tin, Tm, Tout]
    du = zeros(9)
    evaporatorborehole_model!(du, u, p)

    mf  = sf.ρ((Tin + Tm) / 2) * ev.Vs                      # mass flow as in the coupled model
    kin = OptimalGshpFlow.uniform_Tb_coeffs(bh, λground, mf, (Tin + Tout) / 2, sf)
    @test du[9] ≈ kin * Tin - Tout + (1 - kin) * Tb atol = 1e-10
end
