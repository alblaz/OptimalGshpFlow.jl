"""
EvaporatorParameters

`refrigerant::AbstractRefrigerant`         refrigerant properties (spline `Refrigerant` or `ClapeyronRefrigerant`),
`secondaryfluid::SecondaryFluid`           secondary fluid properties,
`plateheatexchanger::PlateHeatExchanger`   plate heat exchanger geometry,
`hei`                                      refrigerant enthalpy at the inlet of the evaporator,
`Vs`                                       volume flow rate of secondary fluid (m^3/s), 
`Vref`                                     volume flow rate of the refrigerant (m^3/s), 
`ΔTsup`                                    refrigerant super heating degree, 
`Ts_h`                                     secondary fluid temperature (in Kelvin) at the inlet of the evaporator
"""
@with_kw struct EvaporatorParameters{R <: AbstractRefrigerant, S <: SecondaryFluid, P <: PlateHeatExchanger, T <: Real}
    refrigerant::R
    secondaryfluid::S
    plateheatexchanger::P
    hei::T
    Vs::T
    Vref::T
    ΔTsup::T
    Ts_h::T
end
# scalars of mixed types (e.g. Int, Float64, ForwardDiff Duals) are promoted to a common type
EvaporatorParameters(r::AbstractRefrigerant, sf::SecondaryFluid, phx::PlateHeatExchanger, x::Vararg{Real, 5}) =
    EvaporatorParameters(r, sf, phx, promote(x...)...)


"""
BoreholeParameters

`secondaryfluid::SecondaryFluid`           secondary fluid properties,
`Vs`                                       volume flow rate of secondary fluid (m^3/s), 
`Tin`                                      borehole inlet temperature
`borehole::Borehole`                       borehole heat exchanger configuration,
`λground`                                  ground conductivity,
`Tb`                                       secondary fluid temperature (in Kelvin) at the inlet of the evaporator
"""
@with_kw struct BoreholeParameters{S <: SecondaryFluid, B <: Borehole, T <: Real}
    secondaryfluid::S
    Vs::T
    Tin::T
    borehole::B
    λground::T
    Tb::T
end
function BoreholeParameters(sf::SecondaryFluid, Vs::Real, Tin::Real, borehole::Borehole, λground::Real, Tb::Real)
    Vs, Tin, λground, Tb = promote(Vs, Tin, λground, Tb)
    return BoreholeParameters(sf, Vs, Tin, borehole, λground, Tb)
end


"""
EvaporatorBoreholeParameters    

`refrigerant::AbstractRefrigerant`         refrigerant properties (spline `Refrigerant` or `ClapeyronRefrigerant`),
`secondaryfluid::SecondaryFluid`           secondary fluid properties,
`plateheatexchanger::PlateHeatExchanger`   plate heat exchanger geometry,
`hei`                                      refrigerant enthalpy at the inlet of the evaporator,
`Vs`                                       volume flow rate of secondary fluid (m^3/s), 
`Vref`                                     volume flow rate of the refrigerant (m^3/s), 
`ΔTsup`                                    refrigerant super heating degree, 
`borehole::Borehole`                       borehole heat exchanger configuration,
`λground`                                  ground conductivity,
`Tb`                                       secondary fluid temperature (in Kelvin) at the inlet of the evaporator
"""
@with_kw struct EvaporatorBoreholeParameters{R <: AbstractRefrigerant, S <: SecondaryFluid, P <: PlateHeatExchanger,
                                              B <: Borehole, T <: Real}
    refrigerant::R
    secondaryfluid::S
    plateheatexchanger::P
    hei::T
    Vs::T
    Vref::T
    ΔTsup::T
    borehole::B
    λground::T
    Tb::T
end
function EvaporatorBoreholeParameters(r::AbstractRefrigerant, sf::SecondaryFluid, phx::PlateHeatExchanger, hei::Real,
                                      Vs::Real, Vref::Real, ΔTsup::Real, borehole::Borehole, λground::Real, Tb::Real)
    hei, Vs, Vref, ΔTsup, λground, Tb = promote(hei, Vs, Vref, ΔTsup, λground, Tb)
    return EvaporatorBoreholeParameters(r, sf, phx, hei, Vs, Vref, ΔTsup, borehole, λground, Tb)
end


"""
    `evaporator_model!(du, u, p::EvaporatorParameters)`

Heat balances and equalities defining the evaporator's physical behavior considering a two region steady state model for the heat exchanger with moving boundary bewteen the two-phase region and the superheated region.

`du`   vector that contains the result of equalities resulting from heat balances and equations, 
`u`    vector of variables, 
`p::EvaporatorParameters`  parameters    

variables:
    `u0 = pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m`

`pme`     refrigerant evaporation pressure, mapped with `pressure_transform(refrigerant)` (with `NestedPressure()`; see `evaporator_transforms`), 
`Qe`      evaporator power exchange in the two phase region, 
`Qsup`    evaporator power exchange in the superheated region, 
`ϵme`     heat exchanger effectiveness in the two-phase region, mapped with `Unit()`, 
`ϵmsup`   heat exchanger effectiveness in the superheat region, mapped with `Unit()`, 
`xms`     coordinate of the moving boundary, mapped with `Unit()`, 
`Tsf_l`   secondary fluid temperature (in Kelvin) at the evaporator outlet, 
`Tsf_m`   secondary fluid temperature (in Kelvin) at the moving boundary, 
"""
function evaporator_model!(du, u, p::EvaporatorParameters, formulation::PropertyFormulation = NestedPressure())
    @unpack refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, Ts_h = p
    rs = evaporator_refrigerant_state!(du, formulation, refrigerant, u, 8, ΔTsup)
    evaporator_residuals!(du, u, Ts_h, rs, refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup)
    return nothing
end

"""
    evaporator_residuals!(du, u, Ts_h, rs, refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup)

Residuals `du[1:8]` of the two-zone (two-phase + superheat) evaporator for the unknowns
`u[2:8] = Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m`, with secondary-fluid inlet temperature `Ts_h`
and the refrigerant state `rs` (see [`evaporator_refrigerant_state!`](@ref), which also handles
`u[1]`). Shared by [`evaporator_model!`](@ref) and [`evaporatorborehole_model!`](@ref).
"""
function evaporator_residuals!(du, u, Ts_h, rs, refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup)
    # mapped unknowns
    x = (Qe = u[2], Qsup = u[3], ϵe = from_R(Unit(), u[4]), ϵsup = from_R(Unit(), u[5]),
         xs = from_R(Unit(), u[6]), Ts_l = u[7], Ts_m = u[8])
    return evaporator_balances!(du, x, Ts_h, rs, refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup)
end

"""
    evaporator_balances!(du, x, Ts_h, rs, refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup)

The eight evaporator residuals (`du[1:8]`, unscaled: W for the heat balances, dimensionless for
the effectiveness equations) on physical values `x = (; Qe, Qsup, ϵe, ϵsup, xs, Ts_l, Ts_m)`.
"""
function evaporator_balances!(du, x, Ts_h, rs, refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup)
    A = plateheatexchanger.A
    (; Qe, Qsup, ϵe, ϵsup, xs, Ts_l, Ts_m) = x

    # refrigerant states
    (; pe, Te, hes, hsup) = rs
    ρref_sup = rs.ρsup

    # two-phase: locally constant capacities and U value
    Tsm_e = (Ts_l + Ts_m)/2
    Ue, hs_e, href_e = evaporator_uvalue_twophase(Tsm_e, Vs, Qe, pe, xs, secondaryfluid, refrigerant, plateheatexchanger)
    Cs_e = Cmin_e = secondaryfluid.ρ(Tsm_e)* secondaryfluid.cp(Tsm_e) * Vs   # capacity of secondary fluid

    # superheat
    Tsm_sup = (Ts_h + Ts_m)/2
    Usup, hs_sup, href_sup = evaporator_uvalue_superheat(Tsm_sup, Vs, Vref, rs, secondaryfluid, plateheatexchanger)
    Cs_sup = secondaryfluid.ρ(Tsm_sup) * secondaryfluid.cp(Tsm_sup) * Vs
    Cref_sup = ρref_sup*Vref*(hsup - hes)/ΔTsup   # mean capacity of the refrigerant in the superheating region
    Cmin_sup = min(Cref_sup, Cs_sup)
    Cmax_sup = max(Cref_sup, Cs_sup)
    Cr = Cmin_sup/Cmax_sup

    # heat balance equations (two-phase)
    du[1] = Cs_e * (Ts_m - Ts_l) - Qe             # heat balance secondary fluid
    du[2] = ρref_sup * Vref * (hes - hei) - Qe    # heat balance refrigerant
    du[3] = ϵe * Cmin_e * (Ts_m - Te) - Qe        # ϵ-NTU method
    du[4] = 1 - exp(-Ue * A * xs / Cmin_e) - ϵe

    # heat balance equations (superheater)
    du[5] = Cs_sup * (Ts_h - Ts_m) - Qsup         # heat balance secondary fluid
    du[6] = Cref_sup * ΔTsup - Qsup
    du[7] = ϵsup * Cmin_sup * (Ts_h - Te) - Qsup  # ϵ-NTU
    du[8] = (1 - exp((-Usup*A*(1 - xs)/Cmin_sup*(1-Cr))))/(1 - Cr *(exp(-Usup*A*(1 - xs)/Cmin_sup*(1-Cr)))) - ϵsup

    return nothing
end


"""
    mixture_boiling_factor(refrigerant, p, q, h_ideal; β_L = 3e-4, B = 1.0)

Ratio of the nucleate boiling heat-transfer coefficient of a zeotropic mixture to the ideal one
`h_ideal` (the pure-fluid correlation evaluated with the mixture properties), at pressure `p` and
heat flux `q` [W/m²] (Thome & Shakir 1987, AIChE Symp. Ser. 83(257), 46):

`F = [1 + (h_ideal/q) ΔT_bp (1 − exp(−B q/(ρ_l Δh_lv β_L)))]⁻¹`

with the boiling range `ΔT_bp = T_dew − T_bubble`, the saturated-liquid density `ρ_l`, the latent
heat `Δh_lv = h_dew − h_bubble` (overall composition), the liquid mass-transfer coefficient
`β_L = 0.0003 m/s` and `B = 1` (all the heat goes into evaporation). The light component
evaporates first and the liquid at the bubble interface is enriched in the heavy one, whose
higher local boiling point lowers the driving temperature difference. 1 for a pure fluid.
"""
mixture_boiling_factor(::AbstractRefrigerant, p, q, h; kw...) = one(h)
mixture_boiling_factor(r::ClapeyronRefrigerant, p, q, h; kw...) = ispure(r) ? one(h) : _thome_shakir(r, p, q, h; kw...)

function _thome_shakir(r::ClapeyronRefrigerant, p, q, h; β_L = 3e-4, B = 1.0)
    Tb, vl = saturation_bubble(r, p)
    Td, vv = saturation_dew(r, p)
    Δh_lv = (Clapeyron.VT_enthalpy(r.eos, vv, Td, r.z) - Clapeyron.VT_enthalpy(r.eos, vl, Tb, r.z)) / r.M
    c = B / (r.M / vl * Δh_lv * β_L)                  # (1 − e^{−c q})/q = c E(c q), finite at q = 0
    E, _ = _glide_EG(c * q)
    return inv(1 + h * (Td - Tb) * c * E)
end

"""
`evaporator_uvalue_twophase(Tsm_e, Vs, Qe, pe, xs, secondaryfluid, refrigerant, plateheatexchanger)`

evaluate U value `Ue` and the heat transfer coefficients `hs_e`, `href_e` in the two-phase region of the heat exchanger

`Tsm_e` mean secondary fluid temperature in Kelvin,
`pe` evaporation pressure,
`Vs` secondary fluid flow rate (m^3/s), 
`Qe` evaporator power exchange in the two phase region, 
`xs` coordinate of the moving boundary,
`secondaryfluid::SecondaryFluid`,
`refrigerant::AbstractRefrigerant`, 
`plateheatexchanger::PlateHeatExchanger`

return values

`Ue`, `hs_e`, `href_e`
"""
function evaporator_uvalue_twophase(Tsm_e, Vs, Qe, pe, xs, secondaryfluid, refrigerant, plateheatexchanger)
    A = plateheatexchanger.A

    # heat_transfer_coefficient secondary fluid side
    secondaryfluid_sm_e =  secondaryfluid.ρ(Tsm_e), secondaryfluid.μ(Tsm_e), secondaryfluid.λ(Tsm_e), secondaryfluid.cp(Tsm_e)
    hs_e = h_platehex_wanniarachchi(Vs, plateheatexchanger, secondaryfluid_sm_e)    

    # heat_transfer_coefficient refrigerant side
    pr = pe/critical_pressure(refrigerant)
    q = abs(Qe)/(A * xs)
    href_e = 1.5 * 55*pr^(0.12)/((-log10(pr))^0.55 * (molar_mass(refrigerant)*1000)^0.5) * q^0.67
    href_e *= mixture_boiling_factor(refrigerant, pe, q, href_e)     # blends: mass-transfer resistance

    Ue =  1/(1/hs_e + 1/href_e) 

    return Ue, hs_e, href_e
end


"""
`evaporator_uvalue_superheat(Tsm_sup, pe, Vs, Vref, ΔTsup, secondaryfluid, refrigerant, plateheatexchanger)

evaluate the U value and the heat transfer coefficients `Usup`, `hs_sup`, `href_sup` in the superheat region of the heat exchanger

`Tsm_e` mean secondary fluid temperature in Kelvin,
`pe` evaporation pressure,
`Vs` secondary fluid flow rate (m^3/s), 
`Vref` refrigerant flow rate (m^3/s), 
`ΔTsup` superheat,
`Qsup` evaporator power exchange in the su[erheat region, 
`secondaryfluid::SecondaryFluid`,
`refrigerant::AbstractRefrigerant`, 
`plateheatexchanger::PlateHeatExchanger`

return values

`Usup`, `hs_sup`, `href_sup`
"""

function evaporator_uvalue_superheat(Tsm_sup, pe, Vs, Vref, ΔTsup, secondaryfluid, refrigerant, plateheatexchanger)
    rs = (ρsup = ρ_superheated(refrigerant, pe, ΔTsup), cpm = cp_superheated(refrigerant, pe, ΔTsup/2),
          μm = μ_superheated(refrigerant, pe, ΔTsup/2), λm = λ_superheated(refrigerant, pe, ΔTsup/2))
    return evaporator_uvalue_superheat(Tsm_sup, Vs, Vref, rs, secondaryfluid, plateheatexchanger)
end

# refrigerant properties from the state `rs`: density at the outlet, cp, μ, λ at mid superheat
function evaporator_uvalue_superheat(Tsm_sup, Vs, Vref, rs, secondaryfluid, plateheatexchanger)

    # heat_transfer_coefficient secondary fluid side
    ρs_sup = secondaryfluid.ρ(Tsm_sup)
    cps_sup = secondaryfluid.cp(Tsm_sup)
    μs_sup = secondaryfluid.μ(Tsm_sup)
    λs_sup = secondaryfluid.λ(Tsm_sup)
    properties_s_sup = ρs_sup, μs_sup, λs_sup, cps_sup

    hs_sup = h_platehex_wanniarachchi(Vs , plateheatexchanger, properties_s_sup)

    # heat_transfer_coefficient refrigerant side
    propertiesR_sup = rs.ρsup, rs.μm, rs.λm, rs.cpm

    href_sup = h_platehex_wanniarachchi(Vref, plateheatexchanger, propertiesR_sup)
    Usup =  1/(1/hs_sup + 1/href_sup)

    return Usup, hs_sup,href_sup
end


"""
     `evaporator_model(u,p)`

unpack, remap and complete the solution obtained by solving evaporator_model!(du,u,p) such that all the parameters of the evaporator are available.
The return value is of type EvaporatorSolution 

"""
function evaporator_model(u,p::EvaporatorParameters)
    pme, Qe, Qsup, ϵme, ϵmsup, xms, Ts_l, Ts_m = u
    @unpack refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, Ts_h = p
    A = plateheatexchanger.A 

    # transform mapped variables
    pe = from_R(pressure_transform(refrigerant), pme)
    ϵe = from_R(Unit(), ϵme)
    ϵsup = from_R(Unit(), ϵmsup)
    xs  = from_R(Unit(), xms)
    
    # two-phase
    Te = T_dew(refrigerant, pe)
    Tsm_e = (Ts_l + Ts_m)/2
    Cs_e = Cmin_e = secondaryfluid.ρ(Tsm_e )* secondaryfluid.cp(Tsm_e) * Vs   #  capacity of secondary fluid 
    Cref_e = Inf
    
    Ue, hs_e, href_e  =  evaporator_uvalue_twophase(Tsm_e, Vs, Qe, pe, xs, secondaryfluid, refrigerant, plateheatexchanger)
    NTU_e = Ue*A*xs/Cmin_e
    twophase_heatexchanger = HeatExchangerSolution(Qe, Ts_m, Ts_l, Te, Te, Cs_e, Cref_e, Vs, Vref, hs_e, href_e, Ue*A*xs, ϵe, NTU_e)

    # superheating region
    
    # superheat
    Tsm_sup = (Ts_h + Ts_m)/2   
    Usup, hs_sup, href_sup = evaporator_uvalue_superheat(Tsm_sup, pe,  Vs, Vref, ΔTsup, secondaryfluid, refrigerant, plateheatexchanger)
    ρs_sup = secondaryfluid.ρ(Tsm_sup)      
    cps_sup = secondaryfluid.cp(Tsm_sup) 
    Cs_sup = ρs_sup * cps_sup * Vs
    hsup = h_superheated(refrigerant, pe, ΔTsup)
    ρref_sup = ρ_superheated(refrigerant, pe, ΔTsup)
    hes = h_dew(refrigerant, pe)
    Cref_sup = ρref_sup*Vref*(hsup - hes)/ΔTsup   #Cref_sup (mean capacity of the refrigerant in the superheating region)
    Cmin_sup = min(Cref_sup,Cs_sup)

    NTU_sup = Usup*A*(1-xs)/Cmin_sup
    
    Tch_sup = Te + ΔTsup
    superheat_heatexchanger = HeatExchangerSolution(Qsup, Ts_h, Ts_m, Tch_sup, Te, Cs_sup, Cref_sup, Vs, Vref, hs_sup, href_sup, Usup*A*(1-xs), ϵsup, NTU_sup)

    return EvaporatorSolution(plateheatexchanger,xs, twophase_heatexchanger, superheat_heatexchanger)
end



"""
    `evaporatorborehole_model!(du, u, p)`

Heat balances and equalities defining the evaporator's physical behavior considering a two region steady state model for the heat exchanger with moving boundary bewteen the two-phase region and the superheated region.

`du`   vector that contains the result of equalities resulting from heat balances and equations, 
`u`    vector of variables, 
`p`    vector of parameters. 

variables:
    `u0 = pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m, Tsf_h`

`pme`     refrigerant evaporation pressure, mapped with `pressure_transform(refrigerant)` (with `NestedPressure()`; see `evaporator_transforms`), 
`Qe`      evaporator power exchange in the two phase region, 
`Qsup`    evaporator power exchange in the superheated region, 
`ϵme`     heat exchanger effectiveness in the two-phase region, mapped with `Unit()`, 
`ϵmsup`   heat exchanger effectiveness in the superheat region, mapped with `Unit()`, 
`xms`     coordinate of the moving boundary, mapped with `Unit()`, 
`Tsf_l`   secondary fluid temperature (in Kelvin) at the evaporator outlet, 
`Tsf_m`   secondary fluid temperature (in Kelvin) at the moving boundary, 
`Ts_h`    secondary fluid temperature (in Kelvin) at the inlet of the evaporator
"""
function evaporatorborehole_model!(du, u, p::EvaporatorBoreholeParameters, formulation::PropertyFormulation = NestedPressure())
    @unpack refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, borehole, λground, Tb = p
    ## u = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m, Tsf_h] (+ formulation unknowns)
    Ts_l, Ts_m, Ts_h = u[7], u[8], u[9]

    rs = evaporator_refrigerant_state!(du, formulation, refrigerant, u, 9, ΔTsup)
    evaporator_residuals!(du, u, Ts_h, rs, refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup)

    # borehole model (Hellström), coefficient at the mean secondary-fluid temperature
    Tsm_e = (Ts_l + Ts_m)/2
    kin = uniform_Tb_coeffs(borehole, λground, secondaryfluid.ρ(Tsm_e)*Vs, (Ts_l + Ts_h)/2, secondaryfluid)
    du[9] = kin*Ts_l - Ts_h + (1 - kin)*Tb

    return nothing
end


"""
     `evaporatorborehole_model(u,p)`

unpack, remap and complete the solution obtained by solving evaporatorborehole_model!(du,u,p) such that all the parameters of the evaporator and of the borehole are available.
The return value is of type EvaporatorBoreholeSolution 

"""
function evaporatorborehole_model(u,p::EvaporatorBoreholeParameters)
    pme, Qe, Qsup, ϵme, ϵmsup, xms, Ts_l, Ts_m, Ts_h = u
    @unpack refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, borehole, λground, Tb = p

    _u = u[1:8]
    _p = EvaporatorParameters(refrigerant=refrigerant, secondaryfluid=secondaryfluid, plateheatexchanger=plateheatexchanger, 
                              hei=hei, Vs=Vs, Vref=Vref, ΔTsup=ΔTsup, Ts_h=u[9])
    

    evaporator_solution = evaporator_model(_u, _p)

    # BoreholeSolution(borehole, λground,Tin,Tout,Tb,Vs,properties)
    
    borehole_solution = BoreholeSolution(borehole, λground,Ts_l, Ts_h, Tb, Vs, secondaryfluid)
    return EvaporatorBoreholeSolution(evaporator_solution, borehole_solution)
end



"""
     `borehole_model!(du,u,p)`

Heat balance at the borehole according to Hellström.

`du`   vector that contains the result of equalities resulting from heat balances and equations, 
`u`    vector of variables, 
`p`    vector of parameters.

variables:
    `u0 = [Tout]`

`Tout`    borehole outlet temperature 
"""
function borehole_model!(du,u,p::BoreholeParameters)
    @unpack secondaryfluid, Vs, Tin, borehole, λground, Tb = p
    kin = uniform_Tb_coeffs(borehole, λground, secondaryfluid.ρ((u[1] + Tin)/2)*Vs, (u[1] + Tin)/2, secondaryfluid)
    du[1] = kin*Tin - 1*u[1] + (1-kin)* Tb
end


"""
     `borehole_model(u,p::BoreholeParameters)`

unpack, remap and complete the solution obtained by solving borehole_model!(du,u,p) such that all the parameters of the borehole are available.
The return value is of type BoreholeSolution 

"""
function borehole_model(u,p::BoreholeParameters)
    @unpack secondaryfluid, Vs, Tin, borehole, λground, Tb = p
    Tout = u[1]
    return BoreholeSolution(borehole, λground, Tin, Tout, Tb, Vs, secondaryfluid)
end


# function borehole_model!(secondaryfluid, Vs, Tin, borehole, λg, Tb)

#     kin = uniform_Tb_coeffs(borehole, λground, secondaryfluid.ρ(Tsm_e )*Vs, (u[7] + u[9]/2), secondaryfluid)
#     Tout =    kin*u[2]  + (1-kin)* Tb 

# end