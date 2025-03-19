"""
EvaporatorParameters

`refrigerant::Refrigerant`                 refrigerant properties,
`secondaryfluid::SecondaryFluid`           secondary fluid properties,
`plateheatexchanger::PlateHeatExchanger`   plate heat exchanger geometry,
`hei`                                      refrigerant enthalpy at the inlet of the evaporator,
`Vs`                                       volume flow rate of secondary fluid (m^3/s), 
`Vref`                                     volume flow rate of the refrigerant (m^3/s), 
`ΔTsup`                                    refrigerant super heating degree, 
`Ts_h`                                     secondary fluid temperature (in Kelvin) at the inlet of the evaporator
"""
@with_kw struct EvaporatorParameters
    refrigerant::Refrigerant
    secondaryfluid::SecondaryFluid
    plateheatexchanger::PlateHeatExchanger
    hei
    Vs
    Vref
    ΔTsup 
    Ts_h
end


"""
BoreholeParameters

`secondaryfluid::SecondaryFluid`           secondary fluid properties,
`Vs`                                       volume flow rate of secondary fluid (m^3/s), 
`Tin`                                      borehole inlet temperature
`borehole::Borehole`                       borehole heat exchanger configuration,
`λground`                                  ground conductivity,
`Tb`                                       secondary fluid temperature (in Kelvin) at the inlet of the evaporator
"""
@with_kw struct BoreholeParameters 
    secondaryfluid::SecondaryFluid
    Vs
    Tin 
    borehole::Borehole
    λground
    Tb
end


"""
EvaporatorBoreholeParameters    

`refrigerant::Refrigerant`                 refrigerant properties,
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
@with_kw struct EvaporatorBoreholeParameters 
    refrigerant::Refrigerant
    secondaryfluid::SecondaryFluid
    plateheatexchanger::PlateHeatExchanger
    hei
    Vs
    Vref
    ΔTsup 
    borehole::Borehole
    λground
    Tb
end


"""
    `evaporator_model!(du, u, p::EvaporatorParameters)`

Heat balances and equalities defining the evaporator's physical behavior considering a two region steady state model for the heat exchanger with moving boundary bewteen the two-phase region and the superheated region.

`du`   vector that contains the result of equalities resulting from heat balances and equations, 
`u`    vector of variables, 
`p::EvaporatorParameters`  parameters    

variables:
    `u0 = pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m`

`pme`     refrigerant evaporation pressure mapped using tanh transformation, 
`Qe`      evaporator power exchange in the two phase region, 
`Qsup`    evaporator power exchange in the superheated region, 
`ϵme`     heat exchanger efficiency in the two phase region mapped using tanh transformation, 
`ϵmsup`   heat exchanger efficiency in the two phase region mapped using tanh transformation, 
`xms`     coordinate of the moving boundary mapped using tanh transformation, 
`Tsf_l`   secondary fluid temperature (in Kelvin) at the evaporator outlet, 
`Tsf_m`   secondary fluid temperature (in Kelvin) at the moving boundary, 
"""
function evaporator_model!(du, u, p::EvaporatorParameters)

    @unpack refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, Ts_h = p
    A = plateheatexchanger.A
    ## u is a vector that contains the unknonws of the problem    
    ## u = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m]

    # first we compute model variables that are considered locally "constant": 
    # These are Cs_e, Ue, Cs_sup, Cref_sup, Cr,  U_sup,

    # two-phase
    pe = tanhmap_R_to_ab(u[1], refrigerant.pr_min, refrigerant.pcrit)
    xs = tanhmap_R_to_01(u[6])
    Qe = u[2]
    Tsm_e = (u[7] + u[8])/2

    Ue, hs_e, href_e  =  evaporator_uvalue_twophase(Tsm_e, Vs, Qe, pe, xs, secondaryfluid, refrigerant, plateheatexchanger)
    Cs_e = Cmin_e = secondaryfluid.ρ(Tsm_e )* secondaryfluid.cp(Tsm_e) * Vs   #  capacity of secondary fluid 

    # superheat
    Tsm_sup = (Ts_h + u[8])/2   
    Usup, hs_sup, href_sup = evaporator_uvalue_superheat(Tsm_sup, pe,  Vs, Vref, ΔTsup, secondaryfluid, refrigerant, plateheatexchanger)
    ρs_sup = secondaryfluid.ρ(Tsm_sup)      
    cps_sup = secondaryfluid.cp(Tsm_sup) 
    Cs_sup = ρs_sup * cps_sup * Vs
    hsup =  refrigerant.superheated_vapor.h(log10(pe),ΔTsup)     # fixed super heating of 5 degrees (find H for that super heating degree)
    ρref_sup  = refrigerant.superheated_vapor.ρ(log10(pe), ΔTsup)
    hes = refrigerant.saturated_vapor.h(log10(pe))
    Cref_sup = ρref_sup*Vref*(hsup - hes)/ΔTsup   #Cref_sup (mean capacity of the refrigerant in the superheating region)
    Cmin_sup = min(Cref_sup,Cs_sup)
    Cmax_sup = max(Cref_sup,Cs_sup)
    Cr = Cmin_sup/Cmax_sup    

    # Heat balance equations (two-phase)
    du[1] = Cs_e * (u[8] - u[7]) - u[2]  # heat balance secondary fluid
    du[2] = refrigerant.superheated_vapor.ρ(log10(tanhmap_R_to_ab(u[1], refrigerant.pr_min, refrigerant.pcrit)), ΔTsup) * Vref *  (refrigerant.saturated_vapor.h(log10(tanhmap_R_to_ab(u[1], refrigerant.pr_min, refrigerant.pcrit)))  - hei) - u[2]   # heat balance refrigerant
    du[3] = tanhmap_R_to_01(u[4]) *Cmin_e*(u[8] - refrigerant.saturated_vapor.T(log10(tanhmap_R_to_ab(u[1], refrigerant.pr_min, refrigerant.pcrit))) ) - u[2]        # ϵ-NTU method 
    du[4] =  1 - exp(-Ue * A * tanhmap_R_to_01(u[6]) / Cmin_e)  -  tanhmap_R_to_01(u[4])

    # heat balance equations (superheater)
    du[5] = Cs_sup* (Ts_h - u[8]) - u[3]  # heat balance secondary fluid
    du[6] = Cref_sup * ΔTsup - u[3]
    du[7] = tanhmap_R_to_01(u[5]) * Cmin_sup * (Ts_h - refrigerant.saturated_vapor.T(log10(tanhmap_R_to_ab(u[1], refrigerant.pr_min, refrigerant.pcrit)))  ) - u[3]     # ϵ-NTU 
    du[8] = (1 - exp((-Usup*A*(1- tanhmap_R_to_01(u[6]))/Cmin_sup*(1-Cr))))/(1 - Cr *(exp(-Usup*A*(1- tanhmap_R_to_01(u[6]) )/Cmin_sup*(1-Cr)))) - tanhmap_R_to_01(u[5])  

    return nothing
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
`refrigerant::Refrigerant`, 
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
    pr = pe/refrigerant.pcrit    
    href_e = 1.5 * 55*pr^(0.12)/((-log10(pr))^0.55 * (refrigerant.M*1000)^0.5) * (abs(Qe)/(A * xs))^0.67

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
`refrigerant::Refrigerant`, 
`plateheatexchanger::PlateHeatExchanger`

return values

`Usup`, `hs_sup`, `href_sup`
"""

function evaporator_uvalue_superheat(Tsm_sup, pe, Vs, Vref, ΔTsup, secondaryfluid, refrigerant, plateheatexchanger)
    
    # heat_transfer_coefficient secondary fluid side
    ρs_sup = secondaryfluid.ρ(Tsm_sup)      
    cps_sup = secondaryfluid.cp(Tsm_sup)   
    μs_sup = secondaryfluid.μ(Tsm_sup)
    λs_sup = secondaryfluid.λ(Tsm_sup)
    properties_s_sup = ρs_sup, μs_sup, λs_sup, cps_sup

    hs_sup = h_platehex_wanniarachchi(Vs , plateheatexchanger, properties_s_sup)    

    # heat_transfer_coefficient refrigerant side
    ρref_sup  = refrigerant.superheated_vapor.ρ(log10(pe), ΔTsup)
    cpref_sup = refrigerant.superheated_vapor.cp(log10(pe), ΔTsup/2)
    μref_sup  = refrigerant.superheated_vapor.μ(log10(pe),ΔTsup/2)
    λref_sup  = refrigerant.superheated_vapor.λ(log10(pe),ΔTsup/2)
    propertiesR_sup = ρref_sup, μref_sup, λref_sup, cpref_sup
    
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
    pe = tanhmap_R_to_ab(pme, refrigerant.pr_min, refrigerant.pcrit)
    ϵe = tanhmap_R_to_01(ϵme)
    ϵsup = tanhmap_R_to_01(ϵmsup)
    xs  = tanhmap_R_to_01(xms)
    
    # two-phase
    Te = refrigerant.saturated_vapor.T(log10(pe))
    Tsm_e = (Ts_l + Ts_m)/2
    Ue, hs_e, href_e  =  evaporator_uvalue_twophase(Tsm_e, Vs, Qe, pe, xs, secondaryfluid, refrigerant, plateheatexchanger)
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
    hsup =  refrigerant.superheated_vapor.h(log10(pe),ΔTsup)     # fixed super heating of 5 degrees (find H for that super heating degree)
    ρref_sup  = refrigerant.superheated_vapor.ρ(log10(pe), ΔTsup)
    hes = refrigerant.saturated_vapor.h(log10(pe))
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

`pme`     refrigerant evaporation pressure mapped using tanh transformation, 
`Qe`      evaporator power exchange in the two phase region, 
`Qsup`    evaporator power exchange in the superheated region, 
`ϵme`     heat exchanger efficiency in the two phase region mapped using tanh transformation, 
`ϵmsup`   heat exchanger efficiency in the two phase region mapped using tanh transformation, 
`xms`     coordinate of the moving boundary mapped using tanh transformation, 
`Tsf_l`   secondary fluid temperature (in Kelvin) at the evaporator outlet, 
`Tsf_m`   secondary fluid temperature (in Kelvin) at the moving boundary, 
`Ts_h`    secondary fluid temperature (in Kelvin) at the inlet of the evaporator
"""
function evaporatorborehole_model!(du, u, p::EvaporatorBoreholeParameters)

    @unpack refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, borehole, λground, Tb = p
    A = plateheatexchanger.A
    ## u is a vector that contains the unknonws of the problem    
    ## u = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m, Tsf_h]

    # first we compute model variables that are considered locally "constant": 
    # These are Cs_e, Ue, Cs_sup, Cref_sup, Cr,  U_sup,

    # two-phase
    pe = tanhmap_R_to_ab(u[1], refrigerant.pr_min, refrigerant.pcrit)
    xs = tanhmap_R_to_01(u[6])
    Qe = u[2]
    Tsm_e = (u[7] + u[8])/2

    Ue, hs_e, href_e  =  evaporator_uvalue_twophase(Tsm_e, Vs, Qe, pe, xs, secondaryfluid, refrigerant, plateheatexchanger)
    Cs_e = Cmin_e = secondaryfluid.ρ(Tsm_e )* secondaryfluid.cp(Tsm_e) * Vs   #  capacity of secondary fluid 

    # superheat
    Tsm_sup = (u[9] + u[8])/2   
    Usup, hs_sup, href_sup = evaporator_uvalue_superheat(Tsm_sup, pe,  Vs, Vref, ΔTsup, secondaryfluid, refrigerant, plateheatexchanger)
    ρs_sup = secondaryfluid.ρ(Tsm_sup)      
    cps_sup = secondaryfluid.cp(Tsm_sup) 
    Cs_sup = ρs_sup * cps_sup * Vs
    hsup =  refrigerant.superheated_vapor.h(log10(pe),ΔTsup)     # fixed super heating of 5 degrees (find H for that super heating degree)
    ρref_sup  = refrigerant.superheated_vapor.ρ(log10(pe), ΔTsup)
    hes = refrigerant.saturated_vapor.h(log10(pe))
    Cref_sup = ρref_sup*Vref*(hsup - hes)/ΔTsup   #Cref_sup (mean capacity of the refrigerant in the superheating region)
    Cmin_sup = min(Cref_sup,Cs_sup)
    Cmax_sup = max(Cref_sup,Cs_sup)
    Cr = Cmin_sup/Cmax_sup    

    # Heat balance equations (two-phase)
    du[1] = Cs_e * (u[8] - u[7]) - u[2]  # heat balance secondary fluid
    du[2] = refrigerant.superheated_vapor.ρ(log10(tanhmap_R_to_ab(u[1], refrigerant.pr_min, refrigerant.pcrit)), ΔTsup) * Vref *  (refrigerant.saturated_vapor.h(log10(tanhmap_R_to_ab(u[1], refrigerant.pr_min, refrigerant.pcrit)))  - hei) - u[2]   # heat balance refrigerant
    du[3] = tanhmap_R_to_01(u[4]) *Cmin_e*(u[8] - refrigerant.saturated_vapor.T(log10(tanhmap_R_to_ab(u[1], refrigerant.pr_min, refrigerant.pcrit))) ) - u[2]        # ϵ-NTU method 
    du[4] =  1 - exp(-Ue * A * tanhmap_R_to_01(u[6]) / Cmin_e)  -  tanhmap_R_to_01(u[4])

    # heat balance equations (superheater)
    du[5] = Cs_sup* (u[9] - u[8]) - u[3]  # heat balance secondary fluid
    du[6] = Cref_sup * ΔTsup - u[3]
    du[7] = tanhmap_R_to_01(u[5]) * Cmin_sup * (u[9] - refrigerant.saturated_vapor.T(log10(tanhmap_R_to_ab(u[1], refrigerant.pr_min, refrigerant.pcrit)))  ) - u[3]     # ϵ-NTU 
    du[8] = (1 - exp((-Usup*A*(1- tanhmap_R_to_01(u[6]))/Cmin_sup*(1-Cr))))/(1 - Cr *(exp(-Usup*A*(1- tanhmap_R_to_01(u[6]) )/Cmin_sup*(1-Cr)))) - tanhmap_R_to_01(u[5])  #1/(1 + exp(-u[5])) # sigmoid_R_to_01(u[5])

    # borehole model 
    # kin = uniform_Tb_coeffs(bh, λground, secondaryfluid.ρ(Tsm_e )*Vs, (u[2] + u[9])/2, secondaryflui)
    # du[9] = kin*u[2] - 1*u[9] + (1-kin)* Tb

    du[9] = uniform_Tb_coeffs(borehole, λground, secondaryfluid.ρ(Tsm_e )*Vs, (u[7] + u[9]/2), secondaryfluid)*u[7] - u[9] + (1-uniform_Tb_coeffs(borehole, λground, secondaryfluid.ρ(Tsm_e )*Vs, (u[7] + u[9]/2), secondaryfluid))*Tb

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