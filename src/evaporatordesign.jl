#=
#fix
# naming
# mapping constrained variables to -∞ -∞
# implement exponentials as interpolations between inlet outlet of evaporator
# volumetric efficiency as function of pressure ratio

# BASIC DESIGN OF EVAPORATOR (SUPER BASIC HYPOTHESIS)
# 1. nominal condensation temperature 40°C
# 2. full condensation, no subcooling
# 3. nominal evaporation temperature 0 °C
# 4. full evaporation, no superheating
# 5. nominal thermal power exchange in the evaporator 5000 W
# 6. other extra assumption, we assume a fix volume flow rate, we do not take into account volumetric efficiency of 1, etc.
=#

struct EvaporatorDesign{T <: Number}
    Qnom::T
    refrigerant::String
    Vref::T
    Te::T
    secondaryfluid::String
    Vs::T
    Thh::T
    Thl::T
    A::T
    U::T
    hei::T
    ΔT_sup::T
end

function EvaporatorDesign(Q2nom,Te,Tc, refrigerant,secondaryfluid; U =2000., ΔTs = 3., ΔT_l = 1.5, ΔT_sub=5., ΔT_sup=5., ηᵥ=0.95)
    A,Vref,Vsf,hei = evaporator_design(Q2nom,Te,Tc, refrigerant,secondaryfluid; U =U, ΔTsf = ΔTs, ΔT_l = ΔT_l, ΔT_sub=ΔT_sub, ΔT_sup=ΔT_sup, ηᵥ=ηᵥ)
    Thl = Te + ΔT_l
    Thh = Thl + ΔTs
    return EvaporatorDesign(Q2nom, refrigerant, Vref, Te, secondaryfluid, Vsf, Thh, Thl, A, U, hei, ΔT_sup)
end


"""
    `evaporator_design(Q2,Te,Tc, refrigerant,secondary_fluid; U =2000., ΔTs = 3., ΔT_l = 1.5)`

- `Q2`                  heat exchanged at evaporator, 
- `Te`                  evaporation temperature in Kelvin, 
- `Tc`                  condensation temperature in Kelvin, 
- `refrigerant`         a valid refrigerant from CoolProp,  
- `secondary_fluid`     a valid incompressible antifreeze mixture from CoolProp, 
- `U`                   design U-value for the heat exchanger (assumed), 
- `ΔTs`                 secondary fluid temperature difference accross the evaporator, 
- `ΔT_l`                pinch point on the low temperature side of the evaporator,
- `ΔT_sub`              subcooling,
- `ΔT_sup`              superheating,
- `ηᵥ=0.95              volumetric efficiency.

Evaluate the design parameters for the evaporator. The design parameter estimated are:
- `A`       the Area of the heat exchanger, 
- `Vref`    the volume flow rate of the refrigerant, 
- `Vsf`     the volume flow rate of the secondary fluid.
- `hei`     enthalpy of the refrigerant at the inlet of the evaporator


# Examples
```julia-repl
julia> using CoolProp
julia> A,Vref,Vsf,hei = evaporator_design(5000, 273.15-5, 273.15+40, "R134a", "INCOMP::MEA-28%")
```
"""
function evaporator_design(Q2,Te,Tc, refrigerant,secondary_fluid; U =2000., ΔTsf = 3., ΔT_l = 1.5, ΔT_sub=5., ΔT_sup=5., ηᵥ=0.95)
    
    hco = PropsSI("H", "T",  Tc-ΔT_sub, "Q", 0, refrigerant) #enthalpy at the outlet of the condenser - nominal conditions
    hei = hco                                        # enthalpy at the inlet of the evaporator - nominal conditions
    hes = PropsSI("H", "P", PropsSI("P", "T",  Te, "Q", 1, refrigerant), "T", Te+ΔT_sup, refrigerant) #enthalpy at the outlet of the evaporator - nominal design conditions
    ρeo = PropsSI("D", "P", PropsSI("P", "T",  Te, "Q", 1, refrigerant), "T", Te+ΔT_sup, refrigerant) #density at the outlet of the evaporator - nominal design conditions

    #The density is used to evaluate the mass flow rate of the refrigerant simply by multiplying by volume flow rate of the refrigerant
    ## VR volume flow rate of refrigerant
    Vref  = Q2/((hes - hei)*ρeo)*ηᵥ

    ## Vs secondary fluid volume flow rate assuming a ΔTs = 3 across the evaporator. I assumed water properties for the secondary fluid
    cps = PropsSI("C","T",Te+ΔT_l + ΔTsf/2,"P",101325,secondary_fluid)
    ρs  = PropsSI("D","T",Te+ΔT_l + ΔTsf/2,"P",101325,secondary_fluid)
    Vsf = Q2/(ΔTsf*ρs*cps)

    ## We assume a U value (default 2000 W/m^2), and a value for ΔT_low  and we compute the area A

    # U = 2000  # W/m^2
    Tsf_l = Te + ΔT_l
    Tsf_h = Tsf_l + ΔTsf
    ΔT_h = Tsf_h - Te
    # LMTD = (ΔT1 - ΔT2)/log(ΔT1/ΔT2)
    A = Q2/U * log(ΔT_h/ΔT_l)/ (ΔT_h - ΔT_l)

    return A,Vref,Vsf,hei
end