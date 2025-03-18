
"""
  struct SecondaryFluid

  Fields\\
  ≡≡≡≡≡≡

  ```
  ρ  :: Real,
  μ  :: Real,
  cp :: Real,
  λ  :: Real,
  ```
"""
struct SecondaryFluid
    name::String
    ρ
    μ
    cp
    λ
 end

"""
Constructor\\
≡≡≡≡≡≡≡≡≡≡≡

    `SecondaryFluid(Trange, fluid_name)`
 
- `Trange` temperature range used to build the interpolatants of the secondary properties
- `fluid_name` is a string with a valid secondary fluid CoolProp name.
"""
function SecondaryFluid(Trange, fluid_name)
    spl = secondary_fluid(Trange, fluid_name)
    return SecondaryFluid(fluid_name, spl.ρ, spl.μ, spl.cp, spl.λ)
end

"""
`secondary_fluid(Trange, secondaryfluid_name; ns = 40 )`

`Trange` is the range of temperature, `secondaryfluid_name` is a string containing an incompressible fluid supported by CoolProp, `ns=40` is the number of sample points for the interpolation. 

Example
```julia-repl
Trange = range(-15. +  273.15, 30 + 273.15, length = ns)
ethanol =  "INCOMP::MEA-28%"
spl_sf = secondary_fluid(Trange, ethanol)
spl_sf.cp(273.15)    #specific heat capacity
```
"""
function secondary_fluid(Trange, secondaryfluid_name)
    ρs_T = map(x -> PropsSI("D","T",x,"P",101325,secondaryfluid_name), Trange)
    μs_T = map(x -> PropsSI("V","T",x,"P",101325,secondaryfluid_name), Trange)
    cps_T = map(x -> PropsSI("C","T",x,"P",101325,secondaryfluid_name), Trange)
    λ_T  = map(x -> PropsSI("CONDUCTIVITY","T",x,"P",101325,secondaryfluid_name), Trange)

    spl_ρs  = cubic_spline_interpolation(Trange, ρs_T, extrapolation_bc = Flat())
    spl_μs  = cubic_spline_interpolation(Trange, μs_T, extrapolation_bc = Flat())
    spl_cps = cubic_spline_interpolation(Trange, cps_T, extrapolation_bc = Flat())
    spl_λs  = cubic_spline_interpolation(Trange, λ_T, extrapolation_bc = Flat())

    spl = (ρ=spl_ρs, μ = spl_μs, cp=spl_cps, λ= spl_λs)
    return spl
end


struct RefrigerantSaturatedVapor
    h
    T
    ρ
    cp
    μ
    λ
 end

 struct RefrigerantSuperHeatedVapor
    h
    ρ
    cp
    μ
    λ
 end

"""
  struct Refrigerant

  Fields\\
  ≡≡≡≡≡≡
  
  - saturated_vapor::RefrigerantSaturatedVapor
  - superheated_vapor::RefrigerantSuperHeatedVapor
  - `pcrit`: critical pressure

"""
struct Refrigerant
    name::String
    saturated_vapor::RefrigerantSaturatedVapor
    superheated_vapor::RefrigerantSuperHeatedVapor
    pcrit
    pr_min
    M
end

"""
Constructor\\
≡≡≡≡≡≡≡≡≡≡≡

    `Refrigerant(refrigerant; superheating_range = 0:1.:40., ns = 100, pr_min = 1e-4)`

Arguments
=========

- `refrigerant::String` temperature range used to build the interpolatants of the secondary properties\\

Keyword arguments
=================

- `superheating_range = 0:1.:40``: range of temperatures used to generate the interpolants
- `ns=100`: number of points used on the log10 pressure axis used 
- `pr_min = 1e-4`: minimum pressure value used for generating the interpolants
"""
function Refrigerant(refrigerant_name; superheating_range = 0:1.:40., ns = 100, pr_min = 1e-4)
    saturated = refrigerant_twophaseregion(refrigerant_name, 1.; ns = 100, pr_min = 1e-4)
    superheated = refrigerant_superheated(refrigerant_name; superheating_range = 0:1.:40., ns = 100, pr_min = 1e-4)
    pcrit = PropsSI.("Pcrit",  refrigerant_name) # critical pressure 
    M = PropsSI.("M",  refrigerant_name)

    return Refrigerant(refrigerant_name, RefrigerantSaturatedVapor(saturated.h, saturated.T, saturated.ρ, saturated.cp, saturated.μ,  saturated.λ),
    RefrigerantSuperHeatedVapor(superheated.h, superheated.ρ, superheated.cp, superheated.μ,  superheated.λ),
    pcrit,
    pr_min,
    M)
end

"""
`refrigerant_twophaseregion(refrigerant, x; ns = 100, pr_min = 1e-4)`

`refrigerant` is a string containing a refrigerant supported by CoolProp, x is the vapor quality, `ns=100` is the number of sample point utilized to create the interpolant and `pr_min` is the reduced pressure used as lower bound for the interpolation.

One dimensional interpolants of relevant properties and state variables of a refrigerant (supported by CoolProp) along iso-vapor-quality lines as function of `log10(p)` where p is pressure. The quantities included are enthalpy (`h`), temperature (`T`), density (`ρ`), specific heat capacity (`cp`) and visccosity (`μ`).

Example:
```julia-repls
spl_sat = refrigerant_twophaseregion("R134a", 1) 
p = 10000.
spl_sat.h(log10(p))  #  enthalpy
spl_sat.T(log10(p))  #  temperature 
spl_sat.ρ(log10(p))  #  density
spl_sat.cp(log10(p)) #  specific heat capacity
spl_sat.μ(log10(p))  #  viscosity
spl_sat.λ(log10(p))  #  conductivity
```
"""
function refrigerant_twophaseregion(refrigerant, x; ns = 100, pr_min = 1e-4)
    # sample points from the p-h Diagram
    pcrit = PropsSI.("Pcrit",  refrigerant)

    pmin,pmax = pr_min*pcrit, pcrit
    ps = map(exp10, range(log10(pmin),log10(pcrit), ns) )
    h1 = map(p-> PropsSI.("H", "P",  p, "Q", x, refrigerant), ps)
    ρ1 = map(p-> PropsSI.("D", "P",  p, "Q", x, refrigerant), ps)
    T1 = map(p-> PropsSI.("T", "P",  p, "Q", x, refrigerant), ps)
    cp1 = map(p-> PropsSI.("C", "P",  p, "Q", x, refrigerant), ps)
    μ1 = map(p-> PropsSI.("V", "P",  p, "Q", x, refrigerant), ps)
    λ1 = map(p-> PropsSI.("CONDUCTIVITY", "P",  p, "Q", x, refrigerant), ps)
    spl_hsat = cubic_spline_interpolation(range(log10(pmin),log10(pmax), ns), h1, extrapolation_bc = Flat())
    spl_Tref = cubic_spline_interpolation(range(log10(pmin),log10(pmax), ns), T1, extrapolation_bc = Flat())
    spl_ρref = cubic_spline_interpolation(range(log10(pmin),log10(pmax), ns), ρ1, extrapolation_bc = Flat())
    spl_cpref = cubic_spline_interpolation(range(log10(pmin),log10(pmax), ns), cp1, extrapolation_bc = Flat())
    spl_μref = cubic_spline_interpolation(range(log10(pmin),log10(pmax), ns), μ1, extrapolation_bc = Flat())
    spl_λref = cubic_spline_interpolation(range(log10(pmin),log10(pmax), ns), λ1, extrapolation_bc = Flat())

    spl_sat = (h = spl_hsat, T = spl_Tref, ρ = spl_ρref, cp = spl_cpref, μ = spl_μref, λ=spl_λref)
    return spl_sat
end


# build spline for superheated region
"""

`refrigerant_superheatedregion(refrigerant; superheating_range = 0:1.:40., ns = 100, pr_min = 1e-4)`

`refrigerant` is a string containing a refrigerant supported by CoolProp, `superheating_range = 0:1.:40.` is the superheating range for the interpolation, `ns=100` is the number of point utilized to sample the pressure to create the interpolant and `pr_min` is the reduced pressure used as lower bound for the interpolation.

Compute 2D interpolants of relevant properties and state variables of a refrigerant (supported by CoolProp) in the superheated region as function of `log10(p)` where p is pressure and of the degrees of superheat. The quantities included are enthalpy (`h`), temperature (`T`), density (`ρ`), specific heat capacity (`cp`),  visccosity (`μ`) and conductivity (`λ`). 

Example:
```
spl = refrigerant_superheated("R134a") 
p, ΔTsup = 10000., 3.
spl.h(log10(p), ΔTsup)  #  enthalpy
spl.ρ(log10(p), ΔTsup)  #  density
spl.cp(log10(p),ΔTsup)  #  specific heat capacity
spl.μ(log10(p), ΔTsup)  #  visccosity
spl.λ(log10(p), ΔTsup)  #  conductivity
```
"""
function refrigerant_superheated(refrigerant; superheating_range = 0:1.:40., ns = 100, pr_min = 1e-4)
    
    pcrit = PropsSI.("Pcrit",  refrigerant)
    pmin,pmax = pr_min*pcrit, pcrit
    ps = map(exp10, range(log10(pmin),log10(pcrit), ns) )

    h1 = map(p-> PropsSI.("H", "P",  p, "Q", 1., refrigerant), ps)
    T1 = map(p-> PropsSI.("T", "P",  p, "Q", 1., refrigerant), ps)
    ρ1 = map(p-> PropsSI.("D", "P",  p, "Q", 1., refrigerant), ps)
    μ1 = map(p-> PropsSI.("V", "P",  p, "Q", 1., refrigerant), ps)
    λ1 = map(p-> PropsSI.("CONDUCTIVITY", "P",  p, "Q", 1., refrigerant), ps)
    cp1 = map(p-> PropsSI.("C", "P",  p, "Q", 1., refrigerant), ps)

    C = [(p...,s) for p in zip(ps,T1), s in superheating_range]

    B_H = similar(ps*superheating_range')
    B_H[:,1] .= h1
    B_H[:, 2:end] .= map(x-> PropsSI.("H", "P",  x[1], "T", x[2] + x[3], refrigerant), C[:,2:end] )

    B_ρ = similar(B_H)
    B_ρ[:,1] .= ρ1
    B_ρ[:, 2:end] .= map(x-> PropsSI.("D", "P",  x[1], "T", x[2] + x[3], refrigerant), C[:,2:end] )

    B_μ = similar(B_H)
    B_μ[:,1] .= μ1
    B_μ[:, 2:end] .= map(x-> PropsSI.("V", "P",  x[1], "T", x[2] + x[3], refrigerant), C[:,2:end])

    B_cp = similar(B_H)
    B_cp[:,1] .= cp1
    B_cp[:, 2:end] .= map(x-> PropsSI.("C", "P",  x[1], "T", x[2] + x[3], refrigerant), C[:,2:end] )

    B_λ = similar(B_H)
    B_λ[:,1] .= λ1
    B_λ[:, 2:end] .= map(x-> PropsSI.("CONDUCTIVITY", "P",  x[1], "T", x[2] + x[3], refrigerant), C[:,2:end] )


    splpsΔT_H = cubic_spline_interpolation( (range(log10(pmin),log10(pmax), ns), superheating_range ) ,B_H, extrapolation_bc = Flat())
    splpsΔT_ρ = cubic_spline_interpolation( (range(log10(pmin),log10(pmax), ns), superheating_range ) ,B_ρ, extrapolation_bc = Flat())
    splpsΔT_μ = cubic_spline_interpolation( (range(log10(pmin),log10(pmax), ns), superheating_range ) ,B_μ, extrapolation_bc = Flat())
    splpsΔT_cp = cubic_spline_interpolation( (range(log10(pmin),log10(pmax), ns), superheating_range ) ,B_cp, extrapolation_bc = Flat())
    splpsΔT_λ = cubic_spline_interpolation( (range(log10(pmin),log10(pmax), ns), superheating_range ) ,B_λ, extrapolation_bc = Flat())


    spl = ( h =splpsΔT_H, ρ = splpsΔT_ρ, μ=splpsΔT_μ, cp = splpsΔT_cp, λ =  splpsΔT_λ  )
    return spl
end