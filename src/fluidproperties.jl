
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
    s
 end

 struct RefrigerantSuperHeatedVapor
    h
    ρ
    cp
    μ
    λ
    s
    # his
 end

"""
  struct Refrigerant

  Fields\\
  ≡≡≡≡≡≡
  
  - saturated_vapor::RefrigerantSaturatedVapor
  - superheated_vapor::RefrigerantSuperHeatedVapor
  - `pcrit`: critical pressure
  - `pr_min`: minimum reduced pressure of the sampling; `M`: molar mass [kg/mol]
  - `pmin`, `pmax`: sampled pressure range [Pa] (the [`pressure_bounds`](@ref))

"""
struct Refrigerant <: AbstractRefrigerant
    name::String
    saturated_vapor::RefrigerantSaturatedVapor
    superheated_vapor::RefrigerantSuperHeatedVapor
    pcrit
    pr_min
    M
    pmin    # sampled pressure range [Pa], see log10p_grid
    pmax
end

"""
Constructor\\
≡≡≡≡≡≡≡≡≡≡≡

    `Refrigerant(refrigerant; superheating_range = 0:1.:40., ns = 100, pr_min = 1e-2, pr_max = 0.95)`

Arguments
=========

- `refrigerant::String` temperature range used to build the interpolatants of the secondary properties\\

Keyword arguments
=================

- `superheating_range = 0:1.:40``: range of temperatures used to generate the interpolants
- `ns=100`: number of points used on the log10 pressure axis used 
- `pr_min = 1e-2`: minimum reduced pressure p/pcrit used for generating the interpolants
  (clamped to just above the triple point pressure). Some fluids need a larger value because
  CoolProp transport properties fail at low pressure, e.g. `Refrigerant("R32"; pr_min = 0.05)`.
- `pr_max = 0.95`: maximum reduced pressure p/pcrit used for generating the interpolants.
  It must stay below 1: at the critical point cp, μ and λ diverge and the cubic spline rings.
"""
function Refrigerant(refrigerant_name; superheating_range = 0:1.:40., ns = 100, pr_min = 1e-2, pr_max = 0.95)
    saturated = refrigerant_twophaseregion(refrigerant_name, 1.; ns = ns, pr_min = pr_min, pr_max = pr_max)
    superheated = refrigerant_superheated(refrigerant_name; superheating_range = superheating_range, ns = ns, pr_min = pr_min, pr_max = pr_max)
    pcrit = PropsSI.("Pcrit",  refrigerant_name) # critical pressure 
    M = PropsSI.("M",  refrigerant_name)
    _, ps = log10p_grid(refrigerant_name; ns = ns, pr_min = pr_min, pr_max = pr_max)

    return Refrigerant(refrigerant_name, RefrigerantSaturatedVapor(saturated.h, saturated.T, saturated.ρ, saturated.cp, saturated.μ,  saturated.λ, saturated.s),
    RefrigerantSuperHeatedVapor(superheated.h, superheated.ρ, superheated.cp, superheated.μ, superheated.λ, superheated.s),
    pcrit,
    pr_min,
    M,
    first(ps),
    last(ps))
end

"""
`log10p_grid(refrigerant; ns = 100, pr_min = 1e-2, pr_max = 0.95)`

Uniform grid in `log10(p)` used to sample `refrigerant` properties, from
`max(pr_min*pcrit, 1.01*ptriple)` to `pr_max*pcrit`. Returns the range of `log10(p)`
and the corresponding pressures. `pr_max < 1` keeps the grid away from the critical
point, where CoolProp fails or properties diverge.
"""
function log10p_grid(refrigerant; ns = 100, pr_min = 1e-2, pr_max = 0.95)
    0 < pr_min < pr_max < 1 || throw(ArgumentError("need 0 < pr_min < pr_max < 1, got pr_min=$pr_min, pr_max=$pr_max"))
    pcrit = PropsSI("Pcrit", refrigerant)
    ptriple = try PropsSI("ptriple", refrigerant) catch; 0. end
    pmin = max(pr_min * pcrit, 1.01 * ptriple)
    pmax = pr_max * pcrit
    log10ps = range(log10(pmin), log10(pmax), ns)
    ps = exp10.(log10ps)
    return log10ps, ps
end

"""
`refrigerant_twophaseregion(refrigerant, x; ns = 100, pr_min = 1e-2, pr_max = 0.95)`

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
function refrigerant_twophaseregion(refrigerant, x; ns = 100, pr_min = 1e-2, pr_max = 0.95)
    # sample points from the p-h Diagram
    log10ps, ps = log10p_grid(refrigerant; ns = ns, pr_min = pr_min, pr_max = pr_max)
    h1 = map(p-> PropsSI.("H", "P",  p, "Q", x, refrigerant), ps)
    ρ1 = map(p-> PropsSI.("D", "P",  p, "Q", x, refrigerant), ps)
    T1 = map(p-> PropsSI.("T", "P",  p, "Q", x, refrigerant), ps)
    cp1 = map(p-> PropsSI.("C", "P",  p, "Q", x, refrigerant), ps)
    μ1 = map(p-> PropsSI.("V", "P",  p, "Q", x, refrigerant), ps)
    λ1 = map(p-> PropsSI.("CONDUCTIVITY", "P",  p, "Q", x, refrigerant), ps)
    s1 = map(p-> PropsSI.("S", "P",  p, "Q", x, refrigerant), ps)
    spl_hsat = cubic_spline_interpolation(log10ps, h1, extrapolation_bc = Flat())
    spl_Tref = cubic_spline_interpolation(log10ps, T1, extrapolation_bc = Flat())
    spl_ρref = cubic_spline_interpolation(log10ps, ρ1, extrapolation_bc = Flat())
    spl_cpref = cubic_spline_interpolation(log10ps, cp1, extrapolation_bc = Flat())
    spl_μref = cubic_spline_interpolation(log10ps, μ1, extrapolation_bc = Flat())
    spl_λref = cubic_spline_interpolation(log10ps, λ1, extrapolation_bc = Flat())
    spl_sref = cubic_spline_interpolation(log10ps, s1, extrapolation_bc = Flat())

    spl_sat = (h = spl_hsat, T = spl_Tref, ρ = spl_ρref, cp = spl_cpref, μ = spl_μref, λ=spl_λref, s = spl_sref)
    return spl_sat
end


# build spline for superheated region
"""

`refrigerant_superheated(refrigerant; superheating_range = 0:1.:40., ns = 100, pr_min = 1e-2, pr_max = 0.95)`

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
function refrigerant_superheated(refrigerant; superheating_range = 0:1.:40., ns = 100, pr_min = 1e-2, pr_max = 0.95)
    log10ps, ps = log10p_grid(refrigerant; ns = ns, pr_min = pr_min, pr_max = pr_max)

    h1 = map(p-> PropsSI.("H", "P",  p, "Q", 1., refrigerant), ps)
    T1 = map(p-> PropsSI.("T", "P",  p, "Q", 1., refrigerant), ps)
    ρ1 = map(p-> PropsSI.("D", "P",  p, "Q", 1., refrigerant), ps)
    μ1 = map(p-> PropsSI.("V", "P",  p, "Q", 1., refrigerant), ps)
    λ1 = map(p-> PropsSI.("CONDUCTIVITY", "P",  p, "Q", 1., refrigerant), ps)
    cp1 = map(p-> PropsSI.("C", "P",  p, "Q", 1., refrigerant), ps)
    s1 = map(p-> PropsSI.("S", "P",  p, "Q", 1., refrigerant), ps)

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

    B_s = similar(B_H)
    B_s[:,1] .= s1
    B_s[:, 2:end] .= map(x-> PropsSI.("S", "P",  x[1], "T", x[2] + x[3], refrigerant), C[:,2:end] )

    splpsΔT_H = cubic_spline_interpolation( (log10ps, superheating_range ) ,B_H, extrapolation_bc = Flat())
    splpsΔT_ρ = cubic_spline_interpolation( (log10ps, superheating_range ) ,B_ρ, extrapolation_bc = Flat())
    splpsΔT_μ = cubic_spline_interpolation( (log10ps, superheating_range ) ,B_μ, extrapolation_bc = Flat())
    splpsΔT_cp = cubic_spline_interpolation( (log10ps, superheating_range ) ,B_cp, extrapolation_bc = Flat())
    splpsΔT_λ = cubic_spline_interpolation( (log10ps, superheating_range ) ,B_λ, extrapolation_bc = Flat())
    splpsΔT_s = cubic_spline_interpolation( (log10ps, superheating_range ) ,B_s, extrapolation_bc = Flat())


    spl = ( h =splpsΔT_H, ρ = splpsΔT_ρ,  cp = splpsΔT_cp,  μ=splpsΔT_μ, λ =  splpsΔT_λ, s =  splpsΔT_s  )
    return spl
end

##################################################
# property interface (see properties/interface.jl) #
##################################################

# The spline backend covers saturated and superheated vapour of a pure fluid only.
# Properties are interpolated in log10(p) (and in the superheat ΔT); outside the sampled range
# the splines extrapolate flat.

# the sampled range: outside it the splines extrapolate flat (zero derivative in p)
pressure_bounds(r::Refrigerant) = (r.pmin, r.pmax)

T_dew(r::Refrigerant, p)  = r.saturated_vapor.T(log10(p))
T_bubble(r::Refrigerant, p) = T_dew(r, p)               # pure fluid
h_dew(r::Refrigerant, p)  = r.saturated_vapor.h(log10(p))
ρ_dew(r::Refrigerant, p)  = r.saturated_vapor.ρ(log10(p))
s_dew(r::Refrigerant, p)  = r.saturated_vapor.s(log10(p))
cp_dew(r::Refrigerant, p) = r.saturated_vapor.cp(log10(p))
μ_dew(r::Refrigerant, p)  = r.saturated_vapor.μ(log10(p))
λ_dew(r::Refrigerant, p)  = r.saturated_vapor.λ(log10(p))

h_superheated(r::Refrigerant, p, ΔT)  = r.superheated_vapor.h(log10(p), ΔT)
s_superheated(r::Refrigerant, p, ΔT)  = r.superheated_vapor.s(log10(p), ΔT)
ρ_superheated(r::Refrigerant, p, ΔT)  = r.superheated_vapor.ρ(log10(p), ΔT)
cp_superheated(r::Refrigerant, p, ΔT) = r.superheated_vapor.cp(log10(p), ΔT)
μ_superheated(r::Refrigerant, p, ΔT)  = r.superheated_vapor.μ(log10(p), ΔT)
λ_superheated(r::Refrigerant, p, ΔT)  = r.superheated_vapor.λ(log10(p), ΔT)

# (p, T) access: vapour only, T at or above the dew point (within the sampled superheat range)
function _spline_superheat(r::Refrigerant, p, T, phase)
    phase === :liquid && throw(ArgumentError("the spline Refrigerant has no liquid-phase properties"))
    return T - T_dew(r, p)
end
h_pT(r::Refrigerant, p, T; phase = :vapor)  = h_superheated(r, p, _spline_superheat(r, p, T, phase))
s_pT(r::Refrigerant, p, T; phase = :vapor)  = s_superheated(r, p, _spline_superheat(r, p, T, phase))
ρ_pT(r::Refrigerant, p, T; phase = :vapor)  = ρ_superheated(r, p, _spline_superheat(r, p, T, phase))
cp_pT(r::Refrigerant, p, T; phase = :vapor) = cp_superheated(r, p, _spline_superheat(r, p, T, phase))
μ_pT(r::Refrigerant, p, T; phase = :vapor)  = μ_superheated(r, p, _spline_superheat(r, p, T, phase))
λ_pT(r::Refrigerant, p, T; phase = :vapor)  = λ_superheated(r, p, _spline_superheat(r, p, T, phase))
