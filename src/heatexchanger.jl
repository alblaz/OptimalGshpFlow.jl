###############################################
# extract heat exchanger temperature profiles #
###############################################

ϵ_counterflow(x,NTU,Cr) = (1 - exp(-NTU*(x)*(1-Cr)))/(1-Cr*exp(-NTU*(x)*(1-Cr)))

# the fluid with lower capacity Cc is flowing according to x (from 0 to 1), while the hot fluid with higher capacity is flowing on the contrary direction (from 1 to 0). This is important to keep in mind if we want to plot the inverse representation where the fluid with lower capacity enter the heat exchanger at x=1. In this case we first would perform the calculation above and then we would remap the the values of the

# function counterflow_heatexchanger(Ch,Thi,Tho,Cc,Tci,Tco,UA; x =0:0.1:1, defaultdirection = true) 
function counterflow_heatexchanger(heatexchanger_solution::HeatExchangerSolution; n=11 , defaultdirection = true) 

    x = range(0,1,n)
    @unpack Cc, Ch, NTU, Tch, Tcl, Thh, Thl, UA = heatexchanger_solution

    Thi, Tho = Thh, Thl
    Tci, Tco = Tcl, Tch

    Cmin = min(Cc,Ch)
    Cmax = max(Cc,Ch)
    # NTU = UA/Cmin
    Cr = Cmin/Cmax

    if Cc > Ch
        Th = @.  Thi - ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)
        Tc = @.  Tco - ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)*Ch/Cc
    else
        ϵx = ϵ_counterflow.(x,NTU,Cr)
        Tc = @.  Tci +  ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)
        Th = @.  Tho +  ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)*Cc/Ch
    end

    if defaultdirection == false
        reverse!(Tc) 
        reverse!(Th)
    end

    return x,Tc,Th
end
"""
    effectiveness_counterflow(NTU, γ, φ = 0)

Effectiveness of a counter-flow zone in which the refrigerant temperature changes along the flow
by the pressure drop and the glide (Qiao 2025, Appl. Therm. Eng. 279, 127896, eq. 56):

`ε = [1 + γ/((1−φ)NTU) − (γ + γ/((1−φ)NTU) + 1) e^{−(1−φ)NTU}] / (1 − φ e^{−(1−φ)NTU})`

with `ε = (T_f,out − T_f,in)/(T_r,in − T_f,in)` (`f` the single-phase fluid, `r` the refrigerant),
`NTU = U A / C_f`, `γ = ΔT_sat/(T_r,in − T_f,in)` (`ΔT_sat` the change of the refrigerant
saturation temperature over the zone caused by the pressure drop) and `φ = C_f/C_r` (glide; 0 for
a pure fluid). Limits: `γ = 0`, `φ = 0` → `1 − e^{−NTU}`; `φ = 0` → Qiao eq. 58;
`φ = 1` (balanced) → `NTU (1 + γ/2)/(1 + NTU)`. In an evaporator `γ > 0` and `ε` can exceed 1 (the
refrigerant enters warmer than it leaves).

Evaluated as `ε = [E + γ G] / (1/NTU + φ E)` with `a = (1 − φ) NTU`, `E = (1 − e^{−a})/a` and
`G = (E − e^{−a})/a` (eq. 56 divided by `a`; series near `a = 0`): finite and smooth through
`φ = 1`.
"""
function effectiveness_counterflow(NTU, γ, φ = zero(NTU))
    a = (1 - φ) * NTU
    E, G = _glide_EG(a)
    return (E + γ * G) / (1/NTU + φ * E)
end

# E = (1 − e^{−a})/a and G = (E − e^{−a})/a, with their series for |a| < 1e-3 (error ≈ a⁴/100)
function _glide_EG(a)
    if abs(a) < 1e-3
        return (1 - a/2 + a^2/6 - a^3/24, 1/2 - a/3 + a^2/8 - a^3/30)
    end
    e = exp(-a)
    E = -expm1(-a) / a
    return (E, (E - e) / a)
end
