#=
Validation of the Clapeyron/EntropyScaling refrigerant properties (`ClapeyronRefrigerant`)
against CoolProp, and timing compared with the spline-based `Refrigerant` and with PropsSI.

Run from an environment that has OptimalGshpFlow and CoolProp:
    include("examples/properties_validation.jl")
=#
using OptimalGshpFlow, CoolProp, ForwardDiff, Printf

relerr(a, b) = a / b - 1
# CoolProp reference value, NaN where CoolProp itself fails (e.g. R32 conductivity below ≈ -30 °C)
ref(args...) = try PropsSI(args...) catch; NaN end
pct(x) = @sprintf("%+7.2f%%", 100x)

#############################################
# 1. thermodynamic and transport deviations #
#############################################
println("Maximum deviation vs CoolProp, superheated vapour (ΔT = 2…20 K) and saturated states, pr = 0.02…0.7")
println(rpad("fluid", 8), rpad("h", 10), rpad("ρ", 10), rpad("cp", 10), rpad("μ (ES)", 10), rpad("λ (ES)", 10))
for f in ("R290", "R134a", "R32", "R1234yf")
    r = ClapeyronRefrigerant(f)
    e = Dict(k => 0.0 for k in (:h, :ρ, :cp, :μ, :λ))
    upd!(k, x) = !isnan(x) && abs(x) > abs(e[k]) && (e[k] = x)
    for pr in (0.02, 0.05, 0.1, 0.2, 0.4, 0.7)
        p = pr * r.pcrit
        Ts = T_dew(r, p)
        upd!(:h, relerr(h_dew(r, p), ref("H", "P", p, "Q", 1, f)))
        for ΔT in (2., 5., 10., 20.)
            T = Ts + ΔT
            upd!(:h,  relerr(h_pT(r, p, T),  ref("H", "P", p, "T", T, f)))
            upd!(:ρ,  relerr(ρ_pT(r, p, T),  ref("D", "P", p, "T", T, f)))
            upd!(:cp, relerr(cp_pT(r, p, T), ref("C", "P", p, "T", T, f)))
            upd!(:μ,  relerr(μ_pT(r, p, T; phase = :vapor), ref("V", "P", p, "T", T, f)))
            upd!(:λ,  relerr(λ_pT(r, p, T; phase = :vapor), ref("CONDUCTIVITY", "P", p, "T", T, f)))
        end
    end
    println(rpad(f, 8), join((rpad(@sprintf("%.1e", abs(e[k])), 10) for k in (:h, :ρ, :cp)), ""),
            rpad(pct(e[:μ]), 10), rpad(pct(e[:λ]), 10))
end

##################################
# 2. zeotropic mixture (R454C)   #
##################################
println("\nR454C (R32/R1234yf 21.5/78.5 wt%): Clapeyron vs CoolProp HEOS mixture")
mx = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785], name = "R454C")
fl = "HEOS::R32[$(mx.z[1])]&R1234yf[$(mx.z[2])]"
println(rpad("p [bar]", 9), rpad("Tbub [°C]", 20), rpad("Tdew [°C]", 20), rpad("glide [K]", 10), "Δh_lv rel. diff.")
for p in (2e5, 4.7e5, 10e5, 20e5)
    tb, td = T_bubble(mx, p), T_dew(mx, p)
    tbc, tdc = PropsSI("T", "P", p, "Q", 0, fl), PropsSI("T", "P", p, "Q", 1, fl)
    dh, dhc = h_dew(mx, p) - h_bubble(mx, p), PropsSI("H", "P", p, "Q", 1, fl) - PropsSI("H", "P", p, "Q", 0, fl)
    println(rpad(p / 1e5, 9), rpad(@sprintf("%.2f (CP %.2f)", tb - 273.15, tbc - 273.15), 20),
            rpad(@sprintf("%.2f (CP %.2f)", td - 273.15, tdc - 273.15), 20), rpad(@sprintf("%.2f", td - tb), 10), pct(relerr(dh, dhc)))
end

#############
# 3. timing #
#############
function timeit(f, n = 500)
    f()
    t = @elapsed for _ in 1:n
        f()
    end
    return t / n * 1e6   # µs
end
r   = ClapeyronRefrigerant("R290")
spl = try
    Refrigerant("R290")
catch       # spline construction for R290 needs the sampling fix (fix/refrigerant-pcrit-sampling)
    nothing
end
p   = 4.7e5
T   = T_dew(r, p) + 5
s   = s_pT(r, p, T)
println("\nTiming (µs per call), R290 at $(p/1e5) bar")
spl === nothing && println("(spline Refrigerant(\"R290\") not constructible on this version: skipped)")
spline_cases = spl === nothing ? () : (
        ("spline  h_dew (Refrigerant)",        () -> spl.saturated_vapor.h(log10(p))),
        ("spline  h superheated",              () -> spl.superheated_vapor.h(log10(p), 5.)))
for (label, f) in (spline_cases...,
        ("Clapeyron T_dew",                    () -> T_dew(r, p)),
        ("Clapeyron h_dew",                    () -> h_dew(r, p)),
        ("Clapeyron h_pT",                     () -> h_pT(r, p, T)),
        ("Clapeyron h_ps (compressor outlet)", () -> h_ps(r, 15e5, s)),
        ("EntropyScaling λ_pT",                () -> λ_pT(r, p, T)),
        ("Clapeyron T_dew (R454C mixture)",    () -> T_dew(mx, p)),
        ("CoolProp PropsSI h(p,T)",            () -> PropsSI("H", "P", p, "T", T, "R290")),
        ("ForwardDiff dh_dew/dp",              () -> ForwardDiff.derivative(q -> h_dew(r, q), p)))
    println(rpad(label, 38), @sprintf("%9.2f", timeit(f)))
end
