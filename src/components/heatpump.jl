#=
Closed heat-pump cycle (step 9d): evaporator, compressor, condenser and borehole as one
`ComponentSystem`, with the standard node names, and an initial guess from a few temperatures.
=#

"""
    heatpump_cycle(; evaporator, compressor, condenser, borehole, Tb, Tw_in) -> ComponentSystem

Closed cycle: a cycle [`Evaporator`](@ref) (`hei = nothing`), a [`Compressor`](@ref), a
[`Condenser`](@ref) and a [`BoreholeHX`](@ref), all for the same refrigerant. The expansion valve is
isenthalpic: the condenser outlet and the evaporator inlet share the node `h_liq`.

Nodes: `T_cold`, `T_warm` (brine), `Tb` (borehole wall, fixed; not the undisturbed ground `Tg`), `p_suc`, `T_suc` (evaporator outlet), `m`
(refrigerant mass flow), `pc` (condensing pressure), `h_dis` (discharge enthalpy), `h_liq`
(liquid line), `Tw_in` (fixed), `Tw_out` (water). Guess with [`heatpump_guess`](@ref); change the
borehole wall or water temperature with [`set_fixed`](@ref).
"""
function heatpump_cycle(; evaporator::Evaporator, compressor::Compressor, condenser::Condenser, borehole::BoreholeHX,
                        Tb, Tw_in)
    evaporator.hei === nothing ||
        throw(ArgumentError("heatpump_cycle: the evaporator needs `hei = nothing` (inlet enthalpy from the condenser)"))
    return ComponentSystem((evaporator = evaporator, compressor = compressor, condenser = condenser, borehole = borehole);
        nodes = (T_cold = (:evaporator => :Ts_out, :borehole => :Tin),
                 T_warm = (:evaporator => :Ts_in, :borehole => :Tout),
                 Tb = (:borehole => :Tb,),
                 p_suc = (:evaporator => :p_suc, :compressor => :p_suc),
                 T_suc = (:evaporator => :T_suc, :compressor => :T_suc),
                 m = (:evaporator => :m, :compressor => :m, :condenser => :m),
                 pc = (:compressor => :pc, :condenser => :pc),
                 h_dis = (:compressor => :h_dis, :condenser => :h_in),
                 h_liq = (:condenser => :h_out, :evaporator => :hei),
                 Tw_in = (:condenser => :Tw_in,),
                 Tw_out = (:condenser => :Tw_out,)),
        fixed = (Tb = Tb, Tw_in = Tw_in))
end

"""
    heatpump_guess(sys; Te, Tc, T_cold, T_warm, Tw_out) -> u0

Initial guess of a [`heatpump_cycle`](@ref) from the dew temperatures of evaporation `Te` and
condensation `Tc`, and the brine and water temperatures: the refrigerant states, mass flow and
zone duties follow from the cycle (compression at the guessed pressures, saturated liquid).
"""
function heatpump_guess(sys::ComponentSystem; Te, Tc, T_cold, T_warm, Tw_out)
    (; evaporator, compressor, condenser) = sys.components
    r = evaporator.refrigerant
    pe, pc = _p_dew(r, Te), _p_dew(r, Tc)
    Tsuc = Te + evaporator.ΔTsup
    m = ρ_pT(r, pe, Tsuc; phase = :vapor) * evaporator.Vref
    h_dis = compression_state(r, pe, Tsuc, compressor.Δp_suction, pc, compressor.ηis).h_out
    h_liq, hd = h_bubble(r, pc), h_dew(r, pc)
    Qe, Qsup = m * (h_dew(r, pe) - h_liq), m * (h_pT(r, pe, Tsuc; phase = :vapor) - h_dew(r, pe))
    Qd, Qc = m * (h_dis - hd), m * (hd - h_liq)
    Tw_in = sys.nodefixed[findfirst(==(:Tw_in), _nodenames(sys))]
    return initial_guess(sys;
        evaporator = (Te = Te, Qe = Qe, Qsup = Qsup, Ts_m = T_warm - (T_warm - T_cold) * Qsup / (Qe + Qsup)),
        condenser = (Qd = Qd, Qc = Qc, Tw_m = Tw_out - (Tw_out - Tw_in) * Qd / (Qd + Qc)),
        nodes = (T_cold = T_cold, T_warm = T_warm, p_suc = pe, T_suc = Tsuc, m = m, pc = pc, h_dis = h_dis,
                 h_liq = h_liq, Tw_out = Tw_out))
end
