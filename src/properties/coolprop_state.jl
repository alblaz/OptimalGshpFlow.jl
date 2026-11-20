#= Cached CoolProp `AbstractState` (HEOS): low-level CoolProp access shared by the refrigerant
backends (transport properties of `ClapeyronRefrigerant`, compressor and saturation pressures of
the spline `Refrigerant`). =#

"""
    CoolPropState(fluid)

Cached CoolProp `AbstractState` (HEOS backend) of a pure fluid, used for the transport properties
of `ClapeyronRefrigerant` and the CoolProp calls of the spline `Refrigerant`.
Updating it and reading outputs costs ≈ 3 µs, against ≈ 100 µs for one `PropsSI` call, with
identical results. A lock serialises access, so a refrigerant can be shared between threads.
"""
mutable struct CoolPropState
    fluid::String
    handle::Clong
    lock::ReentrantLock
    function CoolPropState(fluid::AbstractString)
        cs = new(String(fluid), CoolProp.AbstractState_factory("HEOS", fluid), ReentrantLock())
        finalizer(x -> CoolProp.AbstractState_free(x.handle), cs)
        return cs
    end
end
Base.show(io::IO, cs::CoolPropState) = print(io, "CoolPropState(", cs.fluid, ")")

const COOLPROP_PT = Ref{Clong}(-1)
const COOLPROP_PQ = Ref{Clong}(-1)
const COOLPROP_DT = Ref{Clong}(-1)
const COOLPROP_V  = Ref{Clong}(-1)
const COOLPROP_L  = Ref{Clong}(-1)
const COOLPROP_QT = Ref{Clong}(-1)
const COOLPROP_PS = Ref{Clong}(-1)
const COOLPROP_HP = Ref{Clong}(-1)
const COOLPROP_P  = Ref{Clong}(-1)
const COOLPROP_H  = Ref{Clong}(-1)
const COOLPROP_S  = Ref{Clong}(-1)
const COOLPROP_D  = Ref{Clong}(-1)
const COOLPROP_I  = Ref{Clong}(-1)    # surface tension
const COOLPROP_C  = Ref{Clong}(-1)    # isobaric specific heat (mass)

function _init_coolprop_indices()
    COOLPROP_PT[] = CoolProp.get_input_pair_index("PT_INPUTS")
    COOLPROP_PQ[] = CoolProp.get_input_pair_index("PQ_INPUTS")
    COOLPROP_DT[] = CoolProp.get_input_pair_index("DmassT_INPUTS")
    COOLPROP_V[]  = CoolProp.get_param_index("V")
    COOLPROP_L[]  = CoolProp.get_param_index("CONDUCTIVITY")
    COOLPROP_QT[] = CoolProp.get_input_pair_index("QT_INPUTS")
    COOLPROP_PS[] = CoolProp.get_input_pair_index("PSmass_INPUTS")
    COOLPROP_HP[] = CoolProp.get_input_pair_index("HmassP_INPUTS")
    COOLPROP_P[]  = CoolProp.get_param_index("P")
    COOLPROP_H[]  = CoolProp.get_param_index("Hmass")
    COOLPROP_S[]  = CoolProp.get_param_index("Smass")
    COOLPROP_D[]  = CoolProp.get_param_index("Dmass")
    COOLPROP_I[]  = CoolProp.get_param_index("I")
    COOLPROP_C[]  = CoolProp.get_param_index("Cpmass")
    return nothing
end

function _coolprop_value(cs::CoolPropState, pair, key, x1, x2)
    lock(cs.lock) do
        CoolProp.AbstractState_update(cs.handle, pair[], x1, x2)
        CoolProp.AbstractState_keyed_output(cs.handle, key[])
    end
end

# several outputs of one state update (one flash)
function _coolprop_values(cs::CoolPropState, pair, x1, x2, keys::Tuple)
    lock(cs.lock) do
        CoolProp.AbstractState_update(cs.handle, pair[], x1, x2)
        map(k -> CoolProp.AbstractState_keyed_output(cs.handle, k[]), keys)
    end
end

# CoolProp output for the input pair (x1, x2); for ForwardDiff.Dual inputs the derivatives are
# central differences with relative steps rel1, rel2
function _coolprop(cs::CoolPropState, pair, key, x1, x2, rel1, rel2)
    (x1 isa ForwardDiff.Dual || x2 isa ForwardDiff.Dual) || return _coolprop_value(cs, pair, key, x1, x2)
    a, b = ForwardDiff.value(x1), ForwardDiff.value(x2)
    f = _coolprop_value(cs, pair, key, a, b)
    D = x1 isa ForwardDiff.Dual ? typeof(x1) : typeof(x2)
    ∂ = zero(ForwardDiff.partials(one(D)))
    if x1 isa ForwardDiff.Dual
        h = rel1 * a
        ∂ += (_coolprop_value(cs, pair, key, a + h, b) - _coolprop_value(cs, pair, key, a - h, b)) / (2h) * ForwardDiff.partials(x1)
    end
    if x2 isa ForwardDiff.Dual
        h = rel2 * b
        ∂ += (_coolprop_value(cs, pair, key, a, b + h) - _coolprop_value(cs, pair, key, a, b - h)) / (2h) * ForwardDiff.partials(x2)
    end
    return D(f, ∂)
end
