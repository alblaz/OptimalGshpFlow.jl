#=
Bounded variables as types.

The models solve for unknowns u ∈ ℝ so that unconstrained nonlinear solvers can be used; each
physical variable x is recovered with `from_R(t, u)` and an initial guess is mapped with
`to_R(t, x)`, where `t` describes the admissible interval of x:

| transform          | x ∈        | x = from_R(t, u)                         | used for                      |
|--------------------|------------|------------------------------------------|-------------------------------|
| `Free()`           | ℝ          | u                                        | temperatures, heat flows      |
| `Linear(x0, s)`    | ℝ          | x0 + s u                                 | scaled heat flows             |
| `Unit()`           | (0, 1)     | σ(u)                                     | effectiveness, boundary x     |
| `Box(lo, hi)`      | (lo, hi)   | lo + (hi - lo) σ(u)                      | evaporation temperature       |
| `LogBox(lo, hi)`   | (lo, hi)   | exp(log lo + (log hi - log lo) σ(u))     | pressure                      |
| `Above(lo, scale)` | (lo, ∞)    | lo + scale exp(u)                        | molar volume (above co-volume) |

with σ(u) = 1/(1 + e⁻ᵘ). `LogBox` is a sigmoid in log x: a given step in u is a given relative
change of x, whatever its magnitude, which suits a pressure spanning several decades.
All maps are generic in the number type (ForwardDiff Duals pass through).
=#

"""
    VariableTransform

Map between a bounded physical variable `x` and an unbounded unknown `u ∈ ℝ`: [`Free`](@ref),
[`Linear`](@ref), [`Unit`](@ref), [`Box`](@ref), [`LogBox`](@ref), [`Above`](@ref). Interface: `from_R(t, u)`,
`to_R(t, x)`, `bounds(t)`.
"""
abstract type VariableTransform end

"""
    from_R(t::VariableTransform, u) -> x

Physical value of the unknown `u ∈ ℝ`; `x` lies inside `bounds(t)`.
"""
function from_R end

"""
    to_R(t::VariableTransform, x) -> u

Unknown in ℝ corresponding to the physical value `x` (inverse of [`from_R`](@ref)); `x` must lie
strictly inside `bounds(t)`.
"""
function to_R end

"""
    bounds(t::VariableTransform) -> (lo, hi)

Open interval of the physical variable.
"""
function bounds end

"No bound: `x = u`."
struct Free <: VariableTransform end
from_R(::Free, u) = u
to_R(::Free, x) = x
bounds(::Free) = (-Inf, Inf)

"""
    Linear(x0, scale)

No bound, scaled: `x = x0 + scale·u`, so that `u` is of order one when `x - x0` is of order
`scale` (e.g. a heat flow in units of a nominal duty).
"""
struct Linear{T <: Real} <: VariableTransform
    x0::T
    scale::T
    function Linear(x0::Real, scale::Real)
        x0, scale = promote(float(x0), float(scale))
        isfinite(x0) && 0 < scale < Inf || throw(ArgumentError("Linear needs finite x0 and 0 < scale < ∞, got ($x0, $scale)"))
        return new{typeof(x0)}(x0, scale)
    end
end
from_R(t::Linear, u) = t.x0 + t.scale * u
to_R(t::Linear, x) = (x - t.x0) / t.scale
bounds(t::Linear) = (oftype(t.x0, -Inf), oftype(t.x0, Inf))

"`x ∈ (0, 1)`: `x = σ(u)` (sigmoid)."
struct Unit <: VariableTransform end
from_R(::Unit, u) = sigmoid_R_to_01(u)
to_R(::Unit, x) = sigmoid_01_to_R(x)
bounds(::Unit) = (0.0, 1.0)

"""
    Box(lo, hi)

`x ∈ (lo, hi)`: `x = lo + (hi - lo) σ(u)` (linear scale).
"""
struct Box{T <: Real} <: VariableTransform
    lo::T
    hi::T
    function Box(lo::Real, hi::Real)
        lo, hi = promote(float(lo), float(hi))
        isfinite(lo) && isfinite(hi) && lo < hi || throw(ArgumentError("Box needs finite lo < hi, got ($lo, $hi)"))
        return new{typeof(lo)}(lo, hi)
    end
end
from_R(t::Box, u) = sigmoid_R_to_ab(u, t.lo, t.hi)
to_R(t::Box, x) = sigmoid_ab_to_R(x, t.lo, t.hi)
bounds(t::Box) = (t.lo, t.hi)

"""
    LogBox(lo, hi)

`x ∈ (lo, hi)` with `0 < lo`: sigmoid on `log x`, `x = exp(log lo + (log hi - log lo) σ(u))`.
"""
struct LogBox{T <: Real} <: VariableTransform
    lo::T
    hi::T
    loglo::T
    loghi::T
    function LogBox(lo::Real, hi::Real)
        lo, hi = promote(float(lo), float(hi))
        0 < lo < hi < Inf || throw(ArgumentError("LogBox needs 0 < lo < hi < ∞, got ($lo, $hi)"))
        return new{typeof(lo)}(lo, hi, log(lo), log(hi))
    end
end
from_R(t::LogBox, u) = exp(sigmoid_R_to_ab(u, t.loglo, t.loghi))
to_R(t::LogBox, x) = sigmoid_ab_to_R(log(x), t.loglo, t.loghi)
bounds(t::LogBox) = (t.lo, t.hi)

"""
    Above(lo, scale = lo)

`x ∈ (lo, ∞)`: `x = lo + scale exp(u)`, i.e. `log(x - lo)` on a scale set by `scale > 0`.
"""
struct Above{T <: Real} <: VariableTransform
    lo::T
    scale::T
    function Above(lo::Real, scale::Real = lo)
        lo, scale = promote(float(lo), float(scale))
        isfinite(lo) && 0 < scale < Inf || throw(ArgumentError("Above needs finite lo and 0 < scale < ∞, got ($lo, $scale)"))
        return new{typeof(lo)}(lo, scale)
    end
end
from_R(t::Above, u) = t.lo + t.scale * exp(u)
to_R(t::Above, x) = log((x - t.lo) / t.scale)
bounds(t::Above) = (t.lo, oftype(t.lo, Inf))
