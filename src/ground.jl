#= Ground response: borehole wall temperature from the history of the heat extracted.

First dynamic ground model (ROADMAP step 8): infinite line source (ILS) with temporal
superposition of piecewise-constant loads. Not accurate over several years (no axial effects, no
borehole interaction); enough to drive the steady-state heat pump through a heating season.
=#

"""
    ILS(Fo)
    ILS(t, r, α)

Dimensionless temperature response of the infinite line source, `½ E₁(1/(4 Fo))`, at the Fourier
number `Fo = α t / r²`: a constant heat rate `q` [W/m] per unit length switched on at `t = 0`
changes the temperature at the distance `r` by `q/(2πλ) ILS(α t / r²)`.
"""
ILS(Fo::Number) = 1/2 * expint(1/(4Fo))
ILS(t::T, r::T, α::T) where {T <: Number} = ILS((α*t)/r^2)
ILS(t::Number, r::Number, α::Number) = ILS(promote(t, r, α)...)

"""
    InfiniteLineSource(; λ, α, rb, Tg)

Ground around a borehole as an infinite line source: conductivity `λ` [W/(m K)], thermal
diffusivity `α` [m²/s], borehole radius `rb` [m] (where the wall temperature is taken) and
undisturbed ground temperature `Tg` [K]. See [`borehole_wall_temperature`](@ref).
"""
struct InfiniteLineSource{T <: Real}
    λ::T
    α::T
    rb::T
    Tg::T
end
InfiniteLineSource(; λ, α, rb, Tg) = InfiniteLineSource(promote(float(λ), float(α), float(rb), float(Tg))...)

"""
    borehole_wall_temperature(g::InfiniteLineSource, tq, q, t) -> Tb

Borehole wall temperature [K] at the time `t` [s] for a piecewise-constant heat extraction: `q[i]`
[W per metre of borehole, > 0 when heat is extracted] from the time `tq[i]` on (`tq` increasing;
before `tq[1]` no load). Temporal superposition of the load steps:

`Tb = Tg − Σᵢ (q[i] − q[i-1]) / (2πλ) ⋅ ILS(t − tq[i], rb, α)` over the steps with `tq[i] < t`.

Cost O(number of steps) per call; for long simulations use coarse steps (or load aggregation,
not implemented).
"""
function borehole_wall_temperature(g::InfiniteLineSource, tq, q, t)
    s = zero(promote_type(eltype(q), typeof(g.Tg), typeof(t)))
    qprev = zero(eltype(q))
    for i in eachindex(tq, q)
        tq[i] < t || break
        s += (q[i] - qprev) * ILS(t - tq[i], g.rb, g.α)
        qprev = q[i]
    end
    return g.Tg - s / (2π * g.λ)
end
