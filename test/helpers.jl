# Shared helpers for the tests.
#
# ForwardDiff with one tag: `ForwardDiff.derivative/gradient/jacobian` tag the dual numbers with the
# type of the differentiated function, so every closure gets its own Dual type and the model and
# property code below it (Clapeyron included) is compiled again for each one: ≈ 1 s for a single
# property, much more for a model Jacobian. With one tag and one chunk size the Dual types are the
# same in all tests, and that code is compiled once.
using ForwardDiff

struct TestTag end
const TT{V} = ForwardDiff.Tag{TestTag, V}
const JAC_CHUNK = ForwardDiff.Chunk{7}()      # model unknowns: 7 to 11

"""Derivative of the scalar function `f` at `x`."""
dderivative(f, x::Real) = ForwardDiff.extract_derivative(TT{typeof(x)}, f(ForwardDiff.Dual{TT{typeof(x)}}(x, one(x))))

"""Gradient of the scalar function `f` at `x`."""
dgradient(f, x::AbstractVector) =
    ForwardDiff.gradient(f, x, ForwardDiff.GradientConfig(f, x, ForwardDiff.Chunk{length(x)}(), TT{eltype(x)}()), Val(false))

"""Jacobian of the in-place function `f!(y, x)` at `x` (`chunk`: at most `length(x)`; keep the
default where possible, another chunk size is another Dual type)."""
djacobian(f!, y::AbstractVector, x::AbstractVector, chunk = JAC_CHUNK) =
    ForwardDiff.jacobian(f!, y, x, ForwardDiff.JacobianConfig(f!, y, x, chunk, TT{eltype(x)}()), Val(false))

# Plain Newton with the shared-tag Jacobian: solving one model in several test files compiles it once
# (a NonlinearProblem with a new closure compiles the whole NonlinearSolve stack again).
function newton_solve(f!, u; maxit = 30, tol = 1e-8, chunk = JAC_CHUNK)
    du = similar(u)
    for it in 1:maxit
        J = djacobian(f!, du, u, chunk)
        f!(du, u)
        maximum(abs, du) < tol && return u, it
        u = u - J \ du
    end
    return u, -1
end

# Spline refrigerant (CoolProp sampling, ≈ 3 s per construction): built once for all test files.
# Read-only; fluidproperties.jl tests the constructor itself.
const SPL_R290 = Refrigerant("R290")
