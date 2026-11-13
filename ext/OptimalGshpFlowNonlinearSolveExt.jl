module OptimalGshpFlowNonlinearSolveExt

using OptimalGshpFlow: ComponentSystem, residual!
import OptimalGshpFlow: continuation
import NonlinearSolve: NonlinearProblem, solve
using NonlinearSolve: SciMLBase

"""
    NonlinearProblem(sys::ComponentSystem, u0; kwargs...)

Nonlinear problem `residual!(du, u, sys) = 0` from the initial guess `u0` (see
`initial_guess`); `kwargs` go to `NonlinearProblem`.
"""
NonlinearProblem(sys::ComponentSystem, u0; kwargs...) = NonlinearProblem(residual!, u0, sys; kwargs...)

function _solve(sys, u0, alg; kwargs...)
    sol = try
        solve(NonlinearProblem(sys, u0), alg; kwargs...)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    sol === nothing && return (ok = false, u = u0, iterations = -1)
    ok = SciMLBase.successful_retcode(sol) && all(isfinite, sol.u)
    return (ok = ok, u = sol.u, iterations = sol.stats.nsteps)
end

function continuation(build, λs, u0, alg; secant::Bool = true, max_halvings::Integer = 6, kwargs...)
    λs = collect(λs)
    results = NamedTuple[]
    first_step = _solve(build(λs[1]), u0, alg; kwargs...)
    push!(results, (λ = λs[1], first_step...))
    first_step.ok || return results
    λprev, uprev = λs[1], first_step.u
    λold, uold = λprev, nothing                   # point before λprev, for the secant predictor
    for λ in λs[2:end]
        target, its, step_ok = λ, 0, false
        h, halvings = λ - λprev, 0
        while true
            λn = λprev + h
            abs(λn - target) <= 1e-12 * max(1, abs(target)) && (λn = target)
            guess = (secant && uold !== nothing && λprev != λold) ?
                    uprev .+ (uprev .- uold) .* ((λn - λprev) / (λprev - λold)) : copy(uprev)
            r = _solve(build(λn), guess, alg; kwargs...)
            # a failed secant guess is retried from the previous solution
            if !r.ok && secant && uold !== nothing
                r = _solve(build(λn), copy(uprev), alg; kwargs...)
            end
            its += max(r.iterations, 0)
            if r.ok
                λold, uold, λprev, uprev = λprev, uprev, λn, r.u
                λn == target && (step_ok = true; break)
                h = target - λprev                 # try the rest of the step
            else
                halvings += 1
                halvings > max_halvings && break
                h /= 2
            end
        end
        push!(results, (λ = λ, ok = step_ok, u = uprev, iterations = its))
        step_ok || break
    end
    return results
end

end
