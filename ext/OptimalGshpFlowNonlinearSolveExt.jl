module OptimalGshpFlowNonlinearSolveExt

using OptimalGshpFlow: ComponentSystem, residual!
import NonlinearSolve: NonlinearProblem

"""
    NonlinearProblem(sys::ComponentSystem, u0; kwargs...)

Nonlinear problem `residual!(du, u, sys) = 0` from the initial guess `u0` (see
`initial_guess`); `kwargs` go to `NonlinearProblem`.
"""
NonlinearProblem(sys::ComponentSystem, u0; kwargs...) = NonlinearProblem(residual!, u0, sys; kwargs...)

end
