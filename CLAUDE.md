# OptimalGshpFlow.jl: notes for Claude Code

First-principles, steady-state simulation of a ground-source heat pump: a moving-boundary
evaporator (two-phase and superheat zones, ε-NTU) coupled to a borehole (Hellström model),
solved as a nonlinear system with NonlinearSolve. Bounded unknowns are mapped to ℝ (sigmoid,
including on the pressure) so that unconstrained solvers can be used.

**Status, roadmap, decisions and the plan for the next branch:** `docs/dev/ROADMAP.md`.
Read it at the start of a session and update it when a branch is finished.

## Workflow
- One feature per branch, branched from up-to-date `main`. Never commit to `main`. Concerns
  that are independent in theory (properties, variable transforms, component framework,
  glide, moving boundary, MTK, dynamics) go in separate branches and PRs.
- Bugs: report, add a test that fails on the current code, commit it, then fix it in a
  second commit.
- Open PRs with `/opt/homebrew/bin/gh` (it is not on the PATH of the Bash tool). PR bodies
  end with the Claude Code attribution line.
- Never discard uncommitted work: no `git reset --hard` or `git checkout -- .` with a dirty
  tree; use `git stash` or a separate worktree.
- `assets/docs/` holds reference papers (gitignored): Qiao 2025 on ε-NTU with glide, and an
  article on the laminar-turbulent transition in horizontal pipes (not used yet).
- `legacy` branch: the original code (main at 5c0a87d) plus the sigmoid fix, kept for reference.

## Julia
- Run Julia only through the julia-mcp tool, never `julia` on the command line.
  - Package env: `env_path` = repo root. Test-only packages (`[extras]`) cannot be loaded
    there; use a scratch env (e.g. `$CLAUDE_JOB_DIR/tmp/scratchenv` with the repo `dev`-ed).
  - Revise is not loaded: after editing `src/`, restart the session (`julia_restart`) before
    re-running code.
  - Results sometimes come back one call late; if an output looks like the previous call's,
    send a cheap follow-up call (`println("ping")`) to get the real result.
  - Run tests with `Pkg.test(; io = devnull)` inside the package session.
- Units are SI throughout; refrigerant properties are per kg.
- Keep model and property code generic in the number type, so ForwardDiff Duals pass through.
- Type stability and compile time: keep structs concretely typed (parametric fields, no
  untyped or abstract fields) and refrigerants/components concrete. Inference on an abstract
  type blows up inside Clapeyron (PR #32: untyped fields turned the study's first call into a
  20-minute compile hang). For type instability, compile hangs or long first calls, **start
  with the dedicated tools, before any hand bisecting or sampling**:
  - JET.jl `@report_opt target_modules = (OptimalGshpFlow,)` on the concrete entry points:
    lists runtime dispatch with file:line; `@test_opt` to lock a fix in a test.
  - SnoopCompile.jl (`@snoop_inference`, inference triggers): where compile time goes.
  - Cthulhu.jl `@descend`: interactive, cannot be driven through julia-mcp; suggest the user
    runs it in their REPL.
  - Install them in a scratch env, not the package. Also check lowered code for `Core.Box`
    (captured variables reassigned in closures, e.g. inside `@testset` or keyword functions).

## Code map
- `src/models.jl`: residual functions `evaporator_model!`, `evaporatorborehole_model!`,
  `borehole_model!`, plus constructors of the solution types.
- `src/evaporator_formulations.jl`: `PropertyFormulation` (`NestedPressure` default,
  `SaturationTemperature`, `ExplicitTv`), passed as last argument of the evaporator residuals.
- `src/fluidproperties.jl`: spline-based `Refrigerant` and `SecondaryFluid` (CoolProp sampling).
- `src/properties/interface.jl`: `AbstractRefrigerant` and the property functions the models use.
- `src/properties/clapeyron_refrigerant.jl`: `ClapeyronRefrigerant` (Clapeyron + CoolProp or
  EntropyScaling transport); pass `phase = :vapor/:liquid` whenever the phase is known.
  CoolProp goes through a cached `AbstractState`, never `PropsSI` in hot paths (≈ 100 µs/call).
- `src/properties/tv_state.jl`: explicit (T, v) property evaluations for `ClapeyronRefrigerant`;
  `saturation_T` (explicit p, v_l, v_v from the superancillaries of the CoolProp fluid files).
- `src/transforms.jl`: bounded unknowns as types (`Free`, `Unit`, `Box`, `LogBox`, `Above`;
  `from_R`, `to_R`); `pressure_transform`, `temperature_transform`, `volume_transform`,
  `evaporator_transforms` give the maps the models use.
- `src/components/`: component layer (`ComponentSystem`, `Evaporator`, `BoreholeHX`); ports joined
  by shared node unknowns, scaled residuals; `validity_issues` (brine freezing, pinch, bounds).
  `ext/`: `NonlinearProblem(sys, u0)` and `continuation` (NonlinearSolve).
- `src/utilities.jl`: low-level maps between bounded intervals and ℝ (linear, tanh, sigmoid).
- `src/borehole.jl`, `src/pipeflow.jl`, `src/plateheatexchanger.jl`, `src/compressor.jl`,
  `src/pressuredrop.jl`: component correlations.
- `examples/`: scripts used for the original study (need NonlinearSolve and Makie in the env).
