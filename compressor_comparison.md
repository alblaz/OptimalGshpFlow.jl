# CompressorS vs DumbCompressor Comparison

This document explains why CompressorS works while DumbCompressor had issues.

## Key Factors That Make CompressorS Work

1. **Direct Use of Clapeyron.jl Functions**:
   - CompressorS uses direct calls to Clapeyron functions with clear phase specification
   - No intermediate custom functions that might introduce errors
   
   ```julia
   s_in ~ Clapeyron.entropy(modelref, P_in, T_in, phase="unknown")
   h_in ~ Clapeyron.enthalpy(modelref, P_in, T_in, phase="unknown")
   ```

2. **Well-Balanced Equations and Variables**:
   - The CompressorS model has exactly the right number of equations for its variables
   - No redundant or conflicting constraints

3. **Simple ModelingToolkit Setup**:
   ```julia
   @named C1 = CompressorS()
   C1p = C1 |> complete |> simplify
   prob_comp = NonlinearProblem(C1p, [], [])
   sol_comp = NonlinearSolve.solve(prob_comp)
   ```

4. **Clear Parameter and Variable Separation**:
   - Parameters are defined separately from variables
   - Initial guesses are provided for variables

## Issues Fixed in CompressorS

1. **Removed RefrigerantProps Component**:
   - The original version used `Properties.RefrigerantProps` which relied on EntropyScaling
   - EntropyScaling was causing domain errors during thermal conductivity calculations
   - We simplified to use direct Clapeyron calls only

2. **Removed Molecular Weight Scaling**:
   - Division by Mw was causing unit inconsistencies
   - Now using consistent units throughout

3. **Consistent Phase Specification**:
   - Using "phase='unknown'" consistently in all Clapeyron calls
   - This helps Clapeyron make the right decisions about which phase to use

## How to Fix DumbCompressor

Apply these same principles to DumbCompressor:

1. Use direct Clapeyron calls where appropriate
2. Ensure balanced equations (no redundant constraints)
3. Fix the solver success check: use `ModelingToolkit.ReturnCode.Success` 
   instead of the symbol `:Success`
4. Only constrain what's necessary in test functions

In Level 2 of DumbCompressor, you likely have an extra equation from redundant 
constraints like fixed frequency or redundant isentropic conditions.
