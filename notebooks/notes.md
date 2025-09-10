# Research Notes - Organized

## 1. Hardware & Control Considerations
- Inverter losses comment, trend over the years
- Comment on the inverter
- Maybe justify not Rbstar

## 2. Modeling Paradigm & Philosophy
- Literature is also about the modelling paradigm
- More broadly, components based modelling, why?
- Physics based in our case
- Performance based?
- Focus on what it does vs. focus on how it's implemented
- Implementation details matter -> why did we make these choices

## 3. Software Tools & Platforms Comparison
### Existing Tools:
- EES (Engineering Equation Solver)
- Modelica (main contender)
- Trnsys? type 557 -> "DST"?
- Matlab?
- GSHP loop designer
- Python: pygfunction
- Hatef: EES and Trnsys

### Our Choice - Julia/ModelingToolkit:
- Doing in Julia/ModelingToolkit -> open-source
- Symbolic computation has been a powerful tool for modeling
- Open source in the end
- Flexible language in MTK
- Ready made models already exist
- Need to justify

## 4. Technical Implementation Advantages
- Technical reasons: equations based: symbolics/math. All solved at the same time
- Can also size components
- From the evaporator to the HTC, all is there
- Manage the levels of details
- Integrated in one platform
- Coupling to other tools and so on
- Connecting to other tools: Clapeyron
- ADable, gradient descent
- Surrogatize

## 5. Literature & Comparison Strategy
- Relevant to compare to other potential tools EES
- Complementary with next paper, flexibility

## 6. Incomplete Notes (TODO):
- Don't lose the 
- Comment that now it's the 