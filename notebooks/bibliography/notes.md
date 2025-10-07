
## GSHP System Performance and Control Challenges

______________________Ground source heat pump systems have gained significant attention as a key technology for building decarbonization, with the global market maintaining an annual growth rate of around 10% over the last decade, reaching 77.5 GW of worldwide installed capacity as of 2020 [@Beaudry2022]. However, despite their demonstrated efficiency advantages - showing up to 20% higher efficiency than air source heat pumps under equivalent conditions [@Menegazzo2022] - optimal operation remains a complex challenge, particularly regarding flow control strategies.

The fundamental challenge in GSHP operation lies in the complex interaction between multiple subsystems: the heat pump thermodynamic cycle, the ground heat exchanger, and the circulation pumps. Buildings consume more than 30% of global energy usage [@Ma2020], with heating, ventilation, and air-conditioning systems accounting for about 65% of this consumption in developed countries [@Li2024]. This highlights the critical importance of optimizing GSHP performance, where even small improvements in control strategies can yield substantial energy savings.

## Flow Control Strategies in GSHP Systems

### Constant vs. Variable Flow Operation

Traditional GSHP systems often operate at constant flow rates determined under rated full-load conditions for simplicity and reduced cost [@Li2024]. However, research consistently demonstrates the limitations of this approach. @Madani2011 showed that when on/off controlled GSHP systems are dimensioned to cover only 55% of peak heat demand, the seasonal performance factor (SPF) becomes significantly lower than variable capacity systems due to increased auxiliary heating requirements. Conversely, when systems cover more than 65% of peak demand, performance differences become negligible.

The potential for variable flow control has been extensively documented. @Edwards2015 developed a generalized control strategy for optimal part-load operation, demonstrating SPF improvements of 20-40% when using optimal flow control compared to nominal flow operation. Similarly, @Montagud2013 reported energy savings potential of up to 28% of annual electricity consumption through optimization of water circulation pump frequencies in their experimental GSHP installation.

### Pump Speed and Flow Rate Optimization

Multiple studies have focused specifically on circulation pump optimization. @CerveraVazquez2014 extended the in-situ optimization methodology to multistage heat pump systems, achieving energy savings up to 32% by optimizing water circulation pump frequencies. @DelCol2015 investigated variable speed drives across all system components (compressor, water pumps, and fans), identifying the ratio of compressor power to pump power as a key control variable.

The relationship between flow rates and system performance is complex and non-linear. @Edwards2015 found that optimal flow rates vary non-linearly with capacity, and that internal and external flow rates can be optimized independently while still achieving system optimization. @Zhou2024 demonstrated that variable flow operation in medium-depth systems could reduce total power consumption by 22.1% compared to constant flow operation.

## Advanced Control Strategies

### Model Predictive Control and Optimization

Recent advances have moved beyond simple rule-based control toward sophisticated optimization approaches. @DeRidder2011 developed dynamic programming techniques for borehole thermal energy storage systems, while @Hu2015 implemented extremum seeking control for hybrid GSHP systems, achieving nearly optimal efficiency without requiring detailed plant models.

@Gao2016 and @Gao2022 focused on small-scale GSHP systems, developing optimization strategies that minimize energy consumption while maintaining thermal comfort. Their work demonstrated energy savings of 9.59% in spring conditions and improved temperature control robustness.

### Load-Predictive and Adaptive Control

The importance of adaptive control that responds to varying thermal loads has been emphasized across multiple studies. @Li2024 investigated load-predictive control strategies, achieving 9.9-11.0% energy savings through modulation of flow rates and supply water temperatures. @Beaudry2022, in the case of SCW demonstrated that variable flow rate control maintaining a 2°C temperature difference across the heat exchanger resulted in 8-11% annual energy savings compared to constant flow operation.

### Traditional Simulation Environments

The GSHP modeling landscape has been dominated by several established tools, each with distinct advantages and limitations. **TRNSYS** has been widely adopted for dynamic system simulation, particularly for long-term performance analysis. @Capozza2015 used TRNSYS for multi-year analysis of GSHP systems, while @Bordignon2024 employed it for dual-source heat pump validation against real building data. @Li2024 developed transient GSHP models in TRNSYS for evaluating advanced control strategies.

**Engineering Equation Solver (EES)** has been favored for steady-state analysis and component-level modeling. @Corberan2010 developed quasi-steady state mathematical models using EES, demonstrating its effectiveness for system-level performance prediction and control optimization studies.

**Modelica** represents the object-oriented, equation-based modeling paradigm. @Hu2015 used Modelica Buildings Library for developing dynamic simulation models of hybrid GSHP systems, taking advantage of its acausal modeling capabilities and extensive component libraries.

### Emerging Computational Approaches

More recent work has explored advanced computational frameworks. @Salhein2022 implemented model predictive control using MATLAB's ODE solvers for ground heat exchanger optimization.

**Python-based approaches** have gained traction for their flexibility and open-source nature, though specific GSHP applications remain limited in the literature. Co-simulation approaches combining multiple tools have also been explored to leverage the strengths of different platforms [@Somogyi2016].

### Limitations of Current Tools

@Atam2015a and @Atam2015b provided comprehensive reviews highlighting significant limitations in current modeling approaches. Key bottlenecks include:

- **Computational efficiency**: Traditional tools often struggle with the coupled thermal-hydraulic nature of GSHP systems
- **Symbolic manipulation**: Limited capability for automatic differentiation and symbolic equation handling
- **Licensing costs**: Many commercial tools present barriers to research and development
- **Integration complexity**: Difficulty in coupling detailed component models with system-level optimization

### ModelingToolkit.jl: A Modern Alternative

This study employs **ModelingToolkit.jl**, a Julia-based framework that addresses many limitations of traditional tools. Key advantages include:

- **Open-source accessibility**: No licensing barriers for research or commercial use
- **Symbolic computation**: Full symbolic derivation and automatic differentiation capabilities
- **Modern solvers**: Access to state-of-the-art nonlinear and differential equation solvers
- **Performance**: Julia's just-in-time compilation provides computational efficiency comparable to compiled languages
- **Modularity**: Acausal, component-based modeling similar to Modelica but with enhanced flexibility

The choice of ModelingToolkit.jl enables the detailed coupled thermal-hydraulic analysis presented in this work while maintaining computational efficiency and model transparency.

## Research Gaps and Contributions

Despite extensive research on GSHP optimization, several gaps remain:

1. **Rule-based vs. Optimal Control**: While advanced optimization techniques show promise, most practical implementations still rely on simple rule-based strategies due to complexity and cost considerations [@Menegazzo2022].

2. **Turbulent Transition Effects**: Limited research has specifically addressed the impact of laminar-to-turbulent flow transitions in borehole heat exchangers on optimal control strategies.

3. **Integrated Thermal-Hydraulic Modeling**: Few studies provide detailed coupled analysis of thermal and hydraulic effects with sufficient detail to inform practical control strategies.

4. **Modern Computational Tools**: The GSHP community has been slow to adopt modern computational frameworks that could enable more sophisticated analysis.

Our study addresses these gaps by developing a detailed coupled thermal-hydraulic model using modern computational tools that specifically examines the relationship between flow control strategies and the fundamental physics of heat transfer and fluid flow in GSHP systems. Unlike previous rule-based approaches, our analysis provides insights into the physical mechanisms that determine optimal flow rates, particularly around the critical laminar-turbulent transition regions in borehole heat exchangers.
