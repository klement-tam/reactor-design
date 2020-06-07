# reactor-design 

Design and optimisation of a circulating fluidized bed reactor (CFB) to produce maleic anhydride (MA) through catalytic oxidation of n-butane under steady state operation. 

## Description of project
This repository contains a all-in-one, extensive MATLAB file which designs a catalytic Fludisied Bed Reactor to maximise the MA yield in a CFB reactor, using given inlet conditions and reactor geometry. The key decision variable to be optimised is the catalyst recirculation rate. 

The performance of the reactor is simulated and plots are created for two extreme cases:\
    1. Adiabatic system: Heat exchange with the surroundings equals zero. \
    2. Isothermal system: All the heat generated during the reaction is removed.

For each system, the overall performance were shown through the plots: \ 
    - Flow of maleic anhydride (MA) vs. catalyst recirculation rate (function to be maximized)\
    - Outlet temperature vs. catalyst recirculation rate \
    - n-Butane conversion vs. catalyst recirculation rate \
    - Maleic anhydride selectivity vs. catalyst recirculation rate \
    - Maleic anhydride selectivity vs. n-butane conversion 

Once the optimium catalyst recirculation rate was found, the profiles for this recirculation rate were also plotted to gain a depeer analysis of how the reactor works: \
    - Flow of maleic anhydride, n-butane, CO2, H2O, CO vs. riser height \
    - Temperature vs. riser height \
    - Pressure vs. riser height\
    - 1-ɛ vs. riser height \
    - CO Selectivity vs. extent of reaction 1\
    - CO/CO2 vs. extent of reaction 1 

Sensitivity analysis was also performed performed on various aspects of the design: \
    - Errors in certain kinetic reaction parameters (b1, b2, n)\
    - n-butane:air feed ratio\
    - Inlet temperature\
    - Inlet pressure 

Lastly, multivariable optimisations were also performed to further increase yields.

A series of Ordinary Differential Equations (ODEs) had to be developed for mass & energy balances pre-code. The ODE solver tool ode23 is used for the purpose of this reactor design.

## Usage 
The way this code is written is that instead of plotting everything at once (i.e. Adiabatic, Isothermal, sensitivity analysis for all kinetic parameters etc.), it uses dialogue prompts to allow the user to first select the condition of the system, and second which parameter is desired for sensitivity analysis, thereafter plotting the relevant graphs. This gives the code more flexibility and efficiency. 

## Plots
All plots can be found in results folder.
