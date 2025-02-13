# Compressible Aerodynamics Code
## Overview

This repository contains a set of MATLAB scripts and functions designed to model compressible aerodynamics, including key calculations for shock waves, flow properties, and performance analysis in high-speed flows. The code handles both subsonic and supersonic regimes, offering a range of methods for solving the governing equations of compressible flow.
Features

### Working Features

+ Isentropic Relations: Solve for isentropic/stagnation relations 
+ Normal Shock Relations: Solve for downstream/upstream ratios and conditions 

### Future Capability

+  Positional and oblique shock solving: Compute shock positions, shock angles, and oblique shock relations  [TBD]
+  Static Thermodynamic Values: Compute static thermodynamic values across compressions and expansions [TBD]
+  Succcessive shock/expansion train problem
+  Diamond airfoil problem
+  Rocket nozzle problem

## Usage

### Isentropic Flow Relations
To calculate isentropic flow values, open compressible_flows_main.py, and edit the inputs of the following line to reflect the analysis. 

`config = isentropic_1D_ratios(1.4, {"Mach": 2.1})`

In the above case, 1.4 is equal to 'gamma', the specific heat ratio of the working fluid. 
The object input parameter may take any of the following values:
+ "Mach" - Mach number of the fluid
+ "Pt/P" - Stagnation pressure ratio
+ "Tt/T" - Stagnation temperature ratio
+ "Rhot/Rho" - Stagnation density ratio
+ "A/A*" - Sonic area ratio
+ "Mach_angle" - angle of pressure wave / mach wave propagation
+ "P-M Angle" - prandtl-meyer function angle

### Normal Shock Relations
To calculate normal shock ratios, open compressible_flows_main.py and edit the inputs of the following line: 

`config = normal_shock_ratios(1.4, {"P02/P01": 0.4986428})`

In the above case, 1.4 is equal to 'gamma', the specific heat ratio of the working fluid. 
The object input parameter may take any of the following values:
+ "Mach1" - Upstream mach number 
+ "Mach2" - Downstream mach number
+ "P2/P1" - Static pressure ratio
+ "T2/T1" - Static temperature ratio
+ "Rho2/Rho1" - Static density ratio
+ "P02/P01" - Stagnation pressure ratio
