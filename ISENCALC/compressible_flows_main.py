#------------------------------------------------------------------------------+
#
#                           ==== ISENCALC ====
#                         Ethan Labianca-Campbell
#               Purdue School of Aeronautics & Astronautics
#               Aerodynamics and compressible flow functions 
#                        Last Update: Feb 13, 2025
# 
#------------------------------------------------------------------------------+

#-------Libraries and Dependencies:--------------------------------------------+
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import os
from isentropic_relations import isentropic_1D_ratios
from normal_shock_relations import normal_shock_ratios
from oblique_shock_relations import oblique_shock_ratios
#-------ISENTROPIC PROBLEMS----------------------------------------------------+

'''
Example inputs (all equivalent to M = 2.5)
    "Mach": 2.5
    "Pt/P": 17.0859385118
    "Tt/T": 2.25
    "Rhot/Rho": 7.59375016137
    "A/A*":  2.63671875
    "Mach_angle":  23.5781784
    "P-M Angle": 39.12356401
'''

# RUN ISENTROPIC:
print("\n\n=== 1D ISENTROPIC RELATIONS OUTPUTS === \n")
config = isentropic_1D_ratios(1.4, {"Mach": 2.1})
print(config)

#-------SHOCK PROBLEMS---------------------------------------------------------+

# Normal Shockwave
'''
Example inputs (all equivalent to M1=2.5)
    "Mach1": 2.5
    "Mach2": 0.512989176042577
    "P2/P1": 7.125
    "T2/T1": 2.565
    "Rho2/Rho1": 3.3333333333333335
    "P02/P01": 0.4986428
'''
print("\n\n=== NORMAL SHOCK RELATIONS OUTPUTS === \n")
config = normal_shock_ratios(1.4, {"P02/P01": 0.4986428})
print(config)

# Oblique Shockwave
print("\n\n=== OBLIQUE SHOCK RELATIONS OUTPUTS === \n")
config = oblique_shock_ratios(1.4, {"Mach1": 2.5, "delta": math.radians(20)})
print(config)

