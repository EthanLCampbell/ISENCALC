#------------------------------------------------------------------------------+
#
#                           ==== ISENCALC ====
#                         Ethan Labianca-Campbell
#               Purdue School of Aeronautics & Astronautics
#               Aerodynamics and compressible flow functions 
#                       Last Update: Feb 12, 2025
# 
#------------------------------------------------------------------------------+

#-------Libraries and Dependencies:--------------------------------------------+
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import os
from isentropic_relations import isentropic_1D_ratios

#-------ISENTROPIC PROBLEMS----------------------------------------------------+

# Inputs for isentropic relations (all examples equivalent to M=2.5)
'''
"Mach": 2.5
"Pt/P": 17.0859385118
"Tt/T": 2.25
"Rhot/Rho": 7.59375016137
"A/A*":  2.63671875
"Mach_angle":  23.5781784
"P-M Angle": 39.12356401
'''

# RUN ISENTROPIC:
config = isentropic_1D_ratios(1.4, {"Mach": 2.1})
print(config)

#-------SHOCK PROBLEMS---------------------------------------------------------+

# [TBD]
