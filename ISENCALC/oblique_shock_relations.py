#------------------------------------------------------------------------------+
#
# Oblique Shock Relations
# Computes all shock relations as functions of M, gam, delta (turning angle)
# 
#------------------------------------------------------------------------------+

import math
import numpy as np
from normal_shock_relations import normal_shock_ratios
from scipy.optimize import fsolve

print("\n !!Oblique Shock Module WIP!! \n")

class oblique_shock_ratios:
    def __init__(self, gamma: float, param: dict):
        """
        Initialize the ObliqueShockRelations object.
        
        :param gamma: Specific heat ratio (gamma)
        :param param: Dictionary of known parameters (e.g., Mach number, shock angle, etc.)
        """
        self.gamma = gamma
        self.param = param
        self.results = self.compute_ratios()

    def compute_ratios(self):
        """Compute oblique shock relations based on given gamma and one known parameter."""
        results = {}
        gam = self.gamma
        gm1 = gam - 1
        gp1 = gam + 1

        # input calculation logic here!!!! 
        

        return results

    def __repr__(self):
        result_str = f"ObliqueShockRelations(\n  gamma={self.gamma},\n  param={self.param},\n  results=\n"
        for key, value in self.results.items():
            result_str += f"    {key}: {value}\n"
        result_str += ")"
        return result_str




## RUN HERE:

# delta required for turning angle
# other input can be Mach1....

config = oblique_shock_ratios(1.4, {"Mach1": 2.5, "delta": math.radians(20)})
print(config)






    
