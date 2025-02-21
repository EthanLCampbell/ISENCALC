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

    def shock_angle(self, M, delta, gamma=1.4):
        '''
        function to evaluate shock angle numerically
        '''
        delta_rad = np.radians(delta)
    
        def equation(beta_rad):
            tan_delta = 2 * (1 / np.tan(beta_rad)) * ((M**2 * np.sin(beta_rad)**2 - 1) / (M**2 * (gamma + np.cos(2 * beta_rad)) + 2))
            return tan_delta - np.tan(delta_rad)
        
        beta_guess = np.radians(delta + 10)  # Initial guess
        beta_solution = fsolve(equation, beta_guess)
    
        return np.degrees(beta_solution[0])
        
    def compute_ratios(self):
        """Compute oblique shock relations based on given gamma and one known parameter."""
        results = {}
        gam = self.gamma
        gm1 = gam - 1
        gp1 = gam + 1
        delta = 0
        # check turning angle included
        if "delta" in self.param:
            delta = self.param["delta"]
        else:
            print("No turning angle (delta) input")
            return
        
        if "Mach1" in self.param:        
            M1 = self.param["Mach1"]
            beta = self.shock_angle(M1,delta,1.4)
            beta_r = np.radians(beta)
            M1N = M1*np.sin(beta_r)
            shock_config = normal_shock_ratios(1.4, {"Mach1": M1N})
            M2N = shock_config.results["Mach2"]
            M2 = M2N/np.sin(beta_r - np.radians(delta))
            results["Mach1"] = M1
            results["Mach1n"] = M1N
            results["Mach2n"] = M2N
            results["Mach2"] = M2
            results["delta"] = delta
            results["beta"] = beta
            results["P2/P1"] = (2 * gam * M1**2 * np.sin(beta_r)**2 - gm1) / gp1
            results["T2/T1"] = ((2 * gam * M1**2 * np.sin(beta_r)**2 - gm1) * (gm1 * M1**2 * np.sin(beta_r)**2 + 2)) / (gp1**2 * M1**2 * np.sin(beta_r)**2)
            results["Rho2/Rho1"] = (gp1 * M1**2 * np.sin(beta_r)**2) / (gm1 * M1**2 * np.sin(beta_r)**2 + 2)
            results["P02/P01"] = (gp1 * M1**2 * np.sin(beta_r)**2 / (gm1 * M1**2 * np.sin(beta_r)**2 + 2))**(gam / gm1) * (gp1 / (2 * gam * M1**2 * np.sin(beta_r)**2 - gm1))**(1 / gm1)
        
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
#config = oblique_shock_ratios(1.4, {"Mach1": 2.5, "delta": 20})
#print(config)




