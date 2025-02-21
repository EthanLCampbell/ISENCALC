#------------------------------------------------------------------------------+
#
# 1D Normal Shockwave Relations
# Computes all shock relations as functions of M, gam
# 
#------------------------------------------------------------------------------+

import math
import numpy as np
from scipy.optimize import fsolve

class normal_shock_ratios:
    def __init__(self, gamma: float, param: dict):
        """
        Initialize the NormalShockRelations object.
        
        :param gamma: Specific heat ratio (gamma)
        :param param: Dictionary of known parameters (e.g., Mach number, pressure ratio, etc.)
        """
        self.gamma = gamma
        self.param = param
        self.results = self.compute_ratios()
    
    def compute_ratios(self):
        """Compute normal shock relations based on given gamma and one known parameter."""
        results = {}
        gam = self.gamma
        gm1 = gam - 1
        gp1 = gam + 1

        if "Mach1" in self.param:  # Compute based on upstream Mach number
            M1 = self.param["Mach1"]
            results["Mach1"] = M1
            phi = 1+gm1/2*M1**2
            results["Mach2"] = math.sqrt((1 + (gm1 / 2) * M1**2) / (gam * M1**2 - gm1 / 2))
            results["P2/P1"] = 1 + (2 * gam / (gam + 1)) * (M1**2 - 1)
            results["T2/T1"] = ((1 + (gm1 / 2) * M1**2) * ((2 * gam / (gam + 1) * M1**2 - gm1 / (gam + 1))) / M1**2)
            results["Rho2/Rho1"] = (gp1 * M1**2) / (gm1 * M1**2 + 2)
            results["P02/P01"] = (gp1/2*M1**2 / (phi))**(gam/gm1) * (1/(2*gam/gp1*M1**2 - gm1/gp1))**(1/gm1)

        if "Mach2" in self.param: # Compute based on downstream mach number
            M2 = self.param["Mach2"]
            def get_M1_from_M2(M1):
                """Equation to solve: M2 - computed M2 from M1 = 0."""
                return np.sqrt((1 + gm1 / 2 * M1**2) / (gam * M1**2 - gm1 / 2)) - M2
            M1 = fsolve(get_M1_from_M2, 2.0)
            phi = 1+gm1/2*M1**2
            results["Mach1"] = M1
            results["P2/P1"] = 1 + (2 * gam / (gam + 1)) * (M1**2 - 1)
            results["T2/T1"] = ((1 + (gm1 / 2) * M1**2) * ((2 * gam / (gam + 1) * M1**2 - gm1 / (gam + 1))) / M1**2)
            results["Rho2/Rho1"] = (gp1 * M1**2) / (gm1 * M1**2 + 2)
            results["P02/P01"] = (gp1/2*M1**2 / (phi))**(gam/gm1) * (1/(2*gam/gp1*M1**2 - gm1/gp1))**(1/gm1)

        if "P2/P1" in self.param:  # Compute based on static pressure ratio
            P2P1 = self.param["P2/P1"]
            def find_mach(M):
                return 1 + (2 * gam / gp1) * (M**2 - 1) - P2P1
            M1 = fsolve(find_mach, 2.0)[0]  # Initial guess for supersonic M1
            phi = 1+gm1/2*M1**2
            results["Mach1"] = M1
            results["Mach2"] = math.sqrt((1 + (gm1 / 2) * M1**2) / (gam * M1**2 - gm1 / 2))
            results["P2/P1"] = P2P1
            results["T2/T1"] = ((1 + (gm1 / 2) * M1**2) * ((2 * gam / (gam + 1) * M1**2 - gm1 / (gam + 1))) / M1**2)
            results["Rho2/Rho1"] = (gp1 * M1**2) / (gm1 * M1**2 + 2)
            results["P02/P01"] = (gp1/2*M1**2 / (phi))**(gam/gm1) * (1/(2*gam/gp1*M1**2 - gm1/gp1))**(1/gm1)
            
        if "T2/T1" in self.param: # Compte based on static temperature ratio
            T2oT1 = self.param["T2/T1"]
            def get_M1_from_T2oT1(M1):
                return (1 + (gm1 / 2) * M1**2) * (2 * gam / gp1 * M1**2 - gm1 / gp1) / M1**2 - T2oT1
            M1 = fsolve(get_M1_from_T2oT1,2.0)
            phi = 1+gm1/2*M1**2
            results["Mach1"] = M1
            results["P2/P1"] = 1 + (2 * gam / (gam + 1)) * (M1**2 - 1)
            results["T2/T1"] = T2oT1
            results["Rho2/Rho1"] = (gp1 * M1**2) / (gm1 * M1**2 + 2)
            results["P02/P01"] = (gp1/2*M1**2 / (phi))**(gam/gm1) * (1/(2*gam/gp1*M1**2 - gm1/gp1))**(1/gm1)

        if "Rho2/Rho1" in self.param:
            results["Bro"] = "Why would you literally ever calculate mach number from denisty ratio??"
            results["Funny_joke"] = "Why do norwegian ships have barcodes on them?\n \t\tSo they can scandanavian"

        if "P02/P01" in self.param:
            P02oP01 = self.param["P02/P01"]
            def get_M1_from_P02oP01(M1):
                return (gp1/2*M1**2 / (1+gm1/2*M1**2))**(gam/gm1) * (1/(2*gam/gp1*M1**2 - gm1/gp1))**(1/gm1) - P02oP01
            M1 = fsolve(get_M1_from_P02oP01,2.0)
            phi = 1+gm1/2*M1**2
            results["Mach1"] = M1
            results["Mach2"] = math.sqrt((1 + (gm1 / 2) * M1**2) / (gam * M1**2 - gm1 / 2))
            results["P2/P1"] = 1 + (2 * gam / (gam + 1)) * (M1**2 - 1)
            results["T2/T1"] = ((1 + (gm1 / 2) * M1**2) * ((2 * gam / (gam + 1) * M1**2 - gm1 / (gam + 1))) / M1**2)
            results["Rho2/Rho1"] = (gp1 * M1**2) / (gm1 * M1**2 + 2)
            results["P02/P01"] = P02oP01

    
        return results

    def __repr__(self):
        result_str = f"NormalShockRelations(\n  gamma={self.gamma},\n  param={self.param},\n  results=\n"
        for key, value in self.results.items():
            result_str += f"    {key}: {value}\n"
        result_str += ")"
        return result_str

'''
Example inputs (all equivalent to M1=2.5)
    "Mach1": 2.5
    "Mach2": 0.512989176042577
    "P2/P1": 7.125
    "T2/T1": 2.565
    "Rho2/Rho1": 3.3333333333333335
    "P02/P01": 0.4986428
'''
## RUN HERE:
#config = normal_shock_ratios(1.4, {"T2/T1": 2.565})
#config = normal_shock_ratios(1.4, {"P02/P01": 0.4986428})
#print(config)




