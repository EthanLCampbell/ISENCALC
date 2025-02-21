#------------------------------------------------------------------------------+
#
# 1D Isentropic Ratios
# Computes all isentropic ratios as functions of M, gam
# 
#------------------------------------------------------------------------------+

import math
from scipy.optimize import fsolve

class isentropic_1D_ratios:
    def __init__(self, gamma: float, param: dict):
        """
        Initialize the SimulationConfig object.
        
        :param gamma: Specific heat ratio (gamma)
        :param run_type: Type of simulation run
        :param param: Dictionary of additional parameters
        """
        self.gamma = gamma
        #self.run_type = run_type
        self.param = param
        self.results = self.compute_ratios()
    
    def compute_ratios(self):
        """Compute aerodynamic ratios based on given gamma and one known parameter."""
        results = {}
        gam = self.gamma
        gm1 = gam-1
        gp1 = gam+1
        if "Mach" in self.param: # compute based on mach number
            M = self.param["Mach"]
            phi = (1 + gm1 / 2 * M**2) # flow function
            results["Pt/P"] = phi ** (gam / gm1)
            results["Tt/T"] = phi
            results["Rho_t/Rho"] = (1 + gm1 / 2 * M**2) ** (1 / gm1)
            results["A/A*"] = (gp1/2)**(-gp1/(2*gm1)) * (phi)**(gp1/(2*gm1)) / M
            results["Mach_Angle"] = math.asin(1/M) * 180 / math.pi
            results["P-M Angle"] = (math.sqrt(gp1/gm1)*math.atan(math.sqrt(gm1/gp1*(M**2-1)))-math.atan(math.sqrt(M**2-1))) * 180/math.pi
        if "Pt/P" in self.param: # compute based on pressure ratio
            PtoP = self.param["Pt/P"]
            M = math.sqrt(2/gm1 * ((PtoP)**(gm1/gam)-1))
            results["M"] = M
            phi = (1 + gm1 / 2 * M**2)
            results["Pt/P"] = PtoP
            results["Tt/T"] = phi
            results["Rho_t/Rho"] = (1 + gm1 / 2 * M**2) ** (1 / gm1)
            results["A/A*"] = (gp1/2)**(-gp1/(2*gm1)) * (phi)**(gp1/(2*gm1)) / M
            results["Mach_Angle"] = math.asin(1/M) * 180 / math.pi
            results["P-M Angle"] = (math.sqrt(gp1/gm1)*math.atan(math.sqrt(gm1/gp1*(M**2-1)))-math.atan(math.sqrt(M**2-1))) * 180/math.pi
        if "Tt/T" in self.param: # compute based on temperature ratio
            TtoT = self.param["Tt/T"]
            M = math.sqrt((2 / (gam - 1)) * (TtoT - 1))
            results["M"] = M
            phi = (1 + gm1 / 2 * M**2)
            results["Pt/P"] = phi ** (gam/gm1)
            results["Tt/T"] = phi
            results["Rho_t/Rho"] = (1 + gm1 / 2 * M**2) ** (1 / gm1)
            results["A/A*"] = (gp1/2)**(-gp1/(2*gm1)) * (phi)**(gp1/(2*gm1)) / M
            results["Mach_Angle"] = math.asin(1/M) * 180 / math.pi
            results["P-M Angle"] = (math.sqrt(gp1/gm1)*math.atan(math.sqrt(gm1/gp1*(M**2-1)))-math.atan(math.sqrt(M**2-1))) * 180/math.pi
        if "Rhot/Rho" in self.param: # compute based on density ratio
            RhotoRho = self.param["Rhot/Rho"]
            M = math.sqrt((2 / gm1) * (RhotoRho ** (gm1) - 1))
            results["M"] = M
            phi = (1 + gm1 / 2 * M**2)
            results["Pt/P"] = phi ** (gam/gm1)
            results["Tt/T"] = phi
            results["Rho_t/Rho"] = (1 + gm1 / 2 * M**2) ** (1 / gm1)
            results["A/A*"] = (gp1/2)**(-gp1/(2*gm1)) * (phi)**(gp1/(2*gm1)) / M
            results["Mach_Angle"] = math.asin(1/M) * 180 / math.pi
            results["P-M Angle"] = (math.sqrt(gp1/gm1)*math.atan(math.sqrt(gm1/gp1*(M**2-1)))-math.atan(math.sqrt(M**2-1))) * 180/math.pi
        if "A/A*" in self.param: # we must use a root-finding method, given solving for M is an implicit equation
            AoAstar = self.param["A/A*"]
            def area_ratio_equation(M):
                return (1 / M) * ((2 / gp1) * (1 + gm1 / 2 * M**2))**(gp1 / (2 * gm1)) - AoAstar
            # supersonic case
            M_guess = 2.0
            M_sup = fsolve(area_ratio_equation,M_guess)
            results["M_sup"] = M_sup
            phi = (1 + gm1 / 2 * M_sup**2)
            results["Pt/P_sup"] = phi ** (gam/gm1)
            results["Tt/T_sup"] = phi
            results["Rho_t/Rho_sup"] = phi ** (1 / gm1)
            results["A/A*_sup"] = (gp1/2)**(-gp1/(2*gm1)) * (phi)**(gp1/(2*gm1)) / M_sup
            results["Mach_Angle_sup"] = math.asin(1/M_sup) * 180 / math.pi
            results["P-M Angle_sup"] = (math.sqrt(gp1/gm1)*math.atan(math.sqrt(gm1/gp1*(M_sup**2-1)))-math.atan(math.sqrt(M_sup**2-1))) * 180/math.pi

            # subsonic case
            M_guess = 0.3
            M_sub = fsolve(area_ratio_equation,M_guess)
            results["M_sub"] = M_sub
            phi = (1 + gm1 / 2 * M_sub**2)
            results["Pt/P_sub"] = phi ** (gam/gm1)
            results["Tt/T_sub"] = phi
            results["Rho_t/Rho_sub"] = phi ** (1 / gm1)
            results["A/A*_sub"] = (gp1/2)**(-gp1/(2*gm1)) * (phi)**(gp1/(2*gm1)) / M_sub
            # Mach waves and PM angles are ill-defined for subsonic speeds
            #results["Mach_Angle"] = math.asin(1/M_sub) * 180 / math.pi 
            #results["P-M Angle"] = (math.sqrt(gp1/gm1)*math.atan(math.sqrt(gm1/gp1*(M_sub**2-1)))-math.atan(math.sqrt(M_sub**2-1))) * 180/math.pi
        if "Mach_angle" in self.param: # compute based on mach angle
            Mach_angle = self.param["Mach_angle"]
            M = 1/math.sin(math.pi/180*Mach_angle)
            results["M"] = M
            phi = (1 + gm1 / 2 * M**2)
            results["Pt/P"] = phi ** (gam/gm1)
            results["Tt/T"] = phi
            results["Rho_t/Rho"] = (1 + gm1 / 2 * M**2) ** (1 / gm1)
            results["A/A*"] = (gp1/2)**(-gp1/(2*gm1)) * (phi)**(gp1/(2*gm1)) / M
            results["Mach_Angle"] = math.asin(1/M) * 180 / math.pi
            results["P-M Angle"] = (math.sqrt(gp1/gm1)*math.atan(math.sqrt(gm1/gp1*(M**2-1)))-math.atan(math.sqrt(M**2-1))) * 180/math.pi
        if "P-M Angle" in self.param:
            nu_radians = self.param["P-M Angle"]
            nu = math.radians(nu_radians)
            def prandtl_meyer_fan(M):
                return (math.sqrt(gp1 / gm1) * math.atan(math.sqrt(gm1 * (M**2 - 1) / gp1)) - math.atan(math.sqrt(M**2 - 1))) - nu
            M = fsolve(prandtl_meyer_fan,2.0)
            results["M"] = M
            phi = (1 + gm1 / 2 * M**2)
            results["Pt/P"] = phi ** (gam/gm1)
            results["Tt/T"] = phi
            results["Rho_t/Rho"] = (1 + gm1 / 2 * M**2) ** (1 / gm1)
            results["A/A*"] = (gp1/2)**(-gp1/(2*gm1)) * (phi)**(gp1/(2*gm1)) / M
            results["Mach_Angle"] = math.asin(1/M) * 180 / math.pi
            results["P-M Angle"] = (math.sqrt(gp1/gm1)*math.atan(math.sqrt(gm1/gp1*(M**2-1)))-math.atan(math.sqrt(M**2-1))) * 180/math.pi
            
        return results
    
    def __repr__(self):
        result_str = f"SimulationConfig(\n  gamma={self.gamma},\n  param={self.param},\n  results=\n"
        for key, value in self.results.items():
            result_str += f"    {key}: {value}\n"
        result_str += ")"
        return result_str


# various input parameters (all equivalent to M = 2.5)
# "Mach": 2.5
# "Pt/P": 17.0859385118
# "Tt/T": 2.25
# "Rhot/Rho": 7.59375016137
# "A/A*":  2.63671875
# "Mach_angle":  23.5781784
# "P-M Angle": 39.12356401

## RUN HERE:
#config = isentropic_1D_ratios(1.4, {"Mach": 2.1})
#print(config)
