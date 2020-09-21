"""
Uplift.py
Performs the uplift calculation, will use the parameters as defined in Input.py

@author: nvandervegt
"""

# Import
import numpy as np;
import scipy.stats as sc;
from Input import *;


# The Uplift class, including calculations
class Uplift:
    
    # Init, empty
    def __init__(self):
        pass;
        
    
    # Deterministic calculation
    def DeterministicCalculation(self):
        print("Starting uplift calculation");
        lambda_h = (k.charupper * D.charupper * d.charlower / kh.charlower) ** 0.5;
        lambda_ = (lambda_h / (Lf.charlower + B.charlower + lambda_h)) * np.exp((0.5 * B.charlower - xexit.charlower) / lambda_h);
        phi_exit = hp.charlower + lambda_*(h.charupper - hp.charlower);
        dphi = phi_exit - hp.charlower;
        dphi_cu = d.charlower * ((ysat.charlower - yw.charupper) / yw.charupper);
        
        # Calculating output
        self.FoS = dphi_cu/dphi;
        self.Pf = sc.norm.pdf((np.log(self.FoS/0.48) + 0.27 * self.BetaReq) / 0.46);
        self.beta = f.ProbabilityToReliabilityIndex(self.Pf);
        print("Factor of Safety:", self.FoS);
        print("Probability of Failure:", self.Pf);
        print("Reliability Index", self.beta);
            
    
    # Set reliability index for the levee segment
    def SetBetaReq(self, _betaReq):
        # Todo: create a better function name/structure
        self.BetaReq = _betaReq;
        
    
    # Monte Carlo analysis
    def MonteCarloCalculation(self):
        # Using hardcoded or OpenTurns MC?
        # Hardcoded -> put a limit to stop calculation (x number of samples, a target COV (N = 1/(Target_cov**2 * Pf)) <- First Pf estimation by FORM?)
        pass;