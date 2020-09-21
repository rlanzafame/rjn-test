"""
Uplift.py
Performs the uplift calculation
"""

# Import
import numpy as np
from Input import *


# The Uplift class, including calculations
class Uplift:
    
    # Allow for access to the parameters within the class
    global k; global kh; global D; global d; global Lf; global B; global xexit; global h; global hp; global ysat; global yw; global mu;
   
    
    # Init, empty
    def __init__(self):
        pass
        
    
    # Deterministic calculation
    def DeterministicCalculation(self):
        print("Starting uplift calculation");
        lambda_h = (k.charupper * D.charupper * d.charlower / kh.charlower) ** 0.5
        lambda_ = (lambda_h / (Lf.charlower + B.charlower + lambda_h)) * np.exp((0.5 * B.charlower - xexit.charlower) / lambda_h)
        phi_exit = hp.charlower + lambda_*(h.charupper - hp.charlower)
        dphi = phi_exit - hp.charlower
        dphi_cu = d.charlower * ((ysat.charlower - yw.charupper) / yw.charupper)
        
        if(dphi < dphi_cu) :
            print("Pass")
        else:
            print("Fail")
            
        print("FoS:", dphi_cu/dphi)