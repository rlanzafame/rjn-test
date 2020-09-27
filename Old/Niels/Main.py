"""
Main.py
This is the main file where all parameters are set. The program should be started from this file.

@author: nvandervegt
"""

# Import own classes / modules
import Other.GenericFunctions as f;
import Modules.Uplift as Uplift;
import Modules.Heave as Heave;
#etc...


# Run calculations
# Uplift
uplift = Uplift.Uplift();
targetBetaUplift = f.ProbabilityToReliabilityIndex(1.30*10**-5)
uplift.SetBetaReq(targetBetaUplift);
uplift.DeterministicCalculation();

# Heave
heave = Heave.Heave();