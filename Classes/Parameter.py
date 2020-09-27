"""
Parameter.py
File for the parameter class

@author: nvandervegt

A parameter class/object can be defined as:
k = Parameter(distribution, mean, stdev=0)
And has the following functions:
k.mean -> gives mean value
k.stdev -> gives stdev
k.distribution -> gives distribution (string)
k.upperchar -> gives 5% excedaance characteristic value
k.lowerchar -> gives 95% exceedance characteristic value

The characteristic values will automatically be calculated when a parameter is initiated
The student t-factor, a and number of observation can be changed by using:
    k.tn_1
    k.a
    k.n
After changing, the command:
    k.RecalculateCharacteristics();
Has to be executed to recalculate the characteristic values with the new settings
"""

# Import
import numpy as np;


# Parameter class
class Parameter:
    
    # Init
    def __init__(self, distribution, mean, stdev=0):
        self.distribution = distribution;
        self.mean = mean;
        self.stdev = stdev;
        self.tn_1 = 1.76; # Student t-factor corresponding with the coincidence bound, Assignment 2020
        self.a = 0; # Ratio between the local and regional variation ,Assignment 2020
        self.n = 15; # Number of observations, Assignment 2020
        self.RecalculateCharacteristics();
        
    
    # Calculates the 5%/95% characteristic values
    # This function is automatically called when a parameter is created
    def RecalculateCharacteristics(self):
        if(self.distribution.lower() == 'deterministic'):
            self.charupper = self.mean;
            self.charlower = self.mean;
        elif(self.distribution.lower() == 'normal'):
            self.charupper = self.mean + 1.645 * self.stdev;
            self.charlower = self.mean - 1.645 * self.stdev;
        elif(self.distribution.lower() == 'lognormal'):
            sigmam = (np.log(1 + (self.stdev / self.mean) ** 2)) ** 0.5;
            um = np.log(self.mean) - 0.5 * sigmam ** 2;
            self.charupper = np.exp(um + self.tn_1 * sigmam * ((1 - self.a) + (1 / self.n)) ** 0.5);
            self.charlower = np.exp(um - self.tn_1 * sigmam * ((1 - self.a) + (1 / self.n)) ** 0.5);