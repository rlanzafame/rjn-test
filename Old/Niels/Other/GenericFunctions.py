"""
GenericFunctions.py
Some generic functions

@author: nvandervegt
"""

# Import libraries
import scipy.stats as sc;


# Calculates the reliability index from a probability
def ProbabilityToReliabilityIndex(_pf):
    beta = sc.norm.isf(_pf);
    return beta;


# Calculates the probability from a reliability index
def ReliabilityIndexToProbability(_beta):
    pf = sc.norm.cdf(-_beta);
    return pf;