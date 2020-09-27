"""
Input.py
Input parameters, define the parameters by using the parameter class

@author: nvandervegt
"""

# Import the Parameter class
from Classes.Parameter import Parameter;


# Question 14 - Internal Stability
# Question 14 - Uplift
k = Parameter('lognormal', 7.52 * 10**-4, 3.76 * 10**-4);
kh = Parameter('lognormal', 1.00 * 10**-6, 5.00 * 10**-7);
D = Parameter('lognormal', 25, 0.5);
d = Parameter('lognormal', 4.5, 0.5);
Lf = Parameter('lognormal', 20, 2);
B = Parameter('deterministic', 37.5);
xexit = Parameter('deterministic', 32.5);
h = Parameter('deterministic', 7.5);
hp = Parameter('normal', 3.50, 0.1);
ysat = Parameter('normal', 14.72, 0.736);
yw = Parameter('deterministic', 10);
mu = Parameter('normal', 1, 0.10);

# Question 14 - Heave
ich = Parameter('lognormal', 0.50, 0.1);

# Question 14 - Piping
# Etc..