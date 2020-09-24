# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 15:34:53 2020

@author: rlanzafame
"""


# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 15:12:51 2020

@author: rlanzafame
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openturns as ot #run pip install openturns in the cmd kernal
from scipy.interpolate import interp1d
from scipy import stats as st
from tabulate import tabulate

#%% Parameters

# these are only used for defining Rc and d
# values for exercise:
# hcrest=7.59
# hwaterlevel=7.56
# hbottom=1.86
# # Rc and d are used directly as inputs
# Rc=hcrest-hwaterlevel # freeboard
# d=hwaterlevel-hbottom # water depth

Rc=0.03 # freeboard
d=5.70 # water depth

slope=3.5
anglenormal=12 # bearing of dike normal (pointing away from dike)

direction=12 # wind direction
F=2848
u=34.8 # wind speed at 10 m in m/s

yb=1
yf=1
yv=1

OTcoeff='mean' # use 'mean' or 'kar'

#%% Calculations

g=9.81

# parameters for OT eqn (normal and max q)
#    mu sd kar:
# A   .023  .003  .026
# B  2.7   0.2   2.5
# Am  .09   .013  .103
# Bm 1.5    .15  1.35
if OTcoeff=='mean':
    otA =  .023
    otB = 2.7
    otAm=  .09
    otBm= 1.5
elif OTcoeff=='kar':
    otA =  .026
    otB = 2.5
    otAm=  .103
    otBm= 1.35

alpha=np.arctan(1/slope)

if direction-anglenormal>180:
    beta = abs(direction-anglenormal-360)
elif direction-anglenormal <=180: 
    beta = abs(direction-anglenormal)

yB=max([1-0.0033*abs(beta),1-0.0033*80])
if beta<80:
    betaHm0=1
    betaTm1=1
elif beta>=110:
    betaHm0=0
    betaTm1=0
else:
    betaHm0=(110-beta)/30
    betaTm1=np.sqrt((110-beta)/30)


# Bretschneider equations
 
Ftilde=F*g/u**2
dtilde=d*g/u**2

H1=np.tanh(0.530*dtilde**0.75)
H2=np.tanh(0.0125*Ftilde**0.42/H1)
Htilde=0.283*H1*H2

T1=np.tanh(0.833*dtilde**0.375)
T2=np.tanh(0.077*Ftilde**0.25/T1)
Ttilde=2.4*3.1415*T1*T2

Hs=Htilde*u**2/g
Ts=Ttilde*u/g

Tp=1.08*Ts

# find spectral values
Hm0=Hs
Tm1=Tp/1.1

# adjust for wind direction
Hm0=Hm0*betaHm0
Tm1=Tm1*betaTm1

# Lm1 is the deep water wave length (can find in Eurotop)
# Lm1 = g*Tm1**2/2/pi, or 1.5614*Tm1**2
Lm1=g*Tm1**2/2/3.1415

if beta>=110:
    Irr,qVdM,qmax,q,qOF=0,0,0,0,0
else:
    # Irribarren number:
    #   spilling < plunging < surging/collapsing
    #           0.5         3.3     deep water
    #           0.4         2.0     break point
    # EuroTop: waves are breaking for Irr < 2-3
    Irr=np.tan(alpha)/np.sqrt(Hm0/Lm1)
    
    q0=np.sqrt(g*Hm0**3)
    q1=otA/np.sqrt(np.tan(alpha))*yb*Irr
    q2=otB*Rc/(yb*yf*yB*yv*Irr*Hm0)
    
    qm2=otBm*Rc/(yf*yB*Hm0)
    
    qVdM=q0*q1*np.exp(-q2**1.3)*1000
    qmax=q0*otAm*np.exp(-qm2**1.3)*1000
    q=min(qVdM,qmax)
    qOF=0.54*np.sqrt(g*(abs(Rc))**3)*1000 

#%% Output

print('Bretschneider:')
print('\n Ftilde =',round(Ftilde,3),' dtilde =',round(dtilde,3))
print(' Htilde =',round(Htilde,3),' Ttilde =',round(Ttilde,3))
print(' Hs =',round(Hs,3),'m, Ts =',round(Ts,3),'s, Tp =',round(Tp,3),'s')

print('\nInputs for van der Meer:')
print('\nHm0 =',round(Hm0,3),'m, Tm1 =',round(Tm1,3),'s, Lm1 =',round(Ts,3),'m')
print('\nIrr =',round(Irr,3),'m, beta =',round(beta,3),'s, yB =',round(yB,3),'m')
print('\nRc =',round(Rc,3),'m, slope =',round(slope,2),', depth =',round(d,3),'m')

print('\nOvertopping and overflow:')
print('\n',round(qVdM,1),'l/m/s = overtopping, van der Meer')
print('\n',round(qmax,1),'l/m/s = overtopping, maximum')
print('\n',round(q,1),'l/m/s = overtopping minimum')
print('\n',round(qOF,1),'l/m/s = overflow')

OTdirIND=np.array([d,Rc,slope,anglenormal,yb,yf,yv,OTcoeff])
OTdirDEP=np.array([direction,beta,u,F,Hm0,Tm1,Irr,q])

OTdirINDheaders=['d (m)','Rc (m)','slope (-)','norm (d)','yb (-)','yf (-)','yv (-)','q coeff']
OTdirDEPheaders=['dir (d)','beta (d)','u (m/s)','F (m)','Hm0 (m)','Tm1 (s)','Irr (-)','q (l/m/s)']
print('\n\n')
print(tabulate([OTdirIND],OTdirINDheaders))
print('\n')
print(tabulate([OTdirDEP],OTdirDEPheaders))

# print(tabulate([[OTdirDEP],[OTdirDEP]],OTdirDEPheaders))
