### GENERAL PACKAGES ###
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openturns as ot #run pip install openturns in the cmd kernal
from scipy.interpolate import interp1d
from scipy import stats as st
# get_ipython().run_line_magic('matplotlib', 'inline')
from tabulate import tabulate

### PACKAGES FOR QUESTION 9 ###
from scipy import stats as st
from scipy import special
from scipy.stats import lognorm

## PREPARATORY CALCULATIONS ##

## DEFINITIONS ##
def overtopping(wind_dir,u,F,Rc,d=d_water,slope=alpha,anglenormal=anglenormal,Lberm = Lberm,b_berm = b_berm ,yf=yf,yv=yv,OTcoeff=OTcoeff):   
   
    # parameters for OT eqn (normal and max q)
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

    #ALTERNATIV #
    """muA11 = 0.023
    sigmaA11 = 0.003
    muA12 = 0.09
    sigmaA12 = 0.013

    muB11 = 2.7
    sigmaB11 = 0.2
    muB12 = 1.5
    sigmaB12 = 0.15

    A11 = st.norm.ppf(0.95,0.023,0.003)
    A12 = st.norm.ppf(0.95,0.09,0.013)

    B11 = st.norm.ppf(0.05,2.7,0.2)
    B12 = st.norm.ppf(0.05,1.5,0.15)"""
    
    # obliquety reduction factors 
    if wind_dir-anglenormal>180:
        beta = abs(wind_dir-anglenormal-360)
    elif wind_dir-anglenormal <=180: 
        beta = abs(wind_dir-anglenormal)
        
    yB=max([1-0.0033*beta,1-0.0033*80])
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
     
    wind_fetch_tilde=wind_fetch*g/wind_speed**2
    d_water_tilde=d_water*g/wind_speed**2
    
    H1=np.tanh(0.530*d_water_tilde**0.75)
    H2=np.tanh(0.0125*wind_fetch_tilde**0.42/H1)
    Htilde=0.283*H1*H2
    
    T1=np.tanh(0.833*d_water_tilde**0.375)
    T2=np.tanh(0.077*wind_fetch_tilde**0.25/T1)
    Ttilde=2.4*3.1415*T1*T2
    
    Hs=Htilde*u**2/g
    Ts=Ttilde*u/g
    
    Tp=1.08*Ts
    
    Hm0=Hs
    Tm1=Tp/1.1
    
    # adjust for wind direction
    Hm0=Hm0*betaHm0
    Tm1=Tm1*betaTm1

    # influence of berm
    rb = b_berm / Lberm
    rdb = 0.5 - 0.5 * cos(np.pi * freeboard_berm / (2 * Hm0))
    yb = min(1 - rb * (1 - r_db),0.6)
    
    # Lm1 is the deep water wave length (can find in Eurotop)
    # Lm1 = g*Tm1**2/2/pi, or 1.5614*Tm1**2
    Lm1=g*Tm1**2/2/np.pi
    
    if beta>=110:
        Irr,qVdM,qmax,q=0,0,0,0
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
        
        qm2=otBm*freeboard_dike/(yf*yB*Hm0)
        
        qVdM=q0*q1*np.exp(-q2**1.3)*1000
        qmax=q0*otAm*np.exp(-qm2**1.3)*1000
        q=min(qVdM,qmax)

def Calculation_overflow(Direction = wind_dir,F = wind_fetch,Uwind = wind_speed,hbottom=d_water,hwaterlevel=h_waterlevel,g=g,hcrest=h_crest,yb=yb,yf=yf,yv=yv,anglenormal=anglenormal,slope=slope):
    Rc=hcrest-hwaterlevel
    alpha=np.arctan(1/slope)
    beta=Direction-anglenormal
    if beta>180:
        beta=beta-360
    ybeta=1-0.0033*abs(beta) 
    Depth=hwaterlevel-hbottom
    if Depth<=0:
        Depth=0.01
    Ftilde=g*F/Uwind**2
    dTilde=g*Depth/Uwind**2
    hm0first=0.343*dTilde**1.14
    hm0second=(4.41*10**(-4)*Ftilde**0.79)/np.tanh(hm0first)
    tm0first=0.1*dTilde**2.01
    tm0second=(2.77*10**(-7)*Ftilde**1.45)/np.tanh(tm0first)
    Hm0=0.24*Uwind**2/g*(np.tanh(hm0first)*np.tanh(hm0second))**0.572
    Tp=7.69*Uwind/g*(np.tanh(tm0first)*np.tanh(tm0second))**0.187
    Xim=np.tan(alpha)/np.sqrt(Hm0/(1.56*(0.9*Tp)**2))
    if Rc>0:
        q=(0.023/np.sqrt(np.tan(alpha)))*Xim*np.exp(-((2.7*Rc/(yb*yf*ybeta*yv*Xim*Hm0))**1.3))*1000*np.sqrt(g*Hm0**3)
        qmax=np.sqrt(g*Hm0**3)*0.09*np.exp(-((1.5*Rc/Hm0)**1/3))*1000
        q=min(q,qmax)
    else:
        q=0.54*np.sqrt(g*(abs(Rc))**3)*1000 
    return q

q_tot=0
q_rows = np.zeros([len(wind),10])
for i in range(len(wind)):
    Direction=wind.iloc[:,0][i]
    F=wind.iloc[:,2][i]
    Uwind=wind.iloc[:,3][i]
    pf=wind.iloc[:,4][i]
    q_direction=Calculation_overflow(Direction,F,Uwind)
    print ('Wind direction:', Direction,'deg, F=',F,'m, u=', Uwind,'m/s --> q=',round(q_direction,3),'l/m/s')
    q_tot+=q_direction*pf/100
    q_directionRL,q_rows[i,0:8]=overtopping(Direction,Uwind,F,Rc)
    q_rows[i,8:]=[pf/100,q_rows[i,7]*pf/100]
    
print ('\n total q is ',round(q_tot,1),'l/m/s, critical value, characteristic:', round(qc_char,1), 'l/s/m')

OTdirINDheaders=['d (m)','Rc (m)','slope (-)','norm (d)','yb (-)','yf (-)','yv (-)','q coeff']
OTdirDEPheaders=['dir (d)','beta (d)','u (m/s)','F (m)','Hm0 (m)','Tm1 (s)','Irr (-)','q (l/m/s)','p (-)','q*p (l/m/s)']

OTdirIND=np.array([depth,Rc,slope,anglenormal,yb,yf,yv,OTcoeff])

print('\n\n')
print(tabulate([OTdirIND],OTdirINDheaders))
print('\n')
print(tabulate(q_rows,OTdirDEPheaders))
print ('\n Total expected  value for q is ',round(sum(q_rows[:,-1]),1),'l/m/s, critical value, characteristic:', round(qc_char,1), 'l/s/m')
