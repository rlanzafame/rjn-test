#!/usr/bin/env python
# coding: utf-8

# <figure>
#   <IMG SRC="https://raw.githubusercontent.com/mbakker7/exploratory_computing_with_python/master/tudelft_logo.png" WIDTH=250 ALIGN="right">
# </figure>
# 
# # CIE5314 Flood defences 
# ### Exercise April 2020
# *Developed: HKV
# 
# ## 1	Import
# In chapter 2 starts the exercise description, in this chapter data en functions to be used in the exercise are imported. 

# In[3]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openturns as ot #run pip install openturns in the cmd kernal
from scipy.interpolate import interp1d
from scipy import stats as st
from tabulate import tabulate
# get_ipython().run_line_magic('matplotlib', 'inline')


# ### 1.1 Import data

# In[4]:


data_loc='C:/Users/rlanzafame/Desktop/Exercise/200401 HKV Assignment_attachment/Assignment_in_Python/' #You have to change this to the data folder
cross_sections=pd.read_excel(data_loc+'cross_sections.xlsx')
wind=pd.read_excel(data_loc+'wind.xlsx')
wind_exceedence=pd.read_excel(data_loc+'Exceedence_probability_wind_velocity.xlsx',index_col=0)
waterlevel_CDF=np.loadtxt(data_loc+'Cumulative_density_function_water_level.xpt')


#%%


def plot_cross_sections (name,create_fig=False, create_legend=False):
    if create_fig:
        plt.figure(figsize=(15,5))
        plt.grid()
        plt.xlabel('Distance [m]',fontsize=16)
        plt.ylabel('Surface level [m + NAP]',fontsize=16)
    plt.plot(cross_sections[name+', x'],cross_sections[name+', z'], label='cross section: '+name)
    if create_legend:
        plt.legend(fontsize=12)

def CDF_to_PDF_to_distribution(CDF,N,delta=0.01,plot_cdf=True,plot_pdf=True,legend=True):
    cdf_x=CDF[:,0]
    cdf_y=CDF[:,1] 
    x_hist=np.arange(np.min(cdf_x),np.max(cdf_x),delta)
    func_hist=interp1d(cdf_x,cdf_y)
    y_hist_cdf=func_hist(x_hist)
    y_hist_pdf=(y_hist_cdf[1::]-y_hist_cdf[0:-1])/(x_hist[1::]-x_hist[0:-1])
    hist_dist=st.rv_histogram((y_hist_pdf,x_hist))
    if plot_cdf:
        plt.plot(x_hist,hist_dist.cdf(x_hist),label='CDf')
    if plot_pdf:
        plt.plot(x_hist,hist_dist.pdf(x_hist), label='PDF')
    if legend:
        plt.legend()
        plt.xlabel('Waterlevel')
        plt.title('Distribution based on data')
    return hist_dist.ppf(np.random.random(N))


#%% OpenTurns example


# N=10000 #sample size
# x= ot.Sample(N,1)  #uses a dataset

# for i in range(N):
#     x[i,0]=CDF_to_PDF_to_distribution(waterlevel_CDF,N=1,delta=0.01,plot_cdf=False,plot_pdf=False,legend=False)[0]
    
# x_dis=ot.UserDefined(x,[1/N]*N) # distribution based on dataset
# a_dis=ot.Normal(20,10) # normal distribution with mean=20, standarddeviation=10
# b_dis=ot.LogNormalMuSigma(10,300).getDistribution() # lognormal distributionm with mean=10 and standarddeviation=300

# # Z function
# def example_Z_function(X):
#     x,a,b=X
#     Z=a*x**2-b
#     return [Z]

# # Create the event you want to estimate the probability
# marginals=[x_dis,a_dis,b_dis]
# description=['x','a','b']
# Z=ot.PythonFunction(len(marginals),1,example_Z_function)
# RS = ot.CorrelationMatrix(len(marginals))
# R = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(RS)
# copula = ot.NormalCopula(R)
# distribution = ot.ComposedDistribution(marginals, copula)
# distribution.setDescription(description)
# vect = ot.RandomVector(distribution)
# G = ot.CompositeRandomVector(Z, vect)
# event = ot.ThresholdEvent(G, ot.Less(), 0.0)

# #Define a solver
# optimAlgo = ot.Cobyla()
# optimAlgo.setMaximumEvaluationNumber(100000)
# optimAlgo.setMaximumAbsoluteError(1.0e-10)
# optimAlgo.setMaximumRelativeError(1.0e-10)
# optimAlgo.setMaximumResidualError(1.0e-10)
# optimAlgo.setMaximumConstraintError(1.0e-10)
   
# #run Form
# startvalues=distribution.getMean()
# startvalues[0]=7 #Let the waterlevel start at for example 7
# algo = ot.FORM(optimAlgo, event, distribution.getMean())
# algo.run()
# result = algo.getResult()
# print ('Form probability:',round(result.getEventProbability(),2))

# #run Monte Carlo Simulation
# experiment = ot.MonteCarloExperiment()
# algo = ot.ProbabilitySimulationAlgorithm(event, experiment)
# algo.setMaximumCoefficientOfVariation(0.05)
# algo.setMaximumOuterSampling(int(N))
# algo.run()
# result=algo.getResult()

# print ('Montecarlo probability',round(result.getProbabilityEstimate(),2 ))




qc_mean= 100 #l/s/m 


# #### Question 9
# Calculate the characteristic values with 5% or 95% exceedance probability to be used in the overtopping and overflow calculation.  

#%%


def characteristic_normal (mean,standard,percentage):
    if percentage==5:
        return mean-1.64*standard
    elif percentage==95:
        return mean +1.64*standard
    else:
        print ('warning') 
 
def characteristic_lognormal(mean,standard,percentage):
    sigmam=np.sqrt(np.log(1+(standard/mean)**2))
    mum=np.log(mean)-0.5*sigmam**2
    tn=1.76
    n=15
    a=0
    if percentage==5:
        return np.exp(mum-tn*sigmam*np.sqrt((1-a)+1/n))
    elif percentage==95:
        return np.exp(mum+tn*sigmam*np.sqrt((1-a)+1/n))
    else:
        print ('warning') 
d_mean=2.35
d_sigma=0.3
d_char1=characteristic_normal(d_mean,d_sigma,5)
qc_sigma=120
qc_char=characteristic_lognormal(qc_mean,qc_sigma,5)
print ('Characteristic value for the bottom level is',round(d_char1,2), 'm')
print ('Characteristic value for the critical discharge is',round(qc_char,2), 'm')


# #### Question 10: 
# Is the safety level requirement met with the semi-probabilistic calculation?
# 
# ---
# **Answer**
# 
# `Section 5.3.2, 5.4  CIE5314 Flood defence`<br/>
# The amount of overtopping can be calculated using Van der Meer and Bruce (2014) formula. In which the significant wave height and period are calculated using Young and Verhagen. The freeboard is equal to 0.04 m. No berm or wall is present ($γ_b=1.0$ and $γ_v=1.0$) ) and the slope consist out of grass ($γ_f=1.0$). Now the overtopping discharge can be calculated using the formulas given below:

#%%


plot_cross_sections ('BF098_1',create_fig=True)
xmin=222
xmax=247
hcrest=np.max(cross_sections['BF098_1, z'])
xdata=cross_sections['BF098_1, x'][xmin:xmax]
ydata=cross_sections['BF098_1, z'][xmin:xmax]
a = st.linregress(xdata,ydata)
slope=1/-a[0] 
yfit=a[0]*xdata+a[1]
plt.plot(cross_sections['BF098_1, x'][np.argmax(cross_sections['BF098_1, z'])],hcrest,'o', label='crestheight')
plt.plot(xdata,yfit,color='black',linestyle='--',label='Inner slope')
plt.legend()
print('The slope is 1:',round(slope,1))
print ('The maximum height is',round(hcrest,2),'m')


#%%


# input
hwaterlevel=7.56
g=9.81
anglenormal=12
yb=1
yf=1
yv=1
depth=hwaterlevel-d_char1
OTcoeff='kar'
Rc=hcrest-hwaterlevel # 8.447
def overtopping(direction,u,F,d=depth,Rc=Rc,slope=slope,anglenormal=anglenormal,yb=yb,yf=yf,yv=yv,OTcoeff=OTcoeff):
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
    
    Hm0=Hs
    Tm1=Tp/1.1
    
    # adjust for wind direction
    Hm0=Hm0*betaHm0
    Tm1=Tm1*betaTm1
    
    # Lm1 is the deep water wave length (can find in Eurotop)
    # Lm1 = g*Tm1**2/2/pi, or 1.5614*Tm1**2
    Lm1=g*Tm1**2/2/3.1415
    
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
        
        qm2=otBm*Rc/(yf*yB*Hm0)
        
        qVdM=q0*q1*np.exp(-q2**1.3)*1000
        qmax=q0*otAm*np.exp(-qm2**1.3)*1000
        q=min(qVdM,qmax)
        
    OTdirDEP=np.array([direction,beta,u,F,Hm0,Tm1,Irr,q])
    return q, OTdirDEP

def Calculation_overflow(Direction,F,Uwind,hbottom=d_char1,hwaterlevel=hwaterlevel,g=g,hcrest=hcrest,yb=yb,yf=yf,yv=yv,anglenormal=anglenormal,slope=slope):
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
    if Rc>=0:
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
    q_directionRL,q_rows[i,0:8]=overtopping(Direction,Uwind,F)
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


## convert your array into a dataframe
df = pd.DataFrame (q_rows)
## save to xlsx file
# filepath = 'q_rows.xlsx'
# df.to_excel(filepath, index=False)

#%%


# print(tabulate([['Alice', 24], ['Bob', 19]], headers=['Name', 'Age']))

#%%


# # Fit for the wind
# from scipy.stats import weibull_min
# from scipy.optimize import curve_fit

# def fitweibb(x,c,scale): #x=wind
#     return 1/(1-weibull_min.cdf(x,c,0,scale)) #returns the return period
# winddirection='WNW'
# column=wind_exceedence[winddirection]
# xfit=np.linspace(wind_exceedence.index.min(),wind_exceedence.index.max(),100) #windspeed
# y=1/column
# fitpar=curve_fit(fitweibb,column.index,y,p0=[2,10])
# yfit=fitweibb(xfit,fitpar[0][0],fitpar[0][1])

# # Plot     
# plt.figure(figsize=(15,5))
# plt.suptitle('Wind direciton '+winddirection, size='20')
# plt.subplot(133)
# plt.plot(xfit,weibull_min.pdf(xfit,fitpar[0][0],0,fitpar[0][1]))
# plt.title('PDF')
# plt.xlabel('Wind velocity')
# plt.ylabel('Probability')

# plt.subplot(131)
# plt.semilogx(y,column.index,label='data')
# plt.semilogx(yfit,xfit,label='fit')
# plt.title('Frequency line')
# plt.xlabel('Return period')
# plt.ylabel('Wind velocity')
# plt.legend()

# plt.subplot(132)
# plt.plot(column.index,1-1/y,label='data')
# plt.plot(xfit,1-1/yfit,label='fit')
# plt.title('CDF')
# plt.xlabel('Wind velocity')
# plt.ylabel('Probability')
# plt.legend()    


# #%%


# #Distribution
# N=100000
# N_array=[1/N]*N
# hwaterlevel_dis= ot.Sample(N,1)
# Uwind_dis= ot.Sample(N,1)
# for i in range(N):
#     hwaterlevel_dis[i,0]= CDF_to_PDF_to_distribution(waterlevel_CDF,N=1,delta=0.01,plot_cdf=False,plot_pdf=False,legend=False)[0]
#     Uwind_dis[i,0]= weibull_min.ppf(np.random.random(1),fitpar[0][0],0,fitpar[0][1])[0]
# h_bottom_dist=ot.Normal(d_mean,d_sigma)
# hwaterlevel_dist=ot.UserDefined(hwaterlevel_dis,N_array)
# Uwind_dist=ot.UserDefined(Uwind_dis,N_array)
# qcrit_dist= ot.LogNormalMuSigma(qc_mean,qc_sigma).getDistribution()


# #%%


# #Z-function
# def Z_overflow(X):
#     h_bottom,h_waterlevel,Uwind,qcrit=X
#     Z=qcrit-Calculation_overflow(Direction,F,Uwind,hbottom=h_bottom,hwaterlevel=h_waterlevel,g=g,hcrest=hcrest,yb=yb,yf=yf,yv=yv,anglenormal=anglenormal,slope=slope)
#     return [Z]

# # Create the event we want to estimate the probability
# marginals=[h_bottom_dist,hwaterlevel_dist,Uwind_dist,qcrit_dist]
# description=['Bottom height','Waterlevel','Wind veolocity', 'Critical overtopping discharge']
# Z=ot.PythonFunction(len(marginals),1,Z_overflow)
# RS = ot.CorrelationMatrix(len(marginals))
# R = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(RS)
# copula = ot.NormalCopula(R)
# distribution = ot.ComposedDistribution(marginals, copula)
# distribution.setDescription(description)
# vect = ot.RandomVector(distribution)
# G = ot.CompositeRandomVector(Z, vect)
# event = ot.ThresholdEvent(G, ot.Less(), 0.0)

# #Define a solver
# optimAlgo = ot.Cobyla()
# optimAlgo.setMaximumEvaluationNumber(100000)
# optimAlgo.setMaximumAbsoluteError(1.0e-10)
# optimAlgo.setMaximumRelativeError(1.0e-10)
# optimAlgo.setMaximumResidualError(1.0e-10)
# optimAlgo.setMaximumConstraintError(1.0e-10)
        
# #run MC
# experiment = ot.MonteCarloExperiment()
# algo = ot.ProbabilitySimulationAlgorithm(event, experiment)
# algo.setMaximumCoefficientOfVariation(0.05)
# algo.setMaximumOuterSampling(int(1e6))
# algo.run()
# result=algo.getResult() 
# pf_overtopping=result.getProbabilityEstimate()



