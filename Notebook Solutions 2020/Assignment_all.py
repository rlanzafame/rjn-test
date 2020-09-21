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
# get_ipython().run_line_magic('matplotlib', 'inline')
from tabulate import tabulate

# ### 1.1 Import data

# In[4]:


data_loc='C:/Users/rlanzafame/Desktop/Exercise/200401 HKV Assignment_attachment/Assignment_in_Python/' #You have to change this to the data folder
cross_sections=pd.read_excel(data_loc+'cross_sections.xlsx')
wind=pd.read_excel(data_loc+'wind.xlsx')
wind_exceedence=pd.read_excel(data_loc+'Exceedence_probability_wind_velocity.xlsx',index_col=0)
waterlevel_CDF=np.loadtxt(data_loc+'Cumulative_density_function_water_level.xpt')


# ### 1.1 Import functions

# In[5]:


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


# ### 1.2 Example probabilistic calculation in Python using the OpenTurns module

# In[6]:


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




mortality=1.19/100 # 1.19 percent
LIR=1/100000
degree_of_evacuation=0.55

Ptraject=LIR/(mortality*(1-degree_of_evacuation))


#input
Preq=1/10000
w_overtopping=0.24
w_erosion=0.24
w_slope=0.04

#calculation
P_overtopping=w_overtopping*Preq
P_erosion=w_erosion*Preq
P_slope=w_slope*Preq


#functions 
def Length_effect(a,b,L):
    return 1+a*L/b
def target_prob (Preq,N):
    return Preq/N

#Input
N_overtopping=1
Ltraject=15860

#calculation
N_erosion=Length_effect(0.9,300,Ltraject)
N_slope=Length_effect(0.033,50,Ltraject)
Ptarget_overtopping=target_prob(P_overtopping,N_overtopping)
Ptarget_erosion=target_prob(P_erosion,N_erosion)
Ptarget_slope=target_prob(P_slope,N_slope)



qc_mean= 100 #l/s/m 



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



# In[18]:


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


# In[42]:


# input
hwaterlevel=7.59
g=9.81
anglenormal=12
yb=1
yf=1
yv=1
depth=hwaterlevel-d_char1
OTcoeff='kar'
Rc=hcrest-hwaterlevel
# Rc=0.03
def overtopping(direction,u,F,Rc,d=depth,slope=slope,anglenormal=anglenormal,yb=yb,yf=yf,yv=yv,OTcoeff=OTcoeff):
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


## convert your array into a dataframe
df = pd.DataFrame (q_rows)
## save to xlsx file
# filepath = 'q_rows.xlsx'
# df.to_excel(filepath, index=False)


# In[20]:


# Fit for the wind
from scipy.stats import weibull_min
from scipy.optimize import curve_fit

def fitweibb(x,c,scale): #x=wind
    return 1/(1-weibull_min.cdf(x,c,0,scale)) #returns the return period
winddirection='WNW'
column=wind_exceedence[winddirection]
xfit=np.linspace(wind_exceedence.index.min(),wind_exceedence.index.max(),100) #windspeed
y=1/column
fitpar=curve_fit(fitweibb,column.index,y,p0=[2,10])
yfit=fitweibb(xfit,fitpar[0][0],fitpar[0][1])

# Plot     
plt.figure(figsize=(15,5))
plt.suptitle('Wind direciton '+winddirection, size='20')
plt.subplot(133)
plt.plot(xfit,weibull_min.pdf(xfit,fitpar[0][0],0,fitpar[0][1]))
plt.title('PDF')
plt.xlabel('Wind velocity')
plt.ylabel('Probability')

plt.subplot(131)
plt.semilogx(y,column.index,label='data')
plt.semilogx(yfit,xfit,label='fit')
plt.title('Frequency line')
plt.xlabel('Return period')
plt.ylabel('Wind velocity')
plt.legend()

plt.subplot(132)
plt.plot(column.index,1-1/y,label='data')
plt.plot(xfit,1-1/yfit,label='fit')
plt.title('CDF')
plt.xlabel('Wind velocity')
plt.ylabel('Probability')
plt.legend()    


# In[21]:


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


# In[22]:


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

# print ('The failure probability for the wind direction '+winddirection,'is: 1/',round(1/pf_overtopping), 'per year, which means the dike meets the required safety level of : 1/',round(1/Ptarget_overtopping),'per year')



# In[24]:


print ('Cross-section BF097, because the leakage length is smallest (smaller dike width). ')
plot_cross_sections ('BF097',create_fig=True)
xmin=74
xmax=114
height=3.6
plt.plot([xmin,xmax],[height,height],color='black',linestyle='--',marker='^', label='Dike width')
plt.legend()
B=xmax-xmin
print ('The widht of the this cross section is:', round(B),'m')


# #### Question 14
# Calculate the characteristic values with 5% or 95% exceedance probability to be used in the internal erosion calculation. 
# 
# ---
# 
# **Answer**
# From table 3 and figure 8: <br/>
# Thickness peat: 2.45 m ($\gamma=14.5kN/m$)<br/>
# Thickness clay: 2.00 m ($\gamma=15.0 kN/m³$)

# In[25]:


#Input
k_mean=7.52e-4 #hydraulic conductivity aquifer
k_sig=0.5*k_mean
D_mean=25 # Aquifer thickness
D_sig=0.5
d_mean= 2.45+2.00# blanket thickness
d_sig=0.5
kh_mean= 1.00e-6#hydraulic conductivity aquitard
kh_sig=0.5*kh_mean
Lf_mean= 20#length of effective foreshore
Lf_sig=0.1*Lf_mean
ysat_mean= (2.45*14.5+2*15)/d_mean #saturated volumetric weight of blanket, y
ysat_sig=0.05*ysat_mean
mu_mean=1.0 #model factor critical head differenc3e
mu_sig=0.1
ic_mean= 0.5#critical exit gradient
ic_sig=0.1
d70_mean=2.8e-4
d70_sig=0.12*d70_mean
mp_mean=1.0 #model uncertainty factor
mp_sig=0.12
hp_mean= 3.5#hinterland level
hp_sig =0.1 

#Calculation
k_char=characteristic_normal(k_mean,k_sig,95)
D_char=characteristic_lognormal(D_mean,D_sig,95)
d_char=characteristic_lognormal(d_mean,d_sig,5)
kh_char=characteristic_normal(kh_mean,kh_sig,95)
Lf_char=characteristic_lognormal(Lf_mean,Lf_sig,5)
hp_char=characteristic_normal(hp_mean,hp_sig,5)
ysat_char=characteristic_normal(ysat_mean,ysat_sig,5)
mu_char=characteristic_normal(mu_mean,mu_sig,5)
ic_char=characteristic_lognormal(ic_mean,ic_sig,5)
d70_char=characteristic_lognormal(d70_mean,d70_sig,5)
mp_char=characteristic_normal(mp_mean,mp_sig,5)

print ('Characteristic value for the Hydraulic conductivity of the aquifer is',round(k_char,5), 'm/s')
print ('Characteristic value for the Aquifer thickness is',round(D_char,2), 'm')
print ('Characteristic value for the Blanket thickness at the exit point thickness is',round(d_char,2), 'm')
print ('Characteristic value for the Hydraulic conductivity is',round(kh_char,10), 'm/s')
print ('Characteristic value for the Length of the effective foreshorel is',round(Lf_char,2), 'm')
print ('Characteristic value for the Hinterland phreatic level is',round(hp_char,2), 'm')
print ('Characteristic value for the Saturated volumetric weight of the blanket is',round(ysat_char,2), 'kN/m³')
print ('Characteristic value for the Model factor addressing the uncertainty in the critical head difference is',round(mu_char,2))
print ('Characteristic value for the Critical exit gradient is',round(ic_char,2))
print ('Characteristic value for the 70%-fractile of the grain size distribution is',round(d70_char,6), 'm')
print ('Characteristic value for the Model uncertainty factor is',round(mp_char,2))


# #### Question 15
# Does the levee meet the required safety standard in this semi-probabilistic calculation? What are the safety factors and what is the failure probability? You may assume the mechanisms uplift, heave and piping are uncorrelated and use the Sellmeijer theory for the piping calculation. 
# 

# In[26]:


def  calculation_internal_erosion(hw,k,D,dc,kh,Lf,hp,ysat,mu,ich,d70,mp,B=B,norm=Ptarget_erosion,return_Probability=True,return_Z=True,print_information=True):
    yw=10
    n=0.25
    ys=26.5
    xexit=B/2
    theta=37/360*2*np.pi
    d70m=2.08E-04
    v=1.33E-06
    L=B+Lf
    kh=abs(kh)
    k=abs(k)
    if (hp+0.3*dc)>hw:  # or in log np.max(0, code) 
        hw=hp+0.30001*dc
    beta=-st.norm.ppf(norm, loc=0, scale=1)
    uplift_positive=((dc*(ysat-yw)/yw))   
    uplift_negative=hp+(np.sqrt(k*D*dc/kh)/(Lf+B+np.sqrt(k*D*dc/kh)*np.exp((B/2-xexit)/np.sqrt(k*D*dc/kh)))*(hw-hp))-hp
    if uplift_negative==0:
        uplift_negative=10e-5
    FoS_uplift=uplift_positive/uplift_negative
    Z_uplift=mu*uplift_positive-uplift_negative
    beta_uplift=-(np.log(np.max((FoS_uplift/0.48,10e-12)))+0.27*beta)/0.46
    prob_uplift=st.norm.cdf(beta_uplift,loc=0,scale=1)
    heave_positive=ich
    heave_negative=(((hp+(np.sqrt(np.max((0,k*D*dc/kh)))/(Lf+B+np.sqrt(np.max((0,k*D*dc/kh))))*np.exp((B/2-xexit)/np.sqrt(np.max((0,k*D*dc/kh)))))*(hw-hp))-hp)/dc)
    if heave_negative==0:
        heave_negative=10e-5
    FoS_heave=heave_positive/heave_negative
    Z_heave=heave_positive-heave_negative
    beta_heave=-(np.log(FoS_heave/0.37)+0.3*beta)/0.48
    prob_heave=st.norm.cdf(beta_heave,loc=0,scale=1)
    F1=(n*((ys/yw)-1)*np.tan((theta)))
    F2=(d70m/(v*k*(Lf+B)/g)**(1/3)*(d70/d70m)**(0.4))
    F3=(0.91*(D/(Lf+B))**(0.28/(((D/(Lf+B))**2.8)-1)+0.04))
    piping_positive=F1*F2*F3*L
    piping_negative=(hw-hp-0.3*dc)
    FoS_piping=piping_positive/piping_negative
    Z_piping=mp*piping_positive-piping_negative
    beta_piping=-(np.log(FoS_piping/1.04)+0.43*beta)/0.37
    prob_piping=st.norm.cdf(beta_piping,loc=0,scale=1)
    Z= np.max((Z_uplift,Z_heave,Z_piping))
    #print (Z,Z_uplift,Z_heave,Z_piping)
    probability=np.min((prob_uplift,prob_heave,prob_piping))
    if print_information:
        print ('The factor of Safety for uplift is',round(FoS_uplift,2),'and the failure probility: 1/',round(1/prob_uplift),'per year')
        print ('The factor of Safety for heave is',round(FoS_heave,2),'and the failure probility: 1/',round(1/prob_heave),'per year')
        print ('The factor of Safety for heave is',round(FoS_piping,2),'and the failure probility: 1/',round(1/prob_piping),'per year')
    if return_Probability*return_Z:
        return probability,Z
    if return_Z:
        return Z
    if return_Probability:
        return probability
    
hwaterlevel_10000=7.35 
    
P_internal_erosion=calculation_internal_erosion(hw=hwaterlevel_10000,k=k_char,D=D_char,dc=d_char,kh=kh_char,Lf=Lf_char,hp=hp_char,ysat=ysat_char,mu=mu_char,ich=ic_char,d70=d70_char,mp=mp_char,return_Z=False)
print ('The failure probability is: 1/', round(1/P_internal_erosion), 'year. This does not meet the requirement of 1/',round(1/Ptarget_erosion), 'per year')




# In[29]:


# N=100000
# N_array=[1/N]*N
# hwaterlevel_dis= ot.Sample(N,1)

# for i in range(N):
#     hwaterlevel_dis[i,0]=CDF_to_PDF_to_distribution(waterlevel_CDF,1,delta=0.01,plot_cdf=False,plot_pdf=False,legend=False)[0]
    
# k_dist=ot.Normal(k_mean,k_sig)
# hwaterlevel_dist=ot.UserDefined(hwaterlevel_dis,N_array)
# D_dist=ot.LogNormalMuSigma(D_mean,D_sig).getDistribution()
# kh_dist=ot.Normal(kh_mean,kh_sig)
# d_dist=ot.LogNormalMuSigma(d_mean,d_sig).getDistribution()
# Lf_dist=ot.LogNormalMuSigma(Lf_mean,Lf_sig).getDistribution()
# hp_dist=ot.Normal(hp_mean,hp_sig)
# ysat_dist=ot.Normal(ysat_mean,ysat_sig)
# mu_dist=ot.Normal(mu_mean,mu_sig)
# ic_dist=ot.LogNormalMuSigma(ic_mean,ic_sig).getDistribution()
# d70_dist=ot.LogNormalMuSigma(d70_mean,d70_sig).getDistribution()
# mp_dist=ot.Normal(mp_mean,mp_sig) 

# # Z function
# def Z_function_internal_erosion(X):
#     k_dist,hwaterlevel_dist,D_dist,kh_dist,d_dist,Lf_dist,hp_dist,ysat_dist,mu_dist,ic_dist,d70_dist,mp_dist=X
#     Z=calculation_internal_erosion(hw=hwaterlevel_dist,k=k_dist,D=D_dist,dc=d_dist,kh=kh_dist,Lf=Lf_dist,hp=hp_dist,ysat=ysat_dist,mu=mu_dist,ich=ic_dist,d70=d70_dist,mp=mp_dist,return_Probability=False,print_information=False)
#     return [Z]


# # In[30]:


# # Create the event we want to estimate the probability
# marginals=[k_dist,hwaterlevel_dist,D_dist,kh_dist,d_dist,Lf_dist,hp_dist,ysat_dist,mu_dist,ic_dist,d70_dist,mp_dist]
# description=['k','waterlevel','D','kh','d','Lf','hp','ysat','mu','ic','d70','mp']
# Z=ot.PythonFunction(len(marginals),1,Z_function_internal_erosion)
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
# optimAlgo.setMaximumEvaluationNumber(N)
# optimAlgo.setMaximumAbsoluteError(1.0e-10)
# optimAlgo.setMaximumRelativeError(1.0e-10)
# optimAlgo.setMaximumResidualError(1.0e-10)
# optimAlgo.setMaximumConstraintError(1.0e-10)


# # Monte Carlo

# # In[31]:


# experiment = ot.MonteCarloExperiment()
# algo = ot.ProbabilitySimulationAlgorithm(event, experiment)
# algo.setMaximumCoefficientOfVariation(0.05)
# algo.setMaximumOuterSampling(int(N))
# algo.run()
# result=algo.getResult()

# pf_erosion_mc=result.getProbabilityEstimate()
# print ('The failure probability is: 1/',round(1/pf_erosion_mc), 'per year, which means the dike does not meet the required safety level of : 1/',round(1/Ptarget_erosion),'per year')


# # Form

# # In[32]:


# startvalues=distribution.getMean()
# startvalues[1]=10
# algo = ot.FORM(optimAlgo, event, startvalues)
# algo.run()
# result = algo.getResult()
# pf_erosion_form=result.getEventProbability()

# print ('The failure probability is: 1/',round(1/pf_erosion_form), 'per year, which means the dike does not meet the required safety level of : 1/',round(1/Ptarget_erosion),'per year')




# In[36]:


def calculation_stability (FoS,yd=1.06):
    beta=-(FoS/yd-0.41)/0.15
    return st.norm.cdf(beta,loc=0,scale=1)

P_stability=calculation_stability (FoS=0.97)
print ('The lowest safety factor of the figures is directive. Filling in the calibration formula a failure probability of:',round(1/P_stability),'per year is found. This does not meet the requirement of 1/',round(1/Ptarget_slope),'per year')



# In[38]:

dRc=0
q_tot=0

for j in range(100):
    dRc+=j/100
    q_tot=0
    for i in range(len(wind)):
        Direction=wind.iloc[:,0][i]
        F=wind.iloc[:,2][i]
        Uwind=wind.iloc[:,3][i]
        pf=wind.iloc[:,4][i]
        q_directionRL,blank=overtopping(Direction,Uwind,F,Rc+dRc)
        q_tot+=q_directionRL*pf/100
        
    #print ('\nThe total overtopping discharge is',round(q_tot,0),'l/m/s, so the dike does not meet the required safety level corresponding with a overtopping discharge of', round(qc_char,2), 'l/s/m')
    if q_tot<qc_char:
        print ('The crest level should be: ', round(hcrest+dRc,3), 'm to meet the criterion of the critical overtopping discharge of', round(qc_char,2), 'l/s/m, because the overtopping discharge is', round(q_tot,1),'l/s/m')
        break
        
print ('Possible motivation:\nNo it does not fit in the local conditions, because there is a road present on the crest of the dike. Increasing the dike height involves high costs. Therefore another solution is probably better.')


# In[39]:


Lf_newmean=Lf_mean

for j in range(100):
    Lf_newmean+=j/100
    Lf_new= characteristic_lognormal(Lf_newmean,Lf_newmean*0.1,5)#characteristic value
    q_tot=0
    Pf_new=calculation_internal_erosion(hw=hwaterlevel_10000,k=k_char,D=D_char,dc=d_char,kh=kh_char,Lf=Lf_new,hp=hp_char,ysat=ysat_char,mu=mu_char,ich=ic_char,d70=d70_char,mp=mp_char,return_Z=False,print_information=False)
    #print ('\nThe total overtopping discharge is',round(q_tot,0),'l/m/s, so the dike does not meet the required safety level corresponding with a overtopping discharge of', round(qc_char,2), 'l/s/m')
    if Pf_new<Ptarget_erosion:
        print ('The berm width should be: ', round(Lf_newmean-Lf_mean,1), 'm to meet the target failure probability of 1/', int(1/Ptarget_erosion), 'per year, because the failure probability for internal erosion is is 1/', int(1/Pf_new),'per year')
        break


