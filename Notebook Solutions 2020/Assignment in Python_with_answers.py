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
get_ipython().run_line_magic('matplotlib', 'inline')


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


N=10000 #sample size
x= ot.Sample(N,1)  #uses a dataset

for i in range(N):
    x[i,0]=CDF_to_PDF_to_distribution(waterlevel_CDF,N=1,delta=0.01,plot_cdf=False,plot_pdf=False,legend=False)[0]
    
x_dis=ot.UserDefined(x,[1/N]*N) # distribution based on dataset
a_dis=ot.Normal(20,10) # normal distribution with mean=20, standarddeviation=10
b_dis=ot.LogNormalMuSigma(10,300).getDistribution() # lognormal distributionm with mean=10 and standarddeviation=300

# Z function
def example_Z_function(X):
    x,a,b=X
    Z=a*x**2-b
    return [Z]

# Create the event you want to estimate the probability
marginals=[x_dis,a_dis,b_dis]
description=['x','a','b']
Z=ot.PythonFunction(len(marginals),1,example_Z_function)
RS = ot.CorrelationMatrix(len(marginals))
R = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(RS)
copula = ot.NormalCopula(R)
distribution = ot.ComposedDistribution(marginals, copula)
distribution.setDescription(description)
vect = ot.RandomVector(distribution)
G = ot.CompositeRandomVector(Z, vect)
event = ot.ThresholdEvent(G, ot.Less(), 0.0)

#Define a solver
optimAlgo = ot.Cobyla()
optimAlgo.setMaximumEvaluationNumber(100000)
optimAlgo.setMaximumAbsoluteError(1.0e-10)
optimAlgo.setMaximumRelativeError(1.0e-10)
optimAlgo.setMaximumResidualError(1.0e-10)
optimAlgo.setMaximumConstraintError(1.0e-10)
   
#run Form
startvalues=distribution.getMean()
startvalues[0]=7 #Let the waterlevel start at for example 7
algo = ot.FORM(optimAlgo, event, distribution.getMean())
algo.run()
result = algo.getResult()
print ('Form probability:',round(result.getEventProbability(),2))

#run Monte Carlo Simulation
experiment = ot.MonteCarloExperiment()
algo = ot.ProbabilitySimulationAlgorithm(event, experiment)
algo.setMaximumCoefficientOfVariation(0.05)
algo.setMaximumOuterSampling(int(N))
algo.run()
result=algo.getResult()

print ('Montecarlo probability',round(result.getProbabilityEstimate(),2 ))


# ---
# ## 2	Exercise description
# ### 2.1	Introduction
# In this exercise you are tasked to carry out a safety assessment and, subsequently, design of reinforcement measures for the failure mechanisms overtopping, piping and slope stability for a levee section at Water Authority Rivierenland.
# 
# Based on Dutch legislation, the safety standards of primary flood defences in the Netherlands must be assessed every 12 years. When the safety standard does not meet the requirements, reinforcement must be prepared and implemented. The rules, standards and methods for flood defence assessments and designs are noted in the `Wettelijk beoordelingsinstrumentarium 2017 (WBI2017)` and `Ontwerpinstrumentarium (OI2014v4)`. 

# ### 2.2	Scope
# The levee trajectory of interest is situated along the Lek from Amsterdam-Rijnkanaal to Fort Everdingen with a total length of 15.86 kilometres as shown in Figure 1. You are going to assess a 280 m section of the levee: section BF095.50 – BF098.30. This levee section is part of levee trajectory (dijktraject) 43-1. 
# 
# <img src="Figure1.png"  width="500"/>
# 
#  
# <center>  <i> Figure 1: Location of dike stretch 43-1 </i> </center>
# 
# The investigated levee section is located North West of Culemborg. Figure 2 shows an aerial photo of this section. Figure 3 shows a Street View photo, where you can see the distance between the levee and the houses. The representative cross-sections of this section are presented in Figure 4. These sections are derived from the elevation map, which is shown in Figure 5. 
#  
# <img src="Figure2.png"  width="500"/>
# 
# <center>  <i> Figure 2: Investigated part of the dike; dike section BF095.50 – BF098.30 </i> </center>
# 
# <img src="Figure3.png"  width="500"/>
#  
# <center>  <i> Figure 3: Google Maps Street view picture taken from the crest of the dike of interest ({Hyperlink: "www.google.nl/maps"}, October 2015) </i> </center>

# In[7]:


plot_cross_sections ('BF098_1',create_fig=True)
plot_cross_sections ('BF098_2')
plot_cross_sections ('BF097')
plot_cross_sections ('BF096', create_legend=True)


# <center>  <i> Figure 4: Different cross-sections of section BF095.50 – BF098.30. The river is located at the right side of the levee. Please note that this figure only shows the surface elevation. Depth of water bodies, like ditches or gullies are not displayed.</i> </center> 
# 
# <img src="Figure5.png"  width="500"/>
# 
# <center>  <i> Figure 5: Elevation map based on AHN3 (https://ahn.arcgisonline.nl/ahnviewer/). Please note that the depth of water bodies, like ditches or gullies are not displayed in this figure </i> </center>
# 

# ### 2.3 Safety standard
# An acceptable risk level can be determined based on the consequences of a flood at this location. In the event of a levee breach along trajectory 43-1, part of the Betuwe will flood. Based on Dutch legislation, the safety standard of a levee trajectory is based on a basic safety level for each individual person. This safety standard can be more strict based on a cost benefit analysis or societal risk.
#  
# Each individual in the Netherlands has the right to live in a place where flood safety meets at least the basic safety level. Additionally, this safety level can be increased when the area is of importance for economy, cultural heritage or when large groups of people can be exposed.
# 
# For the basic safety level, the risk of loss of life is expressed as: the local individual risk (LIR). The LIR describes the annual probability for a person to die in the area as a result of a flood. In the Netherlands an upper limit of this risk level has been set, which is called basic safety. For the Netherlands this basic safety level has been set at 1/100.000 per year. The LIR must, by law, at least be lower.
# 
# Given the consequences of a flood, the mortality and the degree of evacuation; the flood probability for this flood defence can be determined in order to meet the basic safety level of 1/100.000 per year.  
# 
# A consequence of a breach for trajectory 43-1 is shown in Figure 6. 
# 
# <img src="Figure6.jpg"  width="500"/> 
#  
# <center>  <i> Figure 6: Consequence dike breach for trajectory 43-1 </i> </center>
# 
# The calculation of the LIR is a combination of the probability of a flood, the mortality (probability of loss of life given the flood characteristics) and the degree of evacuation. 
# 
# $ LIR= P_{flood}×Mortality×(1-degree\_of\_evacuation) $
# 
# To ensure that the entire area meets the basic safety level (LIR of less than 1/100.000 per year); the most vulnerable part of this area must be considered. The Figure above shows an area north of the regional river Linge. When a levee breach occurs, the average water depth in this area is approximately 4.5 meters. 
# 
# Based on this water depth, the mortality can be calculated. 
# 
# Mortality: $F_D (h)= Փ_N (\frac{ln⁡(h)- μ_N}{σ_N});$
# $μ_N=7.60$ and
# $σ_N=2.75 $
# 
# The degree of evacuation for this area has been set at 55 percent. This means that in the event of a possible flood, 55 percent of the inhabitants can be evacuated preventively.
# 
# At this location other levee trajectories might flood the same area as well. Additionally, given this, and an executed cost-benefit analysis; the safety standard for trajectory 43-1 is set to a failure probability of $P = 1/10.000$ per year. 
# 

# ### 2.4	Load and strength
# In the exercise the three most important failure mechanisms for this location are elaborated on: overtopping and overflow, internal erosion and stability. For these failure mechanisms you are asked to carry out a safety assessment and a dike design; this chapter gives the input information. Information on the load and strength parameters must be used in this exercise. Read this chapter carefully.  
# 

# #### 2.4.1 Critical overtopping discharge
# 
# In a semi-probabilistic calculation for overtopping and overflow, the overtopping discharge has to be calculated with the loads (water level and wind) corresponding with the safety level of the dike cross-section for the mechanism. This overtopping discharge must be compared with the maximum allowed overtopping discharge. When the calculated discharge is lower than the maximum allowed overtopping discharge, the levee meets the required safety standard. The maximum allowed overtopping and overflow discharge depends on the wave height and revetment quality. The critical discharge is lognormal distributed with the parameters given in Table 4. 
# 
# ---
# 
# <center>  <i> Table 1: Mean and standard deviation for the critical discharge for different wave height and divot qualities </i> </center>
# 
# ||Closed  |divot|Open |divot|
# |------|------|------|------|------|
# |Wave height |μ [l/s/m]|	σ [l/s/m]|	μ [l/s/m]|	σ [l/s/m]|
# |0 – 1 m| 	225|	250	|100	|120|
# |1 – 2 m|	100|	120	|70	|80|
# |2 – 3 m |	70	|80	|40	|50|
# 

# #### 2.4.2 Hydraulic load parameters
# Hydraulic load is the driving force which can lead to levee failure. The exceedance frequency of the water levels is therefore very important. The exceedance frequency line of the water levels for the levee section that must be assessed is given in Figure 7 and came from Hydra-NL. Hydra-NL is a probabilistic model and calculates the hydraulic load statistics (water level, waves and overtopping) and integrates uncertainty about the median estimate. The water level corresponding with an exceedance probability of 1/10,000 per year is 7.35 m + NAP. And for an exceedance probability of 1/41.667 per year the water level is 7.56 m + NAP. The cumulative density function of the water level corresponding with the frequency line is attached in the file `Cumulative_density_function_water_level.xpt`. 
# 
# <img src="Figure7.png"  width="500"/>
# <center>  <i> Figure 7: Frequency line of the river water level</i> </center>
# 
# In addition to the water level, the wind is an important load parameter. The wind characteristics are divided into 16 wind directions. The wind velocity and fetch data per direction with an exceedance probability of 1/ 41,667 per year are shown in Table 1. The exceedance probability per wind velocity per direction is given in the file `Exceedence_probability_wind_velocity.xlsx`. 
# 
# ---
#  
# <center>  <i> Table 2: Wind and fetch data with an exceedance probability of 1/41667 per year</i> </center>

# In[8]:


wind


# The wave height and period can be calculated from the provided water depths and wind data. Reference is made to the Hydraulic `Structures Manual (CTB3355)`, which can be found on Brightspace under Readings. Please use the `Young and Verhagen formula`. Formulas concerning the mean and maximum overtopping discharge can be found in the lecture notes of this course. When the water level is higher than the crest of the dike, the part of overflowing water can be calculated by:  
# 
# $q_{overflow}=0.54 \sqrt{g|R_c^3 |}$        
# $R_c<0$	
# 
# In which: <br/>
# $R_c$ = Freeboard $(h_{crest}-h)$

# #### 2.4.3 Geotechnical parameters
# 
# <img src="Figure8.png"  width="200"/>
# 
# 
# <center>  <i> Figure 8: Soil boring at toe of the dike with depth in m +NAP (left overview and right detail) </i> </center>
# 
# <center>  <i> Table 3: Mean geotechnical parameters of the soil boring </i> </center>
# 
# |Colour|	Material|	Saturated volumetric weight [kN/m³]|	Dry Volumetric weight[kN/m³]|	Cohesion [kN/m²]|	Friction angle [deg]	|Undrained shear strength ratio [-]|	Strength increase exponent [-]	|POP top [kN/m²]|	POP bottom [kN/m²]|
# |:------|------|------|------|------|------|------|------|------|------|
# | <font color='purple'>purple</font>|Clay (dike material)	|18.00 	|18.00	|0.01|	32	|0.31|	0.90	|30	|30|
# |<font color='green'>green</font>|	Peat|	14.50|	14.50|	0.01|	30|	0.29|	0.90|	24|	24|
# |<font color='pink'>pink</font>|	Clay	|15.00	|15.00|	0.01|	30	|0.29|	0.90|	24	|24|
# |<font color='orange'>orange</font>|	Sand	|18.00	|20.00|	0	|34	|-|	-|	-	|-|

# #### 2.4.4 Parameters and their distributions
# Some geometry and geotechnical parameters and their distributions are provided in Table 4, make reasonable assumption for the remaining ones. For the roughness coefficient assume the entire revetment is covered with a grass layer. The geotechnical length profile is given in the files ‘Geotechnical_length_profile.tif’.  The dike consists completely out of clay and the result of the representative soil boring is shown in Figure 8 and Table 3. The boring is taken close to the levee. More detailed information about the geotechnical parameters are shown in ‘Geotechnical_data.xlsx’. 
# 
# ---
# <center>  <i> Table 4: Geometry and geotechnical parameters (*these parameters need to be estimated from the available graphs and pictures, or from other sources (DINO database, Google Earth, etc.)) </i> </center>
# 
# |    Parameters    |    Distribution    |    Mean    |    Standard   deviation  or Coefficient of variation     |
# |:------|:------|:------|:------|
# |70%-fractile of   grain size distribution [m]    |    Lognormal    |    2.8e-4    |  $V=0.12$      |
# |Angle normal to dike [deg]    |    Deterministic    |    12    |    -    |
# |    Aquifer thickness   [m]    |    Lognormal    |    25     |    $\sigma=0.5$       |
# |    Bedding angle [deg]    |    Deterministic    |    37    |    -    |
# |    Bottom level [m +   NAP]    |    Normal    |    2.35    |   $\sigma=0.3$      |
# |    Cohesion [kN/m²]     |    Lognormal    |    See table 3    |  $V=0.0001$      |
# |    Constant of White   [-]    |    Deterministic    |    0.25    |    -    |
# |    Critical heave   gradient [-]    |    Lognormal    |    0.5    |  $\sigma=0.1$       |
# |    Dike height [m]    |    Deterministic    |    Estimate    |    -    |
# |    Fetch length [m]    |    Deterministic    |    See table 2   |    -    |
# |    Friciton Angle [deg]    |    Lognormal    |    See table 3    | $V=0.000333$ & $ V =0.1$ (for dike material) |
# |    Gravitational   constant [m/s²]    |    Deterministic    |    9.81    |    -    |
# |    Hinterland phreatic   level [m + NAP]    |    Normal    |    3.5     | $\sigma=0.1$        |
# |    Hydraulic   conductivity aquifer [m/s]    |    Lognormal    |    7.52e-4    | $V=0.50$      |
# |    Hydraulic   conductivity aquitard [m/s]    |    Lognormal    |    1.00e-6    |    $V=0.50$    |
# |Intrusion length [m]| Deterministic|1.0|-|
# |    Kinematic viscosity   [m²/s]    |    Deterministic    |        |    -    |
# |    Length (effective)   foreshore [m]    |    Lognormal    |    20     |     $V=0.1$   |
# |    Model factor piping   [-]    |    Normal    |    1.0    | $\sigma=0.12$        |
# |    Model factor uplift   [-]    |    Normal    |    1.0    |   $\sigma=0.1$      |
# |    POP    |    Lognormal    |         |     $V=0.30$ & $ V =0.45$ (for dike material)      |
# |    Reference value of   70%-fractile of grain size distribution [m]    |    Deterministic    |    2.08e-4    |    -    |
# |    Saturated   volumetric weight blanket [kN/m³]    |    Normal    |    Estimate      |    $V=0.05$    |
# |    Slope     |    Deterministic    |    Estimate   |    -    |
# |    Slope width [m]    |    Deterministic    |    Estimate   |    -    |
# |    Strength increase component [-]    |    Lognormal    |    See table 4    |           $V=0.033$ |
# |    Thickness   hinterland blanket [m]    |    Lognormal    |    Estimate   |       $\sigma=0.5$  |
# |    Undrained shear strength ratio    |    Lognormal    |    See table 3    |         $V=0.207$ & $ V =0.199$ (for dike material)  |
# |    Volumetric weight   sand grains [kN/m³]    |    Deterministic    |    16.5    |    -    |
# |    Volumetric weight   water [kN/m³    |    Deterministic    |    10    |    -    |
# |    Width levee [m]    |    Deterministic    |    Estimate     |    -    |
# 
# 
# For a semi-probabilistic assessment, one should keep the 5% and 95% conﬁdence intervals in mind. In the safety assessment the confidence intervals for a lognormal distribution are calculated with the formulas given below:
# 
# $ X_{lower} =\exp{⁡( μ_m-t_{n-1} σ_m  \sqrt{(1-a)+1/n)})}$   and    $X_{upper}=\exp{⁡(μ_m+t_{n-1} σ_m  \sqrt{((1-a)+1/n))}}$
# 
# $σ_m^2=\ln⁡{[1+(σ_X/μ_X )^2]}$
# $μ_m=\ln{[μ_X] }-1/2 σ_m^2$
# 
# In which:<br/>
# $X$ = Confidence bounds<br/>
# $μ_X$=	Mean value <br/>
# $μ_m$=	Mean value of a logarithm <br/>
# $t_{n-1}$=	Student t-factor corresponding with the coincidence bound (taket_(n-1)=1.76)<br/>
# $σ_X$=	Standard deviation<br/>
# $σ_m$=	Standard deviation of a logarithm <br/>
# $a$=	Ratio between the local and regional variation (takea=0)<br/>
# $n$=	Number of observations (taken=15) <br/>
# 

# #### 2.4.5	Slope stability
# The failure mechanism slope stability is analysed qualitatively later on this exercise. In the semi-probabilistic safety assessment, the water level corresponding with the safety level (1/10,000 per year) has to be used. The safety factor for slope stability can be calculated with the software D-Stability. Some results of the calculation are shown in the figures below. Note that the same stability calculation is discussed in class. To see the influence of a design choice on the safety factor, you are free to use the discussed exercise. 
# 
#  <img src="Figure9.png"  width="500"/>
# <center>  <i> Figure 9: Safety factor for inner slope stability = 0.97 </i> </center>
# 
#   <img src="Figure10.png"  width="500"/>
# <center>  <i> Figure 10: Safety factor for inner slope stability = 1.06 </i> </center>
#  
#  <img src="Figure11.png"  width="500"/>
# <center>  <i> Figure 11: Safety factor for inner slope stability = 1.10 </i> </center>
# 
#   <img src="Figure12.png"  width="500"/>
# <center>  <i> Figure 12: Safety factor for inner slope stability = 1.28 </i> </center>

# ### 2.5 Calibration formula’s
# With the use of calibration formulas safety factors can be transformed into failure probabilities for the failure mechanisms uplift, heave, piping and slope stability. The calibration formulas are given below:
# 
# Calibration formula uplifting:
# $P_{f;u}=ϕ[-\frac{\ln⁡{(F_u/0.48)}+0.27β_{req})}{0.46}]$
# 
# Calibration formula heave:
# $P_{f;h}=ϕ[-\frac{\ln⁡{(F_h/0.37)}+0.3 β_{req})}{0.48}]$
# 
# Calibration formula piping:
# $P_{f;p}=ϕ[-\frac{\ln⁡{(F_p/1.04)}+0.43 β_{req})}{0.37}]$
# 
# Calibration formula slope stability:
# $P_{f;i}=ϕ[-\frac{(F_{d,i}/γ_d )-0.41}{0.15}]$
# 
# In which:<br/>
# $P_{f;j}$ =	The failure probability of the failure mechanism j [1/year]<br/>
# $F$ =	Safety factor of the failure mechanism [-]<br/>
# $γ_d$ =	Model factor, take γ_d=1.06  [-]<br/>
# $β_{req}$ =	The reliability index for the dike segment [-]<br/>

# ## 3 Questions
# ### 3.1 Safety assessment
# In this exercise you are tasked to carry out a safety assessment for the given dike trajectory.
# #### Question 1
# What are the relevant load conditions in this area?  

# In[9]:


print ('Section 3.3 CIE5314 Flood defence \nThe forcing characteristics in this area are the hydraulic loads. Because the dike trajectory is located in the upstream part of the river, the hydraulic loads are affected by the river discharge, while in the downstream part both the discharge and the tide affect the water level. The hydraulic loads are also affected by the wind. The wind causes waves, these wind waves are limited compared to the waves at sea because the fetch is limited by the river shape. ')


# #### Question 2
# Calculate the minimum required safety level for this area.

# In[10]:


mortality=1.19/100 # 1.19 percent
LIR=1/100000
degree_of_evacuation=0.55

Ptraject=LIR/(mortality*(1-degree_of_evacuation))

print ('The minimum required safety level for this area is 1/',round (1/Ptraject),'per year')


# #### Question 3	
# Considering the area of interest and the determination of the safety level of 1/10.000 per year. What can be an important factor for this high safety standard?

# In[11]:


print( 'Many people are living in this area. When a flood occurs, large groups of people will be exposed. Therefore the minimum required safety level is stricter than the calculated basic safety level. ')


# #### Question 4
# Calculate the maximum allowed probability of failure for the failure mechanism slope instability, piping and overtopping.
# 
# ---

# **Answer:** <br/>
# `Section 10.2.2 CIE5314 Flood defence`<br/>
# The failure probability per failure mechanism is calculated with the formula given below.<br/>
# <br/>
# $P_{(req|mech)}=ω_j P_{req}$<br/>
# <br/>
# With:<br/>
# $P_{(req|mech)}$ =Required probability of failure per failure mechanism [1/year]<br/>
# $P_{req}$ = Safety standard of the dike stretch [1/year]<br/>
# $ω$ =Contribution of a failure mechanism tot the total failure

# In[12]:


#input
Preq=1/10000
w_overtopping=0.24
w_erosion=0.24
w_slope=0.04

#calculation
P_overtopping=w_overtopping*Preq
P_erosion=w_erosion*Preq
P_slope=w_slope*Preq

print ('This results in the following maximum allowed probability of failure: ')
print ('Overtopping and overflow: 1/',round(1/P_overtopping),'per year')
print ('Internal erosion: 1/',round(1/P_erosion),'per year')
print ('Slope stability: 1/',round(1/P_slope),'per year')


# #### Question 5
# Calculate the target failure probability for a levee section for the failure mechanisms slope instability, piping and overflow. 
# 
# ---
# 
# **Answer** <br/>
# The target failure probability is calculated with the formula given below <br/>
# <br/>
# $P_{req,i,j}=P_{req|mode}/N_j = ω_j/N_j P_{req} $  <br/>
# $N_j=1+(a_j L)/b_j $<br/>
# <br/>
# With:<br/>
# $N_j$	Length effect factor [-]<br/>
# $a_j$	Faction of the trajectory’s length that is sensitive to failure mechanism j [-]<br/>
# $b_j$	Length of a typical independent section for failure mechanism j [m]<br/>
# $P_{req,i,j}$	Required annual failure probability for section i and failure mechanism j [-]<br/>
# $L$		Length of the trajectory [m]<br/>
# 

# In[13]:


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

print ('This results in the following maximum allowed probability of failure per dike cross section:')
print ('Overtopping and overflow: 1/',round(1/Ptarget_overtopping),'per year')
print ('Internal erosion: 1/',round(1/Ptarget_erosion),'per year')
print ('Slope stability: 1/',round(1/Ptarget_slope),'per year')


# #### Question 6
# Comment on the different target failure probabilities per mechanisms based on the distribution over the mechanisms and the length-effect.  

# In[14]:


print ('Section 10.2.3 CIE5314 Flood defence \nThe failure probability also increases with the length of the dike trajectory, caused by uncertainties and variations in the hydraulic loads, geotechnical resistance, etc. \nThe length-effect for wave overtopping is caused by variation of wind direction and orientation of the dike section, cross section profile and types of revetments. The length effect for internal erosion and slope stability can be caused by spatial variation of the soil parameters that are relevant for the resistance (density, permeability, etc.)\nThe length effect is small for failure mechanisms with a high dependency between the cross-sections. This is often the case for overtopping which is characterized by a large auto-correlation in space because dike elevation and orientation along a trajectory will (often) be similar. The other two failure mechanism have cross-section that are little dependent on each other. The additive effect of individual probabilities of the section on the segment reliability is stronger than for overtopping. Therefore the length effect factor for these mechanisms is much higher. ')


# #### Question 7
# Which levee cross section is representative for the failure mechanism overtopping an overflow? Motivate your answer.

# In[15]:


print ('Cross section BF098_1, because the crest height is lowest and the lower the height the larger the amount of overtopping. ')


# #### Question 8 
# Motivate which critical discharge you use in the overtopping and overflow calculation?
# 
# ----
# 
# **Answer**
# 
# For example: 
# 
#     • 100 l/s/m. The waves are lower than 1 meter, and I don’t know the revetment quality. Therefore I choose conservative for the lowest value. 
#     •225 l/s/m. The waves are lower than 1 meter, and on Google Street view I see that the revetment quality is in good condition. 
# 

# In[16]:


qc_mean= 100 #l/s/m 


# #### Question 9
# Calculate the characteristic values with 5% or 95% exceedance probability to be used in the overtopping and overflow calculation.  

# In[17]:


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
for i in range(len(wind)):
    Direction=wind.iloc[:,0][i]
    F=wind.iloc[:,2][i]
    Uwind=wind.iloc[:,3][i]
    pf=wind.iloc[:,4][i]
    q_direction=Calculation_overflow(Direction,F,Uwind)
    print ('The overtopping discharge in wind direction:', Direction,'degrees with fetch length:',F,'m and wind velocity:', Uwind,'m/s\nis equal to',round(q_direction,0),'l/m/s')
    q_tot+=q_direction*pf/100

print ('\nThe total overtopping discharge is',round(q_tot,1),'l/m/s, so the dike does not meet the required safety level corresponding with a overtopping discharge of', round(qc_char,1), 'l/s/m')


# #### Question 11
# Perform a fully probabilistic Monte Carlo calculation. Use a full distribution for the wind speed of one wind direction. You can use your own code, the `OpenTurns`  module of Python (see paragraph 1.2 for an example) or the `Prob2B` software, provided on Blackboard. You can install the OpenTurns model running the next line in the command prompt: <br/>
# `pip install openturns`
# 
# To transform the given CDF of the waterlevel to a distribution you may use the given function specified in the begin of the exercise:
# `CDF_to_PDF_to_distribution`
# 
# ---
# 
# **Answer** <br/>
# Possible motivations for wind direction are:
# 
#     •	WNW: This wind direction has the largest contribution (largest pf×q)
#     •	ENE: This wind direction is the most common (largest pf)
#     •	NE: Overtopping discharge closest to the total overtopping discharge
# 	
# If a Weibull distribution is used the following parameters have to be found (a fit of the tail is most important):

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


#Distribution
N=100000
N_array=[1/N]*N
hwaterlevel_dis= ot.Sample(N,1)
Uwind_dis= ot.Sample(N,1)
for i in range(N):
    hwaterlevel_dis[i,0]= CDF_to_PDF_to_distribution(waterlevel_CDF,N=1,delta=0.01,plot_cdf=False,plot_pdf=False,legend=False)[0]
    Uwind_dis[i,0]= weibull_min.ppf(np.random.random(1),fitpar[0][0],0,fitpar[0][1])[0]
h_bottom_dist=ot.Normal(d_mean,d_sigma)
hwaterlevel_dist=ot.UserDefined(hwaterlevel_dis,N_array)
Uwind_dist=ot.UserDefined(Uwind_dis,N_array)
qcrit_dist= ot.LogNormalMuSigma(qc_mean,qc_sigma).getDistribution()


# In[22]:


#Z-function
def Z_overflow(X):
    h_bottom,h_waterlevel,Uwind,qcrit=X
    Z=qcrit-Calculation_overflow(Direction,F,Uwind,hbottom=h_bottom,hwaterlevel=h_waterlevel,g=g,hcrest=hcrest,yb=yb,yf=yf,yv=yv,anglenormal=anglenormal,slope=slope)
    return [Z]

# Create the event we want to estimate the probability
marginals=[h_bottom_dist,hwaterlevel_dist,Uwind_dist,qcrit_dist]
description=['Bottom height','Waterlevel','Wind veolocity', 'Critical overtopping discharge']
Z=ot.PythonFunction(len(marginals),1,Z_overflow)
RS = ot.CorrelationMatrix(len(marginals))
R = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(RS)
copula = ot.NormalCopula(R)
distribution = ot.ComposedDistribution(marginals, copula)
distribution.setDescription(description)
vect = ot.RandomVector(distribution)
G = ot.CompositeRandomVector(Z, vect)
event = ot.ThresholdEvent(G, ot.Less(), 0.0)

#Define a solver
optimAlgo = ot.Cobyla()
optimAlgo.setMaximumEvaluationNumber(100000)
optimAlgo.setMaximumAbsoluteError(1.0e-10)
optimAlgo.setMaximumRelativeError(1.0e-10)
optimAlgo.setMaximumResidualError(1.0e-10)
optimAlgo.setMaximumConstraintError(1.0e-10)
        
#run MC
experiment = ot.MonteCarloExperiment()
algo = ot.ProbabilitySimulationAlgorithm(event, experiment)
algo.setMaximumCoefficientOfVariation(0.05)
algo.setMaximumOuterSampling(int(1e6))
algo.run()
result=algo.getResult() 
pf_overtopping=result.getProbabilityEstimate()

print ('The failure probability for the wind direction '+winddirection,'is: 1/',round(1/pf_overtopping), 'per year, which means the dike meets the required safety level of : 1/',round(1/Ptarget_overtopping),'per year')


# #### Question 12
# Is there a difference in the safety judgment between the semi-probabilistic and probabilistic calculation? If so, what could be the reason?

# In[23]:


print ('Yes, in the semi-probabilistic calculation you use characteristic values, which is unfavourable conservative combination between the variables. In a Monte Carlo analysis the distributions of the variables are used, instead of the most conservative value of that parameter. ')


# #### Question 13
# Which levee cross section is representative for the failure mechanism internal erosion? Motivate your answer.

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


# #### Question 16
# What is the dominant failure mechanism for internal erosion and what is the cause of the dominant failure mechanism?

# In[27]:


print ('The dominant failure mechanism is piping, because the leakage length is relatively long caused by the impermeable foreland. ')


# #### Question 17
# Do you expect to find larger or lower failure probabilities with a probabilistic calculation?

# In[28]:


print('I expect to find lower failure probabilities with a probabilistic calculation because in the semi-probabilistic calculation you make use of characteristic values, which is the most unfavourable combination between the variables. In a probabilistic analysis the distributions of the variables are used, instead of the most unfavourable value of that parameter. ')


# #### Question 18
# Perform a fully probabilistic calculation Monte Carlo and Form calculation. An example is given in section 1.2
# 
# ---
# **Answer** <br/>
# Distributions

# In[29]:


N=100000
N_array=[1/N]*N
hwaterlevel_dis= ot.Sample(N,1)

for i in range(N):
    hwaterlevel_dis[i,0]=CDF_to_PDF_to_distribution(waterlevel_CDF,1,delta=0.01,plot_cdf=False,plot_pdf=False,legend=False)[0]
    
k_dist=ot.Normal(k_mean,k_sig)
hwaterlevel_dist=ot.UserDefined(hwaterlevel_dis,N_array)
D_dist=ot.LogNormalMuSigma(D_mean,D_sig).getDistribution()
kh_dist=ot.Normal(kh_mean,kh_sig)
d_dist=ot.LogNormalMuSigma(d_mean,d_sig).getDistribution()
Lf_dist=ot.LogNormalMuSigma(Lf_mean,Lf_sig).getDistribution()
hp_dist=ot.Normal(hp_mean,hp_sig)
ysat_dist=ot.Normal(ysat_mean,ysat_sig)
mu_dist=ot.Normal(mu_mean,mu_sig)
ic_dist=ot.LogNormalMuSigma(ic_mean,ic_sig).getDistribution()
d70_dist=ot.LogNormalMuSigma(d70_mean,d70_sig).getDistribution()
mp_dist=ot.Normal(mp_mean,mp_sig) 

# Z function
def Z_function_internal_erosion(X):
    k_dist,hwaterlevel_dist,D_dist,kh_dist,d_dist,Lf_dist,hp_dist,ysat_dist,mu_dist,ic_dist,d70_dist,mp_dist=X
    Z=calculation_internal_erosion(hw=hwaterlevel_dist,k=k_dist,D=D_dist,dc=d_dist,kh=kh_dist,Lf=Lf_dist,hp=hp_dist,ysat=ysat_dist,mu=mu_dist,ich=ic_dist,d70=d70_dist,mp=mp_dist,return_Probability=False,print_information=False)
    return [Z]


# In[30]:


# Create the event we want to estimate the probability
marginals=[k_dist,hwaterlevel_dist,D_dist,kh_dist,d_dist,Lf_dist,hp_dist,ysat_dist,mu_dist,ic_dist,d70_dist,mp_dist]
description=['k','waterlevel','D','kh','d','Lf','hp','ysat','mu','ic','d70','mp']
Z=ot.PythonFunction(len(marginals),1,Z_function_internal_erosion)
RS = ot.CorrelationMatrix(len(marginals))
R = ot.NormalCopula.GetCorrelationFromSpearmanCorrelation(RS)
copula = ot.NormalCopula(R)
distribution = ot.ComposedDistribution(marginals, copula)
distribution.setDescription(description)
vect = ot.RandomVector(distribution)
G = ot.CompositeRandomVector(Z, vect)
event = ot.ThresholdEvent(G, ot.Less(), 0.0)

#Define a solver
optimAlgo = ot.Cobyla()
optimAlgo.setMaximumEvaluationNumber(N)
optimAlgo.setMaximumAbsoluteError(1.0e-10)
optimAlgo.setMaximumRelativeError(1.0e-10)
optimAlgo.setMaximumResidualError(1.0e-10)
optimAlgo.setMaximumConstraintError(1.0e-10)


# Monte Carlo

# In[31]:


experiment = ot.MonteCarloExperiment()
algo = ot.ProbabilitySimulationAlgorithm(event, experiment)
algo.setMaximumCoefficientOfVariation(0.05)
algo.setMaximumOuterSampling(int(N))
algo.run()
result=algo.getResult()

pf_erosion_mc=result.getProbabilityEstimate()
print ('The failure probability is: 1/',round(1/pf_erosion_mc), 'per year, which means the dike does not meet the required safety level of : 1/',round(1/Ptarget_erosion),'per year')


# Form

# In[32]:


startvalues=distribution.getMean()
startvalues[1]=10
algo = ot.FORM(optimAlgo, event, startvalues)
algo.run()
result = algo.getResult()
pf_erosion_form=result.getEventProbability()

print ('The failure probability is: 1/',round(1/pf_erosion_form), 'per year, which means the dike does not meet the required safety level of : 1/',round(1/Ptarget_erosion),'per year')


# #### Question 19
# Is there a difference in the safety judgment between the semi-probabilistic and probabilistic calculation? 

# In[33]:


print('No, both the semi-probabilistic as probabilistic calculation does not meet the required safety standard. ')


# #### Question 20
# Does the result match your expectations? Explain why: 

# In[34]:


print ('The result does not match with my expectations, because the failure probability became higher in the probabilistic calculation. The reason for this is the use of the calibration formula. The formula is calibrated in such a way 20% of the cases has a higher failure probability and 80% a lower failure probability. The calculated case belongs to the 20%. ')


# #### Question 21
# Which dike cross-section is representative for the failure mechanism slope stability? Motivate your answer.

# In[35]:


# Or: BF098_1, this cross-section is lowest causing a saturated dike.  
print('Cross  section BF097, because this dike cross-section is highest and has the steepest slope. ')


# #### Question 22
# Does the dike meet the required safety standard in a semi-probabilistic calculation? 
# 

# In[36]:


def calculation_stability (FoS,yd=1.06):
    beta=-(FoS/yd-0.41)/0.15
    return st.norm.cdf(beta,loc=0,scale=1)

P_stability=calculation_stability (FoS=0.97)
print ('The lowest safety factor of the figures is directive. Filling in the calibration formula a failure probability of:',round(1/P_stability),'per year is found. This does not meet the requirement of 1/',round(1/Ptarget_slope),'per year')


# ### 3.2 Design 
# The assessed dike does not fulfil the required safety standards, so a redesign of the levee must be made to make it fulfil the certain reliability requirement. First you are asked to give different design measures per failure mechanism after that a design have to be made that meets the required safety standard. 
# #### Question 23
# Do you see any optimizations possibilities without reinforcement. 

# In[37]:


print ('Reduction of uncertainties, this will reduce the characteristic values in the analysis. Also in the fully probabilistic analysis the standard deviations are smaller, leading to a reduction of failure probability.')


# #### Question 24
# What should the crest level need to be in order to meet the criterion of a critical overtopping discharge? Argue why this solution fits or does not fit in the local conditions. 

# In[38]:


hcrestnew=hcrest
q_tot=0

for j in range(100):
    hcrestnew+=j/100
    q_tot=0
    for i in range(len(wind)):
        Direction=wind.iloc[:,0][i]
        F=wind.iloc[:,2][i]
        Uwind=wind.iloc[:,3][i]
        pf=wind.iloc[:,4][i]
        q_direction=Calculation_overflow(Direction,F,Uwind,hcrest=hcrestnew)
        q_tot+=q_direction*pf/100
    #print ('\nThe total overtopping discharge is',round(q_tot,0),'l/m/s, so the dike does not meet the required safety level corresponding with a overtopping discharge of', round(qc_char,2), 'l/s/m')
    if q_tot<qc_char:
        print ('The crest level should be: ', round(hcrestnew,2), 'm to meet the criterion of the critical overtopping discharge of', round(qc_char,2), 'l/s/m, because the overtopping discharge is', round(q_tot,1),'l/s/m')
        break
        
print ('Possible motivation:\nNo it does not fit in the local conditions, because there is a road present on the crest of the dike. Increasing the dike height involves high costs. Therefore another solution is probably better.')


# ### Question 25
# Which design  measures could you propose for the failure mechanism overtopping and overflow? Motivate your  suggestions qualitatively and with sketches. Prioritize your alternatives and motivate your answer.
# 
# ---
# 
# **Answer** <br/>
# Chapter 5 CIE5314 Flood defence
# 
#     •	Stronger inner slope- higher critical discharge
#         o	Hard revetments (stone pitching)
#     •	Add a Berm- Wave height reduction 
#     •	Rougher slope- Reduction of run up
#     •	Less steep slope- Lower breaking parameter
#     •	Increase crest- Reduction of overtopping discharge

# In[ ]:





# #### Question 26
# 
# What should the berm width need to be in order to meet the safety standard? Argue why this solution fits or does not fit in the local conditions. 

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


# In[40]:


print ('Possible motivation:\nThis solution does not fit the local conditions, because there are houses very close to the levee.')


# #### Question 27
# What design measures could you propose for the failure mechanism internal erosion? Motivate your suggestions qualitatively and with sketches. Prioritize your alternatives and motivate your answer.
# 
# ---
# 
# **Answer** <br/>
# *Section 7.3 CIE5314  Flood defence*
# 
#     •	Increase seepage path
#         o	Increase dike width
#         o	Construct a berm
#         o	Relief wells
#         o	Seepage screens or cut-of walls
#     •	Reinforce blanket layer.
#     •	Drainage and filters
# 

# In[ ]:





# #### Question 28
# What design measures could you propose for the failure mechanism slope stability? Motivate your suggestions qualitatively and with sketches. Prioritize your alternatives and motivate your answer.
# 
# ---
# 
# **Answer** <br/>
# *Section 6.4 and 6.5 CIE5314 Flood defence*
# 
#     •	Increase weight on the passive side 
#     •	Reduce slope (shallow slopes especially on the land-side)
#     •	Construct a berm
#     •	Increase strength
#     •	Increase the dike width
#     •	Reduce permeability dike body material
#     •	Place structural elements (e.g. a sheet pile)
#     •	Construct drainage elements
# 

# In[ ]:





# #### Questin 29
# Work out one possible design qualitative. Take the local conditions (buildings, trees etc.) into account, in a way your design fits. Motivate your answer

# In[ ]:




