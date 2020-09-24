
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openturns as ot #run pip install openturns in the cmd kernal
from scipy.interpolate import interp1d
from scipy import stats as st
# get_ipython().run_line_magic('matplotlib', 'inline')
from tabulate import tabulate
from scipy import stats as st
from scipy import special
from scipy.stats import lognorm

#### INPUT DATA ####

### FROM EXERCISE ###

data_loc='<insert_drive_path_with_all_files_here>' #You have to change this to the data folder
cross_sections=pd.read_excel(data_loc+'cross_sections.xlsx')
wind=pd.read_excel(data_loc+'wind.xlsx') #0:direction in degree 1:fetch length,2: windspeed,3: occurrencprobability

wind_dir = wind[0]
wind_fetch_length = wind[1]
wind_speed = wind[2]
wind_occurrence = wind[3]

wind_exceedence=pd.read_excel(data_loc+'Exceedence_probability_wind_velocity.xlsx',index_col=0)
waterlevel_CDF=np.loadtxt(data_loc+'Cumulative_density_function_water_level.xpt')

## Q9 ##
# GENERAL #
h_waterlevel = 1
g = 9.81
anglenormal = 12 

# DIKE PROPERTIES #
xmin=222 #need to be visually determined
xmax=247 #need to be visually determined
h_crest=np.max(cross_sections['BF098_1, z']) #adjust name of representative crosssection
xdata=cross_sections['BF098_1, x'][xmin:xmax]
ydata=cross_sections['BF098_1, z'][xmin:xmax]
a = st.linregress(xdata,ydata)
slope=1/-a[0] 
alpha=np.arctan(1/slope)
yfit=a[0]*xdata+a[1]
#plt.plot(cross_sections['BF098_1, x'][np.argmax(cross_sections['BF098_1, z'])],hcrest,'o', label='crestheight')
#plt.plot(xdata,yfit,color='black',linestyle='--',label='Inner slope')
#plt.legend()
print('The slope is 1:',round(slope,1))
print ('The maximum height is',round(h_crest,2),'m')

h_berm = 0 # height of berm (if present)
b_berm = 0 # widht of berm (if present)
freeboard_berm = h_berm - h_waterlevel # (if present)
freeboard_dike = h_crest - h_waterlevel
Lberm = 1 # defined as distance between h_berm - Hm0    and min(h_berm + Hm0,h_crest) , see p.137 EuroTOP (2018)
    

yf=1 # correction factor for permeability and roughness of slope (if any)
yv=1 # correction factor for vertical wall on top (if any)

p_design_bedlevel = 0.05 # Load or Strength parameter (for deterministic either 95% or 5%)
mu_bedlevel = 1 # mean
sigma_bedlevel = 0.2 # std. dev. 
h_bedlevel = st.norm.ppf(p_design_bedlevel,mu_bedlevel,sigma_bedlevel)
d_water = h_waterlevel - h_bedlevel

# Critical discharge #
t_n1 = 1.76
n = 15

mu_q_crit01 = 100
sigma_q_crit01 = 120
sigma_m_q_crit01_squared = np.log(1 + (sigma_q_crit01 / mu_q_crit01)**2)
mu_m_q_crit01 = np.log(mu_q_crit01) - 0.5 * sigma_m_q_crit01_squared

q_crit01_lower= np.exp(mu_m_q_crit01 - t_n1 * np.sqrt(sigma_m_q_crit01_squared) * np.sqrt(1+1/n)) # qcrit is a strength

