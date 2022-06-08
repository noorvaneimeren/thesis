
"""
Created on Thu May 12 11:46:19 2022

@author: Noor van Eimeren
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd 
import plotly.express as px
import os 
import datetime as dt
from functools import reduce
import math

"""This script returns the value of the particle size passing through a particle counter
 To find the corresponding particle size, call the function:
  particlesize_func('Transformed flow rate [m3/s]', 'Potential [Volt]')
The transformed flow rate is defined as (1/Flow Rate [m3/s])
    x = flow rates
    y = potential
    z = particle size """
    
os.chdir('D:/Users/eimerno/Documents/Python Scripts/Interpolation')
#%% Import data 
"""Input datasets"""
# Calibration data (Note: not necessary to change by user)
df_pc_calibration_data = pd.read_excel('Interpolation_Data.xlsx')

# Potentials used for experiment at set flow rate 
# Note: this is only for one potential for simplicity 6
df_pc_particles = pd.read_excel('Book4.xlsx', skiprows = 15, parse_dates=True)

df_pc_particles.rename(columns={'Unnamed: 0' : 'Datetime'}, inplace=True) 
# Flowrates
df_pc_flowrates = pd.read_excel('Flowrates.xlsx')
df_pc_flowrates.columns = ['Datetime', 'Flowrate']

#pc_experiment_voltage = list(pc_experiment_data.columns)


#%% Create new final dataset with all the information and interpolated flowrates. 

df_pc_particles['Datetime'] = pd.to_datetime(df_pc_particles['Datetime'])
df_pc_particles['Datetime_float'] = (df_pc_particles['Datetime'].dt.strftime('%Y%m%d%H%M%S')).astype(float)

df_pc_flowrates['Datetime'] = pd.to_datetime(df_pc_flowrates['Datetime'])
df_pc_flowrates['Datetime_float'] = (df_pc_flowrates['Datetime'].dt.strftime('%Y%m%d%H%M%S')).astype(float)


# Merge  dataframes
dfs = [df_pc_particles, df_pc_flowrates]
df_pc_final = reduce(lambda left, right: pd.merge(left,right, on=['Datetime_float'], how='outer'), dfs)

flow_interpolation_time = (df_pc_flowrates ['Datetime_float']).astype(float).tolist()
flow_interpolation_flow = df_pc_flowrates ['Flowrate'].tolist()
time_interpolation_fun= interpolate.interp1d (flow_interpolation_time, flow_interpolation_flow, kind='linear')

df_pc_final['Interpolated Flowrates [s/m3]'] =  1/((time_interpolation_fun(df_pc_final['Datetime_float']) * 1.667e-8))

# Drop redundant columns
df_pc_final.drop(['Datetime_y', 'Flowrate'], axis=1, inplace=True)



#%% Create interpolation function for particle sizes
"""Create intetpolation function:
        Define data from provided excel sheet containing all known datapoints from calibration graphs"""

flowrate_c =  np.log10(1/ (df_pc_calibration_data['Flow rate [ml/min]'] * 1.667e-8))
potential_c = np.log10(df_pc_calibration_data['Potential'])
particlesize_c = df_pc_calibration_data['Particle Size [micron]']

# Create arbitrary mesh to interpolate on (note: not required to use the function, only to plot)
xx, yy = np.meshgrid(flowrate_c, potential_c)

# Create a function which returns the interpolated value for particle size at speicifc flowrates and potential
particlesize_fun = interpolate.interp2d(flowrate_c, 
                                          potential_c,
                                          particlesize_c, kind='linear')


#%% 
"""Create new matrix with interpolated grainsizes with the same dimesnsions as the output file
# Set potentials from datafile to list to be used for interpolation"""

final_potentials = df_pc_particles.columns.tolist()
final_potentials.remove('Datetime')
final_potentials.remove('Datetime_float')

final_flowrates = df_pc_final['Interpolated Flowrates [s/m3]'].tolist()

# Interpolate for each flowrate and potential to create mxn matrix containing the correct 
# Particle sizes for the corresponding time instances
final_particlesize = np.zeros((len(final_flowrates), len(final_potentials)))                                                               
for j in range (len(final_potentials)):
    for i in range (len(final_flowrates)):
        final_particlesize [i,j] = particlesize_fun(np.log10(final_flowrates[i]), np.log10(final_potentials[j]))

# Round up particle sizes to use more efficiently?
final_particlesize = pd.DataFrame(final_particlesize, columns = final_potentials)
# add np.log10

# Add to final dataset7
df_pc_final = pd.concat([df_pc_final, final_particlesize], axis=1)


#%% 
# """Create set of functions describing the number of particles and their corresponding
# size at a particular time instance (t1, t2 ....)"""

# Rename columns: 
df_pc_final.columns = ['Datetime', 
                    'ch1_n', 'ch2_n', 'ch3_n', 'ch4_n', 'ch5_n', 'ch6_n', 'ch7_n', 'ch8_n',
                    'Datetime_float', 'Interpolated_flowrates',
                    'ch1_p', 'ch2_p', 'ch3_p', 'ch4_p', 'ch5_p', 'ch6_p', 'ch7_p', 'ch8_p' ]


# Retrieve correct data from dataframe

no_particles = np.zeros((len(final_flowrates), len(final_potentials)))  
particle_size = np.zeros((len(final_flowrates), len(final_potentials)))  

for i in range(len(no_particles)):
        no_particles [i] = (df_pc_final.loc[[i],
               ['ch1_n', 'ch2_n', 'ch3_n', 'ch4_n', 'ch5_n', 'ch6_n', 'ch7_n', 'ch8_n']])

        particle_size [i] = (df_pc_final.loc[[i],
               ['ch1_p', 'ch2_p', 'ch3_p', 'ch4_p', 'ch5_p', 'ch6_p', 'ch7_p', 'ch8_p']])


# create partcile size column
df_pc_final ['Measured_particle_size'] = 5

# run the final interpolation
df_pc_final['No_Particles'] = df_pc_final.apply(lambda row:
                        interpolate.interp1d([row.ch1_p, row.ch2_p, row.ch3_p, row.ch4_p, row.ch5_p, row.ch6_p, row.ch7_p, row.ch8_p],
                                             [row.ch1_n, row.ch2_n, row.ch3_n, row.ch4_n, row.ch5_n, row.ch6_n, row.ch7_n, row.ch8_n],
                                  bounds_error=False
                        )(row.Measured_particle_size),
                        axis=1)
    
df_pc_final['Cummulative_Particles'] = df_pc_final['No_Particles'].cumsum()

df_pc_final['Interpolated Flowrates [ml/min]'] =  time_interpolation_fun(df_pc_final['Datetime_float'])

#%% Create new output file
#df_pc_final.to_excel(r'D:/Users/eimerno/Documents/Python Scripts/New_data1.xlsx')

#%% Plotting
plt.figure(1)    
plt.plot(df_pc_final['Datetime'], df_pc_final['No_Particles'], label='Number of Particles')
#plt.plot(df_pc_final['Datetime'], df_pc_final['Cummulative_Particles'], label = 'Cummulative Number of Particles')
plt.xlabel('Time')
plt.ylabel('Number of Particles')
plt.title('Particles Measured by Particle Counter')
plt.legend(loc = 'upper left')

plt.figure(2)
plt.plot(df_pc_final['Datetime'], df_pc_final['Cummulative_Particles'], label = 'Cummulative No. Particles')
plt.xlabel('Time')
plt.ylabel('Number of Particles')
plt.title('Cummulative Number of Particles Measured by Particle Counter')
plt.legend()

plt.figure(3)
plt.plot(df_pc_final['Datetime'], df_pc_final['Interpolated Flowrates [ml/min]'] ,label = 'Flowrates [ml/min]')
plt.xlabel('Time')
plt.ylabel('Flowrate [ml/min]')
plt.title('Flowrate Evolution over Time')