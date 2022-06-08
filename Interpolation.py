# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:46:19 2022

@author: eimerno
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd 
import plotly.express as px
import datetime
from scipy.interpolate import fitpack
import os

# This script returns the value of the particle size passing through a particle counter
# To find the corresponding particle size, call the function:
#     particlesize_func('Transformed flow rate [m3/s]', 'Potential [Volt]')
# The transformed flow rate is defined as (1/Flow Rate [m3/s])

# For Kriging we need x,y,z coordinates: in this example
# x = flow rates
# y = potential
# z = particle size 


os.chdir('D:/Users/eimerno/Documents/Python Scripts/Interpolation')

#%% Input datasets
# Calibration data (Note: not necessary to change by user)
calibration_data = pd.read_excel('Interpolation_Data.xlsx')

# Potentials used for experiment at set flow rate 
experiment_data = pd.read_excel('Book3.xlsx', skiprows = 14, parse_dates=True)

# Flowrates
flowrates = pd.read_excel('Flowrates.xlsx')


#%% Create intetpolation function
# Define data from provided excel sheet containing all known datapoints from calibration graphs
flowrate_c = (1/(calibration_data['Flow rate [ml/min]'] * 1.667e-8)).tolist()
potential_c = calibration_data['Potential'].tolist()
particlesize_c = calibration_data['Particle Size [micron]'].tolist()


#Create a function which returns the interpolated value for particle size at speicifc flowrates and potential
# particlesize_fun = interpolate.interp2d(flowrate_c, 
#                                           potential_c,
#                                           particlesize_c, kind='linear')

#%%
def interpolation(flowrate, potential):
        # Need to write a norm function that calculates distance from a rib...

        # Segfaults... Problems with the way scipy is compiled?
       tck = interpolate.bisplrep(flowrate_c, potential_c, particlesize_c, task = 1, s=20)
       zi = interpolate.bisplev(flowrate, potential, tck, dx=0, dy=0)

       return zi
   
#%% plotting
# Create arbitrary mesh to interpolate on (note: not required to use the function, only to plot)

x_new = np.linspace (1e6, 3e7, 1000)
y_new = np.linspace (0.2, 8.7, 1000)

X, Y = np.meshgrid(x_new, y_new)

z_new = interpolation(x_new, y_new)
fig = plt.pcolormesh(X, Y, z_new)
plt.xlabel ('Flowrate [s/m3]')
plt.ylabel ('Potential') 
plt.colorbar()




# #%% Use experiment data to find new particle size
# # potential_initial_ex = experiment_data.iloc [0,:]
# # potential_initial_ex.dropna(inplace = True)
# # potential_ex = potential_initial_ex.tolist()

# potential_ex = experiment_data.columns.tolist()
# potential_ex.remove('Unnamed: 0')

# # Calculates the actual particle size for a given potential. 
# # These values are stores in the list: 'final_particlesize'
# final_particlesize = np.zeros(len(potential_ex))
# for i in range(len(potential_ex)):
#     final_particlesize [i] = interpolation(transformed_flowrate, potential_ex[i])

# # The correct particle sizes are translated back to the dataframe
# experiment_data.iloc[0,1] = final_particlesize[0]
# experiment_data.iloc[0,2] = final_particlesize[1]
# experiment_data.iloc[0,3] = final_particlesize[2]
# experiment_data.iloc[0,4] = final_particlesize[3]
# experiment_data.iloc[0,5] = final_particlesize[4]
# experiment_data.iloc[0,6] = final_particlesize[5]
# experiment_data.iloc[0,7] = final_particlesize[6]
# experiment_data.iloc[0,8] = final_particlesize[7]

# # Now the dataframa is equipt with the correct particle sizes and the number of particles at specific time isntances.
# #%% Create new output file
# experiment_data.to_excel(r'D:/Users/eimerno/Documents/Python Scripts/New_data1.xlsx')


# #%% Plotting the time series

# experiment_data.set_index('Unnamed: 0', inplace = True)
# fig = px.line(experiment_data, x=experiment_data.index)