# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:46:19 2022

@author: eimerno
"""
import numpy as np
#import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd 

# This script returns the value of the particle size passing through a particle counter
# To find the corresponding particle size, call the function:
#     particlesize_func('Transformed flow rate [m3/s]', 'Potential [Volt]')
# The transformed flow rate is defined as (1/Flow Rate [m3/s])

# For Kriging we need x,y,z coordinates: in this example
# x = flow rates
# y = potential
# z = particle size 


#%% Input parameters
actual_flowrate = 10 * 1.6667e-8  # [m3/s]
transformed_flowrate = 1 / actual_flowrate

#%% Input datasets
# Calibration data
calibration_data = pd.read_excel('Interpolation_Data.xlsx')

# Potentials used for experiment at set flow rate
experiment_data = pd.read_excel('Book2.xlsx', skiprows = 14)


#%% Create intetpolation function
# Define data from provided excel sheet containing all known datapoints from calibration graphs

flowrate_c = calibration_data.iloc[:,2]
potential_c = calibration_data.iloc[:,4]
particlesize_c = calibration_data.iloc[:,3]

# Create arbitrary mesh to interpolate on (note: not required to use the function, only to plot)
xx, yy = np.meshgrid(flowrate_c, potential_c)

# Create a function which returns the interpolated value for particle size at speicifc flowrates and potential
particlesize_fun = interpolate.interp2d(flowrate_c, 
                                          potential_c,
                                          particlesize_c, kind='linear')

#%% Use experiment data to find new particle size
potential_initial_ex = experiment_data.iloc [0,:]
potential_initial_ex.dropna(inplace = True)
potential_ex = potential_initial_ex.tolist()


final_particlesize = np.zeros(len(potential_ex))
for i in range(len(potential_ex)):
    final_particlesize [i] = particlesize_fun(transformed_flowrate, potential_ex[i])

print ('Final Particlesize is:', final_particlesize)

