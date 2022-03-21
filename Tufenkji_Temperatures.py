# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 14:48:50 2022

@author: eimerno
"""

# Tufenkji Single Contact Collector Efficiency Temperature Dependency
# Temperature varies, viscosity varies as a function of temperature 
import numpy as np
import matplotlib.pyplot as plt


#%% Parameters
d_p = 0.01 * 10**-6  # particle diameter [m]
d_c = 0.6 * 10**-3  # collector diameter [m] 
U = 9.0 * 10**-6       # fluid approach velocity [m/s]
A = 1 * 10**-20      # Hamaker constant [J]
rho_p = 1050.0        # particle density [kg/m3]
rho_f = 1000.0        # fluid density [kg/m3]
g = 9.81            # gravitational acceleration [m2/s]
#T = 353            # temperature [K]
n =0.39             # porosity
k = 1.3805 * 10**-23 # Boltzmann constant [J/K]
#mu = 0.0014176   # Dynamic viscosity water [Pa s]
a_p = 0.5* d_p      # Particle radius [m]
pi = np.pi



#%% Temperature Dependent Dimensionless parameters
T1 = np.linspace(273, 363)
Dinf = np.zeros_like(T1)
Nvdw = np.zeros_like(T1)
Ngr = np.zeros_like(T1)
mu = np.zeros_like(T1)
Na = np.zeros_like(T1)
N_g = np.zeros_like(T1)

for i in range (len(T1)):
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T1[i]-140))
    Na [i] = A / (12 * pi * mu[i] * a_p**2 * U) 
    Dinf[i] = (k*T1[i]) / (6*pi*mu[i]*a_p)
    Nvdw[i] = A / (k*T1[i])
    Ngr[i] = (4/3) * ((pi*(a_p**4)*(rho_p - rho_f)*g)/(k*T1[i]))
    N_g[i] = (2/9) * ((a_p**2*(rho_p - rho_f)*g)/(mu[i]*U)) 

#Dinf = (k*T) / (6*pi*mu*a_p)
Nr = d_p / d_c          # Aspect ratio
Npe = (U * d_c) /  Dinf
#Nvdw = A / (k*T)
#Ngr = (4/3) * ((pi*a_p**4*(rho_p - rho_f)*g)/(k*T))
Na = A / (12 * pi * mu * a_p**2 * U)
N_g = (2/9) * ((a_p**2*(rho_p - rho_f)*g)/(mu*U)) 

#%% Single grain collector efficiency
gamma = (1-n)**(1/3)
As = (2*(1-gamma**5))/(2-(3*gamma)+(3*gamma**5)-(2*gamma**6))

nD = 2.4 * As**(1/3) * Nr**-0.081 * Npe**-0.715 * Nvdw**0.052
nI = 0.55 * As * Nr**1.55 * Npe**-0.125 * Nvdw**0.125
nG = 0.22 * Nr**-0.24 * N_g**1.11 * Nvdw**0.053 

n0 = nD + nI + nG
print (n0)  

plt.plot(T1,n0)
plt.xlabel('Absolute temperature [K]')
plt.ylabel('n0')
plt.title('Temperature dependent single-collector contact efficiency (n0)');