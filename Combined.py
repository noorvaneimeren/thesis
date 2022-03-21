# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 09:46:53 2022

@author: eimerno
"""

# Combination of different model to formulate overview of forces acting on particle before detachment

import numpy as np
import matplotlib.pyplot as plt

#%% Parameters
h = 3.99e-9             # seperation distance (delta in paper)
H = 1.4e-20             # Hamaker constant
dp = 0.884e-6           # particle diameter
rs = 0.442e-6           # particle radius
psi01 = -30e-3          # surface potential 1 [V]
psi02 = -50e-3          # surface potential 2 [V]
eps0 = 8.854e-12        # free space permittivity [C-2 J-1 m-1]
De = 78.0               # dielectric constant 
Na =  6.022e23          # Avogadros numebr [mol -1] 
I = 0.3e-3              #ionic concentration   
#Sulfate and chloride were always present, often in large amounts, with concentrations ranging from 0.11 to 62.93 mg/L a
k = 1.3805 * 10**-23    # Boltzmann constant [J/K]
e = 1.6 * 10**-19       # electron charge [C]
chi = 89.5              # correction factor in lifting formula
omega = 60.0            # correction factor in equation for drag force
u = 0.01               # average flow velocity [m/s]
#T = 273
T = np.linspace(273,363)

#%%

#%% FlvdW = LvdW Force
FlvdW1 = (H * dp) / (12 * h**2)
FlvdW = np.ones(len(T))*FlvdW1


#%% Fedl = EDL force

kappa = np.zeros (len(T))
a = np.zeros(len(T))
b = np.zeros(len(T))
d = np.zeros(len(T))
for i in range(len(T)):
    kappa[i] = np.sqrt ((2 * e**2 * Na * I)/(De * eps0 * k * T[i]))
    a[i] = 2* np.pi * eps0 * De * (psi01**2+psi02**2) * kappa[i] * np.exp(-kappa[i] * h)
    b[i] = (1-np.exp(-2*kappa[i]*h))
    c = (2*psi01*psi02)/(psi01**2+psi02**2)
    d[i] = np.exp(-kappa[i] * h)

Fedl= (a/b) * (c-d)
#%% Net electrostatic adhesive force 

Fadh = FlvdW + Fedl

# Note: Adhesive force becomes more repellent with increasing temperatures

#plt.figure (1)
#plt.plot(T,Fadh,label='Net Adhesive Force (FvdW+FEDL)')
#plt.title ('Net Adhesive Force (FvdW+FEDL)')
#plt.plot(T,Fl, label = 'LvdW')
#plt.xlabel ('Temperature [K]')
#plt.ylabel ('F [N]')
#plt.legend()mu = np.zeros(len(T1))
#%% Fd Drag force 
mu = np.zeros(len(T))
for i in range (len(T)):
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    Fd = (omega*np.pi*mu*(rs)**2*u)/h 

#%% Fl = Lift Force
rho_f = np.zeros(len(T))
Fl = np.zeros(len(T))

for i in range (len(T)): 
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    rho_f [i] = 1e3 + 6.76e-2*(T[i]-273) - 8.99e-3*(T[i]-273)**2 + 9.14e-5*(T[i]-273)**3
    Fl[i]= chi * rs**3 * np.sqrt((mu[i]*rho_f[i]*u**3)/(h**3))
    
#%% Force Balance
Fmob = Fl + Fd
Fnet = Fadh+Fl+Fd
#%% Plotting
plt.figure (1)
plt.plot(T, Fnet, label = 'Net Force')
#plt.plot(T,Fadh,label='Net Adhesive Force (FvdW+FEDL)')
#plt.plot (T, Fnet, label = 'Net Detachment Force (Fd+Fl)')
plt.title ('Net Force Detachment')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.legend()


