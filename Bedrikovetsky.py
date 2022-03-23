# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 12:23:32 2022

@author: eimerno
"""
# Bedrikovetsky, Electrostatic Potentials and Force

import numpy as np
import matplotlib.pyplot as plt


#%% Parameters
hi = 3.99e-9             # seperation distance (delta in paper)
A123 = 2e-20             # Hamaker constant
dp = 0.884e-6           # particle diameter
rs = 0.442e-6           # particle radius
psi01 = -30e-3          # surface potential 1 [V]
psi02 = -50e-3          # surface potential 2 [V]
eps0 = 8.854e-12        # free space permittivity [C-2 J-1 m-1]
De = 78.0               # dielectric constant 
Na =  6.022e23          # Avogadros numebr [mol -1] 
I = 100e-3              #ionic concentration   
#Sulfate and chloride were always present, often in large amounts, with concentrations ranging from 0.11 to 62.93 mg/L a
k = 1.3805e-23          # Boltzmann constant [J/K]
e = 1.6 * 10**-19       # electron charge [C]
chi = 89.5              # correction factor in lifting formula
omega = 60.0            # correction factor in equation for drag force
u = 0.01                # average flow velocity [m/s]
sigLJ = 0.5e-9          # Lennard-Jones potential

T = 293


#%%  Determine Functions V/derV

kappa = np.sqrt((2 * e**2 * Na * I)/(De * eps0 * k * T))
#kappa = 40000000
def FunV (h):
    Z = h/rs
    Vbr = (A123/7560) * (sigLJ/rs)**6 * (((8+Z)/(2+Z)**2) + ((6-Z)/Z**7))
    Vlva = -(A123/6) * ((2*(1+Z))/(Z*(2+Z)) + np.log(Z/(2+Z)))
    Vdlr1 = 2 * psi01 * psi02 * np.log((1+np.exp(-kappa*h))/(1-np.exp(-kappa*h)))
    Vdlr2 = (psi01**2 + psi02**2) * np.log(1-np.exp(-2*kappa*h)) 
    Vdlr = (eps0*De*rs/4) * (Vdlr1-Vdlr2)
    V = Vdlr + Vlva + Vbr
    return V
print(FunV(3.99e-4))

def DerFunV (h, dh):
    V1 = FunV(h + dh)
    V2 = FunV(h - dh)
    derV = (V1 - V2) / dh
    return derV  

#%%
h = np.linspace(1e-10, 5e-9)
V = np.zeros(len(h))
Fel = np.zeros(len(h))

for i in range (len(h)):
    V[i] = FunV(h[i])               # Electrostatic Potential 
    Fel[i] = - DerFunV(h[i], 1e-11)   # Electrostatic Force is dV/dh
    


plt.figure(1, figsize = (10,6) )
plt.plot(h,V, 'r', label = 'Electrostatic Potential')
plt.xlabel ('surface seperation h [m]')
plt.ylabel ('electrostatic potential V [J]')
#plt.title ('Electrostatic Potential as a function of the Seperation Distance' )
plt.grid()
plt.ylim (-8e-19, 3e-19)
plt.legend()

plt.figure(2, figsize = (10,6) )
plt.plot(h,Fel,'b', label='Electrostatic Force')
plt.xlabel ('surface seperation h [m]')
plt.ylabel ('electrostatic Force Fel [N]')
#plt.title ('Electrostatic Potential as a function of the Seperation Distance' )
plt.grid()
plt.ylim (-1.5e-9, 2.5e-9)
plt.legend()