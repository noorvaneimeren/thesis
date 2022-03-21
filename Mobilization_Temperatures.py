# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 09:46:53 2022

@author: eimerno
"""

# Combination of different model to formulate overview of forces acting on particle before detachment
# 1. Bedrikovetsky (Modified Particle Detachment Model for Colloidal Transport in Porous Media): Lift force, Drag force
# 2. Ren Bai, Chi Ten (Particle Detachment in Deep bed Filtration): EDL/vdW forces
# 3. Van Esch: Temperature dependency density and viscosity 
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
#u = 10mday           # average flow velocity [m/s]
#u = np.linspace(0.01,10)
#T = 273
T = np.linspace(273,363)
u = [0.001,0.01]
#%% FvdW = VdW Force
FvdW1 =  H*dp/(12*h**2)
FvdW = np.ones(len(T)) * FvdW1 

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


#%% Fd Drag force 

mu=np.zeros(len(T))
Fd = np.zeros (len(T))

for jVel in range (len(u)):
    plt.figure(1)    
    for iTemp in range (len(T)):
        mu[iTemp] = 2.4318*10**-5 * 10**(247.8/(T[iTemp]-140))
        Fd [iTemp]= (omega*np.pi*mu[iTemp]*(rs)**2*u[jVel])/h
    plt.plot(T,Fd, label = u[jVel])
    plt.ylabel ('Force[N]')
    plt.legend()
    plt.title ('Drag Force versus Temperature for different intrinsic flow Velocities')
    
    
  
#%% Fl = Lift Force
rho_f = np.zeros(len(T))
Fl = np.zeros(len(T))

for jVel in range (len(u)):
    plt.figure(2)
    for iTemp in range (len(T)): 
        mu[iTemp] = 2.4318*10**-5 * 10**(247.8/(T[iTemp]-140))
        rho_f [iTemp] = 1e3 + 6.76e-2*(T[iTemp]-273) - 8.99e-3*(T[iTemp]-273)**2 + 9.14e-5*(T[iTemp]-273)**3
        Fl[iTemp]= chi * rs**3 * np.sqrt((mu[iTemp]*rho_f[iTemp]*u[jVel]**3)/(h**3))
    plt.plot(T,Fl, label= u[jVel])
    plt.xlabel ('Temperature [K]')
    plt.ylabel ('Force[N]')
    plt.legend()
    plt.title ('Lift Force versus Temperature for different intrinsic flow Velocities')
    
    
    
#%% Force Balance
Fmob= np.zeros (len(T))


for jVel in range (len(u)):
    plt.figure(3)
    for iTemp in range (len(T)):
        mu[iTemp] = 2.4318*10**-5 * 10**(247.8/(T[iTemp]-140))
        rho_f [iTemp] = 1e3 + 6.76e-2*(T[iTemp]-273) - 8.99e-3*(T[iTemp]-273)**2 + 9.14e-5*(T[iTemp]-273)**3
        Fd [iTemp]= (omega*np.pi*mu[iTemp]*(rs)**2*u[jVel])/h
        Fl[iTemp]= chi * rs**3 * np.sqrt((mu[iTemp]*rho_f[iTemp]*u[jVel]**3)/(h**3))
        Fmob[iTemp] = Fd[iTemp]+Fl[iTemp] + FvdW[iTemp] + Fedl[iTemp] 
    plt.plot(T,Fmob, label= u[jVel])
    plt.xlabel ('Temperature [K]')
    plt.ylabel ('Force[N]')
    plt.legend()
    plt.title ('Net Mob Force versus Temperature for different intrinsic flow Velocities')




#%% Plotting
plt.figure (4)
plt.plot(T, FvdW , label = 'Fedl')
#plt.plot(T, Fl, label = 'Fl')
#plt.plot(T,Fd,label='Fd')
#plt.plot (T, FvdW, label = 'FvdW')
#plt.title ('Forces')
#plt.xlabel ('Temperature [K]')
#plt.ylabel ('Force [N]')
#plt.legend()


#plt.figure(4)
#plt.plot (T, Fmob, label='Fnet')
#plt.plot (u, Fl, label='Fl')
#plt.title ('Net Force on Particle')
#plt.xlabel ('Temperature [K]')
#plt.ylabel ('Force [N]')
#plt.legend()

plt.figure (5)
plt.plot(T, Fedl, label = 'Fedl')
plt.plot(T,FvdW, label='FvdW')
plt.plot (T, Fd, label='Fd')
plt.plot (T, Fl, label='Fl')
plt.title ('Particle Forces')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.legend()