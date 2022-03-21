# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 10:55:54 2022

@author: eimerno
"""

import numpy as np
#import sympy as sp
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
#u = 10             # average flow velocity [m/s]
u = np.linspace(0.01,10)
#T = 273
T = np.linspace(273,363)
#%% Flift
mu = np.zeros(len(T))
Fd1 = np.zeros(len(T))

u1 = 0.01
u2 = 0.5
u3 = 1.0
u4 = 1.75
u5 = 3
u6 = 4.5
u7 = 6
u8 = 8 

rho_f = np.zeros(len(T))


Fl1 = np.zeros(len(T))
for j in range (len(u)):
    plt.figure()
    for i in range (len(T)): 
        mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
        rho_f [i] = 1e3 + 6.76e-2*(T[i]-273) - 8.99e-3*(T[i]-273)**2 + 9.14e-5*(T[i]-273)**3
        Fl1[i]= chi * rs**3 * np.sqrt((mu[i]*rho_f[i]*u[j]**3)/(h**3))
    plt.plot(T,Fl1)
    plt.xlabel ('Temperature [K]')
    plt.ylabel ('Force[N]')











plt.figure(1)
plt.plot(T,Fl1, label = 'u = 0.01 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Lift versus Temperature') 

Fl2 = np.zeros(len(T))
for i in range (len(T)): 
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    rho_f [i] = 1e3 + 6.76e-2*(T[i]-273) - 8.99e-3*(T[i]-273)**2 + 9.14e-5*(T[i]-273)**3
    Fl2[i]= chi * rs**3 * np.sqrt((mu[i]*rho_f[i]*u2**3)/(h**3))
    
plt.plot(T,Fl2, label = 'u = 0.5 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Lift Force versus Temperature')
plt.legend()

Fl3 = np.zeros(len(T))
for i in range (len(T)): 
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    rho_f [i] = 1e3 + 6.76e-2*(T[i]-273) - 8.99e-3*(T[i]-273)**2 + 9.14e-5*(T[i]-273)**3
    Fl2[i]= chi * rs**3 * np.sqrt((mu[i]*rho_f[i]*u3**3)/(h**3))
    
plt.plot(T,Fl3, label = 'u = 1.0 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Lift Force versus Temperature')
plt.legend()

Fl4 = np.zeros(len(T))
for i in range (len(T)): 
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    rho_f [i] = 1e3 + 6.76e-2*(T[i]-273) - 8.99e-3*(T[i]-273)**2 + 9.14e-5*(T[i]-273)**3
    Fl2[i]= chi * rs**3 * np.sqrt((mu[i]*rho_f[i]*u4**3)/(h**3))
    
plt.plot(T,Fl4, label = 'u = 1.75 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Lift Force versus Temperature')
plt.legend()

Fl5 = np.zeros(len(T))
for i in range (len(T)): 
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    rho_f [i] = 1e3 + 6.76e-2*(T[i]-273) - 8.99e-3*(T[i]-273)**2 + 9.14e-5*(T[i]-273)**3
    Fl2[i]= chi * rs**3 * np.sqrt((mu[i]*rho_f[i]*u5**3)/(h**3))
    
plt.plot(T,Fl2, label = 'u = 3.0 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Lift Force versus Temperature')
plt.legend()

Fl6 = np.zeros(len(T))
for i in range (len(T)): 
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    rho_f [i] = 1e3 + 6.76e-2*(T[i]-273) - 8.99e-3*(T[i]-273)**2 + 9.14e-5*(T[i]-273)**3
    Fl6[i]= chi * rs**3 * np.sqrt((mu[i]*rho_f[i]*u6**3)/(h**3))
    
plt.plot(T,Fl6, label = 'u = 4.5 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Lift Force versus Temperature')
plt.legend()

Fl7 = np.zeros(len(T))
for i in range (len(T)): 
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    rho_f [i] = 1e3 + 6.76e-2*(T[i]-273) - 8.99e-3*(T[i]-273)**2 + 9.14e-5*(T[i]-273)**3
    Fl7[i]= chi * rs**3 * np.sqrt((mu[i]*rho_f[i]*u7**3)/(h**3))
    
plt.plot(T,Fl7, label = 'u = 6.0 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Lift Force versus Temperature')
plt.legend()

Fl8 = np.zeros(len(T))
for i in range (len(T)): 
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    rho_f [i] = 1e3 + 6.76e-2*(T[i]-273) - 8.99e-3*(T[i]-273)**2 + 9.14e-5*(T[i]-273)**3
    Fl8[i]= chi * rs**3 * np.sqrt((mu[i]*rho_f[i]*u8**3)/(h**3))
    
plt.plot(T,Fl8, label = 'u = 8.0 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Lift Force versus Temperature')
plt.legend()


#%% Fd

mu = np.zeros(len(T))


u1 = 0.01
u2 = 0.5
u3 = 1.0
u4 = 1.75
u5 = 3
u6 = 4.5
u7 = 6
u8 = 8 

rho_f = np.zeros(len(T))

mu=np.zeros(len(T))


Fd1 = np.zeros (len(T))
for i in range (len(T)):
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    Fd1 [i]= (omega*np.pi*mu[i]*(rs)**2*u1)/h
    
plt.figure(2)
plt.plot(T,Fd1, label = 'u = 0.01 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Drag Force versus Temperature') 

Fd2 = np.zeros (len(T))
for i in range (len(T)):
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    Fd2 [i]= (omega*np.pi*mu[i]*(rs)**2*u2)/h
    

plt.plot(T,Fd2, label = 'u = 0.5 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Drag Force versus Temperature') 

Fd3 = np.zeros (len(T))
for i in range (len(T)):
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    Fd2 [i]= (omega*np.pi*mu[i]*(rs)**2*u3)/h
    

plt.plot(T,Fd3, label = 'u = 1.0 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Drag Force versus Temperature') 

Fd4 = np.zeros (len(T))
for i in range (len(T)):
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    Fd4 [i]= (omega*np.pi*mu[i]*(rs)**2*u4)/h
    

plt.plot(T,Fd4, label = 'u = 1.75 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Drag Force versus Temperature') 

Fd5 = np.zeros (len(T))
for i in range (len(T)):
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    Fd5 [i]= (omega*np.pi*mu[i]*(rs)**2*u5)/h
    

plt.plot(T,Fd5, label = 'u = 3 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Drag Force versus Temperature') 

Fd6 = np.zeros (len(T))
for i in range (len(T)):
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    Fd6 [i]= (omega*np.pi*mu[i]*(rs)**2*u6)/h
    

plt.plot(T,Fd6, label = 'u = 4.5 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Drag Force versus Temperature') 


Fd7 = np.zeros (len(T))
for i in range (len(T)):
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    Fd7 [i]= (omega*np.pi*mu[i]*(rs)**2*u7)/h
    

plt.plot(T,Fd7, label = 'u =6.0 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Drag Force versus Temperature') 

Fd8 = np.zeros (len(T))
for i in range (len(T)):
    mu[i] = 2.4318*10**-5 * 10**(247.8/(T[i]-140))
    Fd8 [i]= (omega*np.pi*mu[i]*(rs)**2*u8)/h
    

plt.plot(T,Fd8, label = 'u =8.0 m/s')
plt.xlabel ('Temperature [K]')
plt.ylabel ('Force [N]')
plt.title ('Drag Force versus Temperature') 
plt.legend()