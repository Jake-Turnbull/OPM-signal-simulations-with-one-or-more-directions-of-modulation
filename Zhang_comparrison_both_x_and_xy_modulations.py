# -*- coding: utf-8 -*-
"""
Created on Mon May 20 22:23:17 2024

@author: jaket
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.special
from scipy.integrate import odeint
from matplotlib.pyplot import Slider

"""Initialise the values needed for the simulation"""
gamma = 2 * 7 * np.pi
nsteps = 170
w1 = 10 * gamma
w = 2 * np.pi * 1000
ncycles = 100# Define the larmour frequency gamma
P0 = 1 #starting polarisation
R = 10 * gamma #Relax[0]ation
Bmax= 25
q= 1
dt = 0.000001
t = (np.arange(0, ncycles*nsteps, 1)+0.7 )*dt
ttest=(np.linspace(0,0.017,100))*dt
size = nsteps*ncycles

val_step = Bmax/50
B = np.arange(-Bmax, Bmax, val_step)
By=0
Bz=0
Wz = Bz * gamma
Wy = By * gamma
"""Define the inner part of our bessel functions, since it only requires the magnitude of B_mod
and not the time component we can express it simply as 10nT"""
def Beta(gamma,Bm,q,omega):
    return (gamma*Bm)/(q*omega)

Bxm = w1
"""store our bessel functions as distinct values to be qouted later"""
Jv0 =scipy.special.jv(0,Beta(gamma,Bxm,q,w))
Jv1 =scipy.special.jv(1,Beta(gamma,Bxm,q,w))

"""create two new variables, response and response_null. Response will have its value changed by
the sliders in the GUI made later and Response_null is our reference to calculate amplitude difference"""

Response = 0.82*Jv0*Jv1*((R*gamma*B+Jv0**2*By*Bz))/(R**2+(gamma*B)**2+Jv0**2*gamma*By**2+Jv0**2*gamma**2*Bz**2)
Response_null = 0.82*Jv0*Jv1*((R*gamma*B+Jv0**2*0*0))/(R**2+(gamma*B)**2+Jv0**2*gamma*0**2+Jv0**2*gamma**2*0**2)
ss = np.zeros(nsteps * ncycles)
cc = np.zeros(nsteps * ncycles)

"""Initialise empty arrays"""

Pxval = np.zeros(100)
Pyval = np.zeros(100)
Pzval = np.zeros(100)
valzs = np.zeros(100)
valzc = np.zeros(100)

Pxval2 = np.zeros(100)
Pyval2 = np.zeros(100)
Pzval2 = np.zeros(100)
valzs2 = np.zeros(100)
valzc2 = np.zeros(100)

"""Simulate our reference signal"""
for j in range(np.size(t)):
        ss[j] = np.sin(w*t[j])
        cc[j] = np.cos(w*t[j])
"""Define figure size and slider axes"""
fig1,ax = plt.subplots(1,2,figsize =(16,16))
By_slider = fig1.add_axes([0.93, 0.2, 0.01, 0.6])
Bz_slider = fig1.add_axes([0.95, 0.2, 0.01, 0.6])


"""Define our ODE solver"""
def ODE(vec,t,R,gamma,Wx,Wy,Wz):
    gamma = 2 * 7 * np.pi # Define the larmour frequency gamma
    R = 10 * gamma #Relax[0]ation
    #Defined constants^
    #assign each ODE to vector element
    Px = vec[0]
    Py = vec[1]
    Pz = vec[2]
    #Initialise the equations:      
    dPxdt = ((0*np.sin(w*t)+Wy) * Pz - Wz * Py - R * Px)
    dPydt = (Wz * Px - (w1*np.cos(w*t)+Wx) * Pz - R * Py)
    dPzdt = ((w1*np.cos(w*t)+Wx) * Py - (0*np.sin(w*t)+Wy) * Px - R * (Pz - 1))
    
    Out = np.array([dPxdt, dPydt, dPzdt])
    
    return Out 

def ODE2(vec,t,R,gamma,Wx,Wy,Wz):
    gamma = 2 * 7 * np.pi # Define the larmour frequency gamma
    R = 10 * gamma #Relax[0]ation
    #Defined constants^
    #assign each ODE to vector element
    Px = vec[0]
    Py = vec[1]
    Pz = vec[2]
    #Initialise the equations:      
    dPxdt = ((w1*np.sin(w*t)+Wy) * Pz - Wz * Py - R * Px)
    dPydt = (Wz * Px - (w1*np.cos(w*t)+Wx) * Pz - R * Py)
    dPzdt = ((w1*np.cos(w*t)+Wx) * Py - (w1*np.sin(w*t)+Wy) * Px - R * (Pz - 1))
    
    Out = np.array([dPxdt, dPydt, dPzdt])
    
    return Out 

"""Set up loop to calculate values of ODE"""
for x_val in range(100):
    print(x_val)   
    Wx = B[x_val] * gamma
    denom = Wx**2 + Wy**2 + Wz**2 + R**2
    
    Pz_initial = (Wz**2 + R**2) / denom
    Py_initial = (R * Wx + Wz * Wy) / denom
    Px_initial = (Wz * Wx - Wy * R) / denom
    
    vec0 = np.array([Px_initial,Py_initial,Pz_initial])
    #test=ODE(vec0,t,R,gamma)
    val = odeint(ODE, vec0, t, args=(R,gamma,Wx,Wy,Wz)) 
    
    Pxval[x_val] = np.mean(val[(size-5000):size,0])
    Pyval[x_val] = np.mean(val[(size-5000):size,1])
    Pzval[x_val] = np.mean(val[(size-5000):size,2])

    z_modulated = val[:,2]
    valzs[x_val] = 2*np.mean(ss[(size-5000):size] * z_modulated[(size-5000):size])
    valzc[x_val] = 2*np.mean(cc[(size-5000):size] * z_modulated[(size-5000):size])
    
    #test=ODE(vec0,t,R,gamma)
    val2 = odeint(ODE2, vec0, t, args=(R,gamma,Wx,Wy,Wz)) 
    
    Pxval2[x_val] = np.mean(val[(size-5000):size,0])
    Pyval2[x_val] = np.mean(val[(size-5000):size,1])
    Pzval2[x_val] = np.mean(val[(size-5000):size,2])
    
    z_modulated2 = val2[:,2]
    
    valzs2[x_val] = 2*np.mean(ss[(size-5000):size] * z_modulated2[(size-5000):size])
    valzc2[x_val] = 2*np.mean(cc[(size-5000):size] * z_modulated2[(size-5000):size])
    
ax[0].plot(B,Response,'.',label = 'Predicted form')
ax[0].plot(B,valzs,label='Simulated form')
ax[0].legend(loc='upper right')
ax[0].set_ylabel('Response',fontsize=12)
ax[0].set_xlabel('Magnetic field Bx (nT)',fontsize=12)

ax[1].set_ylabel('Response difference between null case',fontsize=12)
ax[1].plot(B,((valzs-Response_null)),label=f'Bz={Bz}, By={By}')
ax[1].plot(B,((valzs2-Response_null)),'--',label=f'Bz={round(Bz,2)}, By={round(By,2)}')
ax[1].legend(loc='upper left')
ax[1].set_xlabel('Magnetic field Bx (nT)',fontsize=12)

#%%

"""Define two sliders to respectively change By and Bz and plot the results and their
amlitude differences"""

def sliderCallbacky(i):
    global  Wz, Bz, By
    ax[0].clear()
    By = B[i]
    Wy = By * gamma
    
    for x_val in range(100):
        print(x_val)
        Wx = B[x_val] * gamma
        
        denom = Wx**2 + Wy**2 + Wz**2 + R**2
        Pz_initial = (Wz**2 + R**2) / denom
        Py_initial = (R * Wx + Wz * Wy) / denom
        Px_initial = (Wz * Wx - Wy * R) / denom
            
        vec0 = np.array([Px_initial,Py_initial,Pz_initial])
        #test=ODE(vec0,t,R,gamma)
        val = odeint(ODE, vec0, t, args=(R,gamma,Wx,Wy,Wz)) 
        
        Pxval[x_val] = np.mean(val[(size-5000):size,0])
        Pyval[x_val] = np.mean(val[(size-5000):size,1])
        Pzval[x_val] = np.mean(val[(size-5000):size,2])
        z_modulated = val[:,2]
        
        valzs[x_val] = 2*np.mean(ss[(size-5000):size] * z_modulated[(size-5000):size])
        valzc[x_val] = 2*np.mean(cc[(size-5000):size] * z_modulated[(size-5000):size])
        val2 = odeint(ODE2, vec0, t, args=(R,gamma,Wx,Wy,Wz)) 
        
        Pxval2[x_val] = np.mean(val[(size-5000):size,0])
        Pyval2[x_val] = np.mean(val[(size-5000):size,1])
        Pzval2[x_val] = np.mean(val[(size-5000):size,2])
        
        z_modulated2 = val2[:,2]
        
        valzs2[x_val] = 2*np.mean(ss[(size-5000):size] * z_modulated2[(size-5000):size])
        valzc2[x_val] = 2*np.mean(cc[(size-5000):size] * z_modulated2[(size-5000):size])
        
    Response = 0.82*Jv0*Jv1*((R*gamma*B+Jv0**2*By*Bz))/(R**2+(gamma*B)**2+Jv0**2*gamma*By**2+Jv0**2*gamma**2*Bz**2)
    ax[0].plot(B,Response,'.',label = 'Predicted form')
    ax[0].plot(B,valzs,label='Simulated form')
    ax[0].legend(loc='upper right')
    ax[0].set_ylabel('Response',fontsize=12)
    ax[0].set_xlabel('Magnetic field Bx (nT)',fontsize=12)
    
    ax[1].set_ylabel('Response difference between null case',fontsize=12)
    ax[1].plot(B,((valzs-Response_null)),label=f'Bz={round(Bz,2)}, By={round(By,2)}')
    ax[1].plot(B,((valzs2-Response_null)),'--',label=f'Bz={round(Bz,2)}, By={round(By,2)}')
    ax[1].legend(loc='upper left')
    ax[1].set_xlabel('Magnetic field Bx (nT)',fontsize=12)

    
def sliderCallbackz(i):
    global Wy, By, Bz
    ax[0].clear()
    Bz = B[i]
    Wz = Bz * gamma
    # Pxval = np.zeros(100)
    # Pyval = np.zeros(100)
    # Pzval = np.zeros(100)
    # valzs = np.zeros(100)
    # valzc = np.zeros(100)
    for x_val in range(100):
        print(x_val)
        Wx = B[x_val] * gamma
        
        denom = Wx**2 + Wy**2 + Wz**2 + R**2
        Pz_initial = (Wz**2 + R**2) / denom
        Py_initial = (R * Wx + Wz * Wy) / denom
        Px_initial = (Wz * Wx - Wy * R) / denom
            
        vec0 = np.array([Px_initial,Py_initial,Pz_initial])
        #test=ODE(vec0,t,R,gamma)
        val = odeint(ODE, vec0, t, args=(R,gamma,Wx,Wy,Wz)) 
        
        Pxval[x_val] = np.mean(val[(size-5000):size,0])
        Pyval[x_val] = np.mean(val[(size-5000):size,1])
        Pzval[x_val] = np.mean(val[(size-5000):size,2])
        z_modulated = val[:,2]
        
        valzs[x_val] = 2*np.mean(ss[(size-5000):size] * z_modulated[(size-5000):size])
        valzc[x_val] = 2*np.mean(cc[(size-5000):size] * z_modulated[(size-5000):size])
        
        val2 = odeint(ODE2, vec0, t, args=(R,gamma,Wx,Wy,Wz)) 
        
        Pxval2[x_val] = np.mean(val[(size-5000):size,0])
        Pyval2[x_val] = np.mean(val[(size-5000):size,1])
        Pzval2[x_val] = np.mean(val[(size-5000):size,2])
        
        z_modulated2 = val2[:,2]
        
        valzs2[x_val] = 2*np.mean(ss[(size-5000):size] * z_modulated2[(size-5000):size])
        valzc2[x_val] = 2*np.mean(cc[(size-5000):size] * z_modulated2[(size-5000):size])
        
    Response = 0.82*Jv0*Jv1*((R*gamma*B+Jv0**2*By*Bz))/(R**2+(gamma*B)**2+(Jv0**2*gamma*By**2)+(Jv0**2*gamma**2*Bz**2))
    ax[0].plot(B,Response,'.',label = 'Predicted form')
    ax[0].plot(B,valzs,label='Simulated form')
    ax[0].legend(loc='upper right')
    ax[0].set_ylabel('Response',fontsize=12)
    ax[0].set_xlabel('Magnetic field Bx (nT)',fontsize=12)
    
    ax[1].set_ylabel('Response difference between null case',fontsize=12)
    ax[1].plot(B,((valzs-Response_null)),label=f'Bz={round(Bz,2)}, By={round(By,2)}')
    ax[1].plot(B,((valzs2-Response_null)),'--',label=f'Bz={round(Bz,2)}, By={round(By,2)}')
    ax[1].legend(loc='upper left')
    ax[1].set_xlabel('Magnetic field Bx (nT)',fontsize=12)

By_slider_touch = Slider(By_slider, 
                 label='Y value',
                 valmin=0, 
                 valmax=100,
                 valstep=1,
                 valinit=50,
                 orientation="vertical")
By_slider_touch .on_changed(sliderCallbacky)

Bz_slider_touch  = Slider(Bz_slider, 
                 label='Z_value',
                 valmin=0, 
                 valmax=100,
                 valstep=1,
                 valinit=50,
                 orientation="vertical")
Bz_slider_touch .on_changed(sliderCallbackz)
