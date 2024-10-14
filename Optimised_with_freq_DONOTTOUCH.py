# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 11:41:36 2024

@author: jaket
"""

""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint

#%%
"""Set up constants and starting values as below"""
plt.close('all')
gamma = 2 * 7 * np.pi
nsteps = 1000
w1 = 10 * gamma
w = 2 * np.pi * 1000
ncycles = 100# Define the larmour frequency gamma
P0 = 1 #starting polarisation
R = 10 * gamma #Relaxation

B = np.arange(-100, 100, 2)
B2 = np.arange(-100,100, 2)

dt = 0.000001
By=50   
t = (np.arange(0, 100000, 1)+0.5 )* dt
#%%

fig1 = plt.figure(figsize =(16,12))
ax1 = fig1.add_axes([0.1,0.15,0.8,0.8])
ax1.set_ylabel('Px')
ax1.set_xlabel('Bx(nT)')
ax1.set_title('Polarisation in x against Bx')

fig2 = plt.figure(figsize =(16,12))
ax2 = fig2.add_axes([0.1,0.15,0.8,0.8])
ax2.set_ylabel('Py')
ax2.set_xlabel('Bx(nT)')
ax2.set_title('Polarisation in y against Bx')

fig3 = plt.figure(figsize =(16,12))
ax3 = fig3.add_axes([0.1,0.15,0.8,0.8])
ax3.set_ylabel('Pz')
ax3.set_xlabel('Bx(nT)')
ax3.set_title('Polarisation in z against Bx')

fig4 = plt.figure(figsize =(16,12))
ax4 = fig4.add_axes([0.1,0.15,0.8,0.8])
ax4.set_ylabel('Pz')
ax4.set_xlabel('Bx(nT)')
ax4.set_title('Sine term = diff B_z value')

fig5 = plt.figure(figsize =(16,12))
ax5 = fig5.add_axes([0.1,0.15,0.8,0.8])
ax5.set_ylabel('Pz')
ax5.set_xlabel('Bx(nT)')
ax5.set_title('Cosine term = diff B_z value')
#%%
Pxval= np.zeros(100)
Pyval= np.zeros(100)
Pzval= np.zeros(100)

ss = np.zeros(nsteps * ncycles)
cc = np.zeros(nsteps * ncycles)

Pxxval= np.zeros((100,100))
Pyyval= np.zeros((100,100))
Pzzval= np.zeros((100,100))

rows, cols = (100, 100)
valzs = np.zeros((rows, cols))
valzc = np.zeros((rows, cols))

for j in range(np.size(t)):
        #ic = i + nc
        #test.append(ic)
        # if ic==999:
        #     print('what')
        ss[j] = np.sin(w*t[j])
        cc[j] = np.cos(w*t[j])
#%%    
t = (np.arange(0, 100000, 1)+0.5) * dt
Wy=0

for z_ in range(10):
    z_val = z_*10
    
    print(z_val)
    Bz = B[z_val]
    Wz = gamma * Bz
    
    for x_val in range(100):
        
        Bx = B2[x_val]
        Wx = gamma * Bx
        
        denom = Wx**2 + Wy**2 + Wz**2 + R**2
        
        Pz_initial = (Wz**2 + R**2) / denom
        Py_initial = (R * Wx + Wz * Wy) / denom
        Px_initial = (Wz * Wx - Wy * R) / denom
        
        vec0 = np.array([Px_initial,Py_initial,Pz_initial])
          
        Wpx=w1*np.cos(w*t)+Wx

        Wpy=w1*np.sin(w*t)+Wy
        
        def ODE(vec,t,R,gamma):
        
            gamma = 2 * 7 * np.pi # Define the larmour frequency gamma
            R = 10 * gamma #Relaxation
            #Defined constants^
            #assign each ODE to vector element
            Px = vec[0]
            Py = vec[1]
            Pz = vec[2]
            #Initialise the equations:    
            
            dPxdt = ((w1*np.sin(w*t)+Wy) * Pz - gamma * Bz * Py - R * Px)
            dPydt = (gamma * Bz * Px - (w1*np.cos(w*t)+Wx) * Pz - R * Py)
            dPzdt = ((w1*np.cos(w*t)+Wx) * Py - (w1*np.sin(w*t)+Wy) * Px - R * (Pz - 1))
            
            #Define our array
            
            Out = np.array([dPxdt, dPydt, dPzdt])
        
            return Out          
        #test=ODE(vec0,t,R,gamma)
        
        val = odeint(ODE, vec0, t, args=(R,gamma))
            
        Pxval[x_val] = np.mean(val[:,0])
        Pyval[x_val] = np.mean(val[:,1])
        Pzval[x_val] = np.mean(val[:,2])
        
        z_modulated = val[:,2]
        
        valzs[z_val,x_val] = 2*np.mean(ss[90000:100000] * z_modulated[90000:100000])
        valzc[z_val,x_val] = 2*np.mean(cc[90000:100000] * z_modulated[90000:100000])
        
    Pxxval[z_val,:] = Pxval
    Pyyval[z_val,:] = Pyval
    Pzzval[z_val,:] = Pzval
    
    ax1.plot(B2,Pxval, label= f'Bz={B[z_*10]}nT')
    ax1.legend(bbox_to_anchor=(0.98, 1), loc='upper left', borderaxespad=0)
    
    ax2.plot(B2,Pyval, label= f'Bz={B[z_*10]}nT')
    ax2.legend(bbox_to_anchor=(0.98, 1), loc='upper left', borderaxespad=0)
    
    ax3.plot(B2,Pzval, label= f'Bz={B[z_*10]}nT')
    ax3.legend(bbox_to_anchor=(0.98, 1), loc='upper left', borderaxespad=0)
    
    ax4.plot(B2,valzs[z_val,:], label= f'Bz={B[z_*10]}nT')
    ax4.legend(bbox_to_anchor=(0.98, 1), loc='upper left', borderaxespad=0)
    
    ax5.plot(B2,valzc[z_val,:], label= f'Bz={B[z_*10]}nT')
    ax5.legend(bbox_to_anchor=(0.98, 1), loc='upper left', borderaxespad=0)
    
    plt.pause(0.5)