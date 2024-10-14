# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 18:48:55 2024

@author: jaket
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint

#%%
"""Set up constants and starting values as below"""
plt.close('all')
gamma = 2 * 7 * np.pi #Define gamma (larmour frequency)
nsteps = 200 #Define number of steps for sine and cosine calculations
w1 = 10 * gamma
w = 2 * np.pi * 1000
ncycles = 100

R = 10 * gamma #Relaxation

Range = 2
Range2 = 5
B_step = Range/50
B_step2 = Range2/50


B = np.arange(-Range, Range, B_step)
B2 = np.arange(-Range2,Range2, B_step2)

dt = 0.000001 
t = (np.arange(0, ncycles*nsteps, 1)) * dt

#%%
#This section just creates the figures and labels them

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

#Set up our arrays, Pxval and its respectives are one dimensional arrays that are updated for each value of Bz then over
#written for each new run with a new Bz value

Pxval= np.zeros(100)
Pyval= np.zeros(100)
Pzval= np.zeros(100)

#set up our arrays to define sine and cosine
ss = np.zeros(nsteps * ncycles)
cc = np.zeros(nsteps * ncycles)
ns = (ncycles - 10) * nsteps + 1
nf = ncycles * nsteps

#These arrays are contain the complete graphs for all values of Bx and Bz
Pxxval= np.zeros((100,100))
Pyyval= np.zeros((100,100))
Pzzval= np.zeros((100,100))

rows, cols = (100, 100)
valzs = np.zeros((rows, cols))
valzc = np.zeros((rows, cols))

#Removing the sine and cosine components so that the modulation can be extracted from Pz allows the code to run faster
#as values don't change each time it loops
for j in range(np.size(t)):
        #ic = i + nc
        #test.append(ic)
        # if ic==999:
        #     print('what')
        ss[j] = np.sin(w*t[j])
        cc[j] = np.cos(w*t[j])
#%%    

def ODE(vec,t,R,gamma):

    gamma = 2 * 7 * np.pi # Define the larmour frequency gamma
    R = 10 * gamma #Relaxation
    #Defined constants^
    
    #assign each ODE to vector element
    Px = vec[0]
    Py = vec[1]
    Pz = vec[2]
    
    #Initialise the equationf-5000:    
    dPxdt = ((0*np.sin(w*t)+Wy) * Pz - Wz * Py - R * Px)
    dPydt = (Wz * Px - (w1*np.cos(w*t)+Wx) * Pz - R * Py)
    dPzdt = ((w1*np.cos(w*t)+Wx) * Py - (0*np.sin(w*t)+Wy) * Px - R * (Pz - 1))
    
    #Define our output array in a managable form
    Out = np.array([dPxdt, dPydt, dPzdt])

    return Out
#Set the static field in y direction to zero
Wy = 0
#Create a loop that steps through values of Bz 
for z_ in range(10):
    z_val = z_*10
    
    print(z_val)
    Bz = B2[z_val]
    Wz = gamma * Bz
    #Create a second nested loop to calculate over 100 values of Bx for each value of By
    
    for x_val in range(100):
        
        Bx = B[x_val]
        Wx = gamma * Bx
        
        denom = Wx**2 + Wy**2 + Wz**2 + R**2
        Pz_initial = (Wz**2 + R**2) / denom
        Py_initial = (R * Wx + Wz * Wy) / denom
        Px_initial = (Wz * Wx - Wy * R) / denom
        
        vec0 = np.array([Px_initial,Py_initial,Pz_initial])
        
        #test=ODE(vec0,t,R,gamma)
       
        #Plug into the ODE solver 
        val = odeint(ODE, vec0, t, args=(R,gamma))
        
        #Extract each row of the array and take the mean of the
        #last 10,000 points to find the exact value of each polarisation
        
        Pxxval[z_val,x_val] = np.mean(val[:,0][nf-5000:nf])
        Pyyval[z_val,x_val] = np.mean(val[:,1][nf-5000:nf])
        Pzzval[z_val,x_val] = np.mean(val[:,2][nf-5000:nf])
        
        z_modulated = val[:,2]
        
        #Take the last 10,000 values of the sine and cosine forms and multiply with the polarization in z

        valzs[z_val,x_val] = 2*np.mean(ss[nf-5000:nf] * z_modulated[nf-5000:nf])
        valzc[z_val,x_val] = 2*np.mean(cc[nf-5000:nf] * z_modulated[nf-5000:nf])
    
    
    #Assign each of the arrays to a values of Pz
    
    ax4.plot(B2,valzs[z_val,:], label= f'Bz={round(B[z_*10])}nT')
    ax4.legend(bbox_to_anchor=(0.98, 1), loc='upper left', borderaxespad=0)
    
    ax5.plot(B2,valzc[z_val,:], label= f'Bz={round(B[z_*10])}nT')
    ax5.legend(bbox_to_anchor=(0.98, 1), loc='upper left', borderaxespad=0)
    
    #Plot the respective sine and cosine forms of Pz for the current value of Bz

#%%      
for num in range(10):
    #Plot the values of polarization in x, y and z directions, removed from the loop to increase efficiency
    
    ax1.plot(B2,Pxxval[num*10,:], label= f'Bz={B[num*10]}nT')
    ax1.legend(bbox_to_anchor=(0.98, 1), loc='upper left', borderaxespad=0)
    
    ax2.plot(B2,Pyyval[num*10,:], label= f'Bz={B[num*10]}nT')
    ax2.legend(bbox_to_anchor=(0.98, 1), loc='upper left', borderaxespad=0)
    
    ax3.plot(B2,Pzzval[num*10,:], label= f'Bz={B[num*10]}nT')
    ax3.legend(bbox_to_anchor=(0.98, 1), loc='upper left', borderaxespad=0)

#%%
#The below code plots the ideal solution for the same range in Bx and Bz as a way of checking the form of Pz

B_x = np.arange(-Range, Range, B_step)
B_y = 0#np.arange(-Range,Range,2)
B_z = np.arange(-Range, Range, B_step)
     
for i in range(10):   
    gamma = 2*7*np.pi;
    freq = 1000*2*7*gamma;
    
    R_op = 1
    R = gamma*10
    P_0 = 1
    relax = (R/gamma)**2;
    
    P_z = (relax + B_z[i*10]**2)/(relax + B_x**2 + B_y**2 +B_z[i*10]**2)*P_0
    ax3.plot(B_x,P_z, 'xb', markersize=1)

#%%
#The following section calculates the gradients of each Bz value used around low values of Bx and plots an estimate
#for the shape of the line at these values as well as the actual values for comparisson, this is done to allow for manual adjustment of the range over whicg the linear equation should be caclulated
plt.pause(1.0)

fig6 = plt.figure(figsize =(16,16))
ax6 = fig6.add_axes([0.1,0.15,0.8,0.8])
ax6.set_ylabel('Px')
ax6.set_xlabel('Bx(nT)')
ax6.set_title('Linear range of Sine components of Pz')

fig7 = plt.figure(figsize =(16,16))
ax7 = fig7.add_axes([0.1,0.15,0.8,0.8])
ax7.set_ylabel('Px')
ax7.set_xlabel('Bx(nT)')
ax7.set_title('Linear range of Cosine components of Pz')

for linear_num in range(10):
    slope_s, intercept_s = np.polyfit(B, valzs[linear_num*10,:], 1)
    slope_c, intercept_c = np.polyfit(B[20:80], valzc[linear_num*10,20:80], 1)
    #Plot
    
    ax6.plot(np.unique(B[20:80]), np.poly1d(np.polyfit(B[20:80], valzs[linear_num*10,20:80], 1))(np.unique(B[20:80])),label=f'Bz={round(B2[linear_num*10],2)}nT : Pz={round(slope_s,6)}Bx+{round(intercept_s,7)}')
    ax6.plot(B,valzs[linear_num*10,:],'xb',markersize=1)
    ax6.legend(loc='upper right', fontsize=12)         
    ax7.plot(np.unique(B[20:80]), np.poly1d(np.polyfit(B[20:80], valzc[linear_num*10,20:80], 1))(np.unique(B[20:80])),label= f'Bz={round(B2[linear_num*10],2)}nT : Pz={round(slope_c,6)}Bx+{round(intercept_s,7)}')
    ax7.plot(B,valzc[linear_num*10,:],'xb',markersize=1)
    ax7.legend(loc='upper right', fontsize=12)
    
    plt.show()
#%%
plt.plot(t,ss,t,2*val[:,2])
plt.xlabel('Time(s)',fontsize=14)
plt.ylabel('Signal (Arb)',fontsize=14)