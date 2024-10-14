# -*- coding: utf-8 -*-
"""
Created on Sun May 12 11:20:47 2024

@author: jaket
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint

plt.close('all')

gamma = 2 * 7 * np.pi
nsteps = 120
w1 = 10 * gamma
w = 2 * np.pi * 1000
ncycles = 100# Define the larmour frequency gamma
P0 = 1 #starting polarisation
R = 10 * gamma #Relaxation

size = nsteps*ncycles

Bmax = 2
val_step = Bmax/50
B = np.arange(-Bmax, Bmax, val_step)
B2 = np.arange(-Bmax,Bmax, val_step)

dt = 0.000001
   
t = (np.arange(0, ncycles*nsteps, 1)+0.5 )*dt

#%%

pxxval_x = np.load('"C:\Users\jaket\OneDrive\Desktop\Physics Year 4 project\Project Code and figures\data\pxxval_xyz_x_100_in_all.npy"')
pxxval_xy = np.load('')
pxxval_xyz = np.load('')

pyyval_x = np.load('')
pyyval_xy = np.load('')
pyyval_xyz = np.load('')

pzzval_x = np.load('')
pzzval_xy = np.load('')
pzzval_xyz = np.load('')
