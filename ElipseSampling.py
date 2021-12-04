#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 14:07:05 2021

@author: karan
"""

#Uniformly distributed points inside 3-elispe

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
#equation: ax^2 + by^2 + cz^2 = 1

a, b, c = 0.25, 0.5, 1

n  = 1000

Xp, Yp, Zp = [], [], []

X = np.random.uniform(-1/np.sqrt(a), 1/np.sqrt(a), n) 
Y = np.random.uniform(-1/np.sqrt(b), 1/np.sqrt(b), n) 
Z = np.random.uniform(-1/np.sqrt(c), 1/np.sqrt(c), n) 


for i in range(1000):
    if a*X[i]**2 + b*Y[i]**2 + c*Z[i]**2 < 1:
        Xp.append(X[i])
        Yp.append(Y[i])
        Zp.append(Z[i])
        


fig = plt.figure()
ax = plt.axes(projection='3d')

ax = plt.axes(projection='3d')
#zline = np.linspace(-1.5, 1.5, 1000)
#xline = np.sin(zline)
#yline = np.cos(zline)
#ax.plot3D(xline, yline, zline, 'gray')
#        

ax.scatter3D(Xp, Yp, Zp, c=Zp);
fig.show()
fig.savefig('ellispesample')