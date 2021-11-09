# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 09:18:00 2021

@author: Karan
"""

import cvxpy as cp
import numpy as np

nc1 = 3
nc3, nc4, nc5 = 2,2,2
nc2 = 4
nc = 2
C1 = cp.Variable((nc1,nc1), symmetric=True)
C3 = cp.Variable((nc3,nc3), symmetric=True)
C4 = cp.Variable((nc4,nc4), symmetric=True)
C5 = cp.Variable((nc5,nc5), symmetric=True)
C2 = cp.Variable((nc2,nc2), symmetric=True)
C = cp.Variable((nc,nc))

#print(C.shape())
Epsilon = 0
constraints = [C1 >> Epsilon * np.identity(nc1)]
constraints += [C3 >> Epsilon * np.identity(nc3)]
constraints += [C4 >> Epsilon * np.identity(nc4)]
constraints += [C5 >> Epsilon * np.identity(nc5)]
constraints += [C2 >> Epsilon * np.identity(nc2)]
constraints += [C >> Epsilon * np.identity(nc)]
constraints += [C[0,0] == 4*C1[0,2] + 2*C1[1,1], C[0,1] == 12*C1[1,2], C[1,1] == 12*C1[2,2], C[1,0] == 0]
constraints += [C2[0,0] == 1-C1[0,0] + 2.9*C3[0,0] - 0.9*C[0,0] + 2.61*C5[0,0]]
constraints += [C2[1,1] + C2[0,2] == -2*C1[0,2] - C1[1,1] + 6*C3[0,1] + 2.9*C3[1,1] + C4[0,0] - 0.9*C4[1,1] - 2.9*C5[0,0] + 5.4*C5[0,1] + 2.61*C5[1,1]]
constraints += [C2[2,2] + C2[1,3] == -C1[2,2] + C4[1,1] - 6*C5[0,1] - 2.9*C5[1,1]]
constraints += [C2[0,1] == -2*C1[0,1] + 3*C3[0,0] + 5.8*C3[0,1] - 1.8*C4[0,1] + 2.7*C5[0,0] + 5.22*C5[0,1]]
constraints += [C2[1,2] + C2[0,3] == -2*C1[1,2] + 3*C3[1,1] + 2*C4[0,1] - 3*C5[0,0] - 5.8*C5[0,1] + 2.7*C5[1,1]]
constraints += [C2[2,3] == -3*C5[1,1]]
constraints += [C2[3,3] == 0]

#prob = cp.Problem(cp.Minimize(-cp.log(12*C1[2,2]*(4*C1[0,2] + 2*C1[1,1]))), constraints)
#prob = cp.Problem(cp.Minimize(-cp.log_det(C)), constraints)
prob = cp.Problem(cp.Minimize(1), constraints)

prob.solve(verbose = True)
print("status:", prob.status)
print("The optimal value is", prob.value)
print("g is : ", C.value)






