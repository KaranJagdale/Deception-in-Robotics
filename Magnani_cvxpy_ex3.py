# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 09:18:00 2021

@author: DELL
"""
#Demonstrates magnani's algorithm to find intersection of zero sublevel sets of polynomials
# -x and x-1

import cvxpy as cp
import numpy as np

nc1 = 2
nc3, nc4, nc5 = 2,2,2
nc2 = 3    #as the degree of constraint equation is 4
nc = 2
C1 = cp.Variable((nc1,nc1), symmetric=True)
C3 = cp.Variable((nc3,nc3), symmetric=True)
C4 = cp.Variable((nc4,nc4), symmetric=True)
C5 = cp.Variable((nc5,nc5), symmetric=True)
C2 = cp.Variable((nc2,nc2), symmetric=True)
C = cp.Variable()
#C3:n, C4:o, C5:p, C1:l, 
#print(C.shape())
Epsilon = 0
Epsilon1 = 1
constraints = [C1 >> Epsilon * np.identity(nc1)]
constraints += [C3 >> Epsilon * np.identity(nc3)]
constraints += [C4 >> Epsilon * np.identity(nc4)]
constraints += [C5 >> Epsilon * np.identity(nc5)]
constraints += [C2 >> Epsilon * np.identity(nc2)]
#constraints += [C >> Epsilon1 * np.identity(nc)]
constraints += [C== 2*C1[1,1]]
constraints += [C2[0,0] == 1-C1[0,0] - C4[0,0] ]
constraints += [C2[1,1] + C2[0,2] + C2[2,0]== -C1[1,1] - 2*C3[0,1] + 2*C4[0,1] - C4[1,1] - C5[0,0] + 2*C5[0,1] ]
constraints += [C2[2,2] == -C5[1,1]]
constraints += [C2[0,1] + C2[1,0] == -2*C1[0,1] - C3[0,0] + C4[0,0] - 2*C4[0,1] + C5[0,1]]
constraints += [C2[1,2] + C2[2,1]  == -C3[1,1] + C4[1,1] - 2*C5[0,1] + C5[1,1]]


#prob = cp.Problem(cp.Minimize(-cp.log(12*C1[2,2]*(4*C1[0,2] + 2*C1[1,1]))), constraints)
#prob = cp.Problem(cp.Minimize(-cp.log_det(C)), constraints)
prob = cp.Problem(cp.Maximize(C), constraints)

prob.solve(verbose = True) #max_iters=100000)
print("status:", prob.status)
print("The optimal value is", prob.value)
Cmat = C.value
# print(np.linalg.det(Cmat))
# print(np.log(np.linalg.det(Cmat)))
# print(cp.log_det(Cmat))
print(Cmat)
print(C1.value)
print(C2.value)
print(C3.value)
print(C4.value)
print(C5.value)



