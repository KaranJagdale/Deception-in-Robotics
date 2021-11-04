# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 10:52:51 2021

@author: DELL
"""

from sympy import *
from sympy import MatrixSymbol, Matrix, Transpose
from SumOfSquares import SOSProblem, poly_opt_prob
import math
import numpy as np
from numpy import linalg as lin

b11, b22, b33, b12, b13, b23 = symbols('b11, b22, b33, b12, b13, b23')
b = Matrix([[b11, b12, b13], [b12, b22, b23], [b13, b23, b33]])

print(b)

l11, l12, l13, l12, l22, l23, l13, l23, l33 = symbols('l11 l12 l13 l12 l22 l23 l13 l23 l33')

m11, m12, m13, m14, m12, m22, m23, m24, m13, m23, m33, m34, m14, m24, m34, m44 = symbols('m11 m12 m13 m14 m12 m22 m23 m24 m13 m23 m33 m34 m14 m24 m34 m44')

n11, n12, n12, n22 = symbols('n11 n12 n12 n22')

o11, o12, o12, o22 = symbols('o11 o12 o12 o22')

p11, p12, p12, p22 = symbols('p11 p12 p12 p22')

x = symbols('x')

h1 = Matrix([1, x, x**2])

C1 = Matrix([[l11, l12, l13], [l12, l22, l23], [l13, l23, l33]])

C2 = Matrix([[m11, m12, m13, m14], [m12, m22, m23, m24], [m13, m23, m33, m34], [m14, m24, m34, m44]])

C3 = Matrix([[n11, n12], [n12, n22]])

C4 = Matrix([[o11, o12], [o12, o22]])

C5 = Matrix([[p11, p12], [p12, p22]])

h2 = Matrix([1, x, x**2, x**3])

h3 = Matrix([1, x])

h4, h5 = h3, h3

C = Matrix([[4*l13+2*l22, 12*l23],[0, 12*l33]])

# we will define below the polynomial which are constrained to be SOS
g = expand(Transpose(h1)*C1*h1)
w1 = expand(Transpose(h3)*C3*h3)
w2 = expand(Transpose(h4)*C4*h4)
w3 = expand(Transpose(h5)*C5*h5)

f1 = x+7/4
f2 = 4*x + 7
print('Karan')
print(g,w1[0],w2[0],w3[0])

p = expand((1-g[0]) + w1[0]*f1 + w2[0]*f2 - w3[0]*f1*f2)

#print(p)

#Defining the SOS constraints
prob = SOSProblem()
prob.add_sos_constraint(p, [x])
prob.add_sos_constraint(w1[0], [x])
prob.add_sos_constraint(g[0], [x])

prob.add_sos_constraint(w2[0], [x])
prob.add_sos_constraint(w3[0], [x])

#C1v, C2v, C3v, C4v, C5v = prob.sym_to_var(C1), prob.sym_to_var(C2), prob.sym_to_var(C3), prob.sym_to_var(C4), prob.sym_to_var(C5)
l11, l12, l13, l12, l22, l23, l13, l23, l33 = prob.sym_to_var(l11), prob.sym_to_var(l12), prob.sym_to_var(l13), prob.sym_to_var(l12), prob.sym_to_var(l22), prob.sym_to_var(l23), prob.sym_to_var(l13), prob.sym_to_var(l23), prob.sym_to_var(l33)
m11, m12, m13, m14, m12, m22, m23, m24, m13, m23, m33, m34, m14, m24, m34, m44 = prob.sym_to_var(m11), prob.sym_to_var(m12), prob.sym_to_var(m13), prob.sym_to_var(m14), prob.sym_to_var(m12), prob.sym_to_var(m22), prob.sym_to_var(m23), prob.sym_to_var(m24), prob.sym_to_var(m13), prob.sym_to_var(m23), prob.sym_to_var(m33), prob.sym_to_var(m34), prob.sym_to_var(m14), prob.sym_to_var(m24), prob.sym_to_var(m34), prob.sym_to_var(m44)
n11, n12, n12, n22 = prob.sym_to_var(n11), prob.sym_to_var(n12), prob.sym_to_var(n12), prob.sym_to_var(n22)

o11, o12, o12, o22 = prob.sym_to_var(o11), prob.sym_to_var(o12), prob.sym_to_var(o12), prob.sym_to_var(o22)

p11, p12, p12, p22 = prob.sym_to_var(p11), prob.sym_to_var(p12), prob.sym_to_var(p12), prob.sym_to_var(p22)

prob.set_objective('max', 12*l33*(4*l13+2*l22))
prob.solve()




#print(simplify(Transpose(h1)*C1*h1))
