# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 10:26:16 2017

@author: Peter Wilson
"""

import numpy as np
import sympy as sp
print("Sympy version: ",sp.__version__)

x = sp.Symbol('x')
y = sp.Symbol('y')

W1 = sp.Symbol('W1')
W2 = sp.Symbol('W2')
W3 = sp.Symbol('W3')
PX1 = sp.Symbol('PX1')
PX2 = sp.Symbol('PX2')
PX3 = sp.Symbol('PX3')
PY1 = sp.Symbol('PY1')
PY2 = sp.Symbol('PY2')
PY3 = sp.Symbol('PY3')

#Nodal coords
x1 = sp.Symbol('x1')
y1 = sp.Symbol('y1')
x2 = sp.Symbol('x2')
y2 = sp.Symbol('y2')
x3 = sp.Symbol('x3')
y3 = sp.Symbol('y3')

xi = sp.Symbol('xi')
eta = sp.Symbol('eta')

xi1 = sp.Symbol('xi1')
xi2 = sp.Symbol('xi2')
xi3 = sp.Symbol('xi3')
eta1 = sp.Symbol('eta1')
eta2 = sp.Symbol('eta2')
eta3 = sp.Symbol('eta3')

a = sp.Symbol('a') #x2 - x1
b = sp.Symbol('b') #y2 - y1
c = sp.Symbol('c') #y3 - y1
d = sp.Symbol('d') #x3 - x1
detJ = sp.Symbol('detJ')
A = sp.Symbol('A')

#Shape functions
N1 = (1.0-xi-eta)
N2 = xi
N3 = eta


def evalGammaXi(expression,lowerLimit,upperLimit):
    return expression.subs([(xi,upperLimit)]) - expression.subs([(xi,lowerLimit)])

def evalGammaEta(expression,lowerLimit,upperLimit):
    return expression.subs([(eta,upperLimit)]) - expression.subs([(eta,lowerLimit)])

# Displacement field
W_test = N1*W1 + N2*W2 + N3*W3

# Rotation field
phi_xi = N1*PX1 + N2*PX2 + N3*PX3 
phi_eta = N1*PY1 + N2*PY2 + N3*PY3 

# General shear gap formulae
sg_gamma1 = W_test + a*((xi-xi*xi/2.0-xi*eta)*PX1 + 0.5*xi*xi*PX2 + xi*eta*PX3) + b*((xi-xi*xi/2.0-xi*eta)*PY1 + 0.5*xi*xi*PY2 + xi*eta*PY3)

sg_gamma2 = W_test + d*((eta-eta*eta/2.0-xi*eta)*PX1 + xi*eta*PX2 + 0.5*eta*eta*PX3) + c*((eta-eta*eta/2.0-xi*eta)*PY1 + xi*eta*PY2 + 0.5*eta*eta*PY3)

# Sub xi and eta as a,b
#sg_gamma1_n1 = evalGammaXi(sg_gamma1,0,0)
#sg_gamma1_n2 = evalGammaXi(sg_gamma1,0,a)
#sg_gamma1_n3 = evalGammaXi(sg_gamma1,0,d)
#sg_gamma2_n1 = evalGammaEta(sg_gamma1,0,0)
#sg_gamma2_n2 = evalGammaEta(sg_gamma1,0,b)
#sg_gamma2_n3 = evalGammaEta(sg_gamma1,0,c)

# Sub xi and eta as parametric 0's and 1's
sg_gamma1_n1 = evalGammaXi(sg_gamma1,0,0)
sg_gamma1_n2 = evalGammaXi(sg_gamma1,0,1)
sg_gamma1_n3 = evalGammaXi(sg_gamma1,0,0)
sg_gamma2_n1 = evalGammaEta(sg_gamma1,0,0)
sg_gamma2_n2 = evalGammaEta(sg_gamma1,0,0)
sg_gamma2_n3 = evalGammaEta(sg_gamma1,0,1)

print("\n\n Shear gaps evaluated at each node =:")

print("\n\n sg_xi_n1 =:",sp.simplify(sg_gamma1_n1.subs([(xi,0),(eta,0)])))

print("\n\n sg_xi_n2 =:",sp.simplify(sg_gamma1_n2.subs([(xi,1),(eta,0)])))

print("\n\n sg_xi_n3 =:",sp.simplify(sg_gamma1_n3.subs([(xi,0),(eta,1)])))



print("\n\n sg_eta_n1 =:",sp.simplify(sg_gamma2_n1.subs([(xi,0),(eta,0)])))

print("\n\n sg_eta_n2 =:",sp.simplify(sg_gamma2_n2.subs([(xi,1),(eta,0)])))

print("\n\n sg_eta_n3 =:",sp.simplify(sg_gamma2_n3.subs([(xi,0),(eta,1)])))



# Assemble shear gap field ----------------------------------------------------
sg_xi_sf = N1*sg_gamma1_n1 + N2*sg_gamma1_n2 + N3*sg_gamma1_n3
sg_eta_sf = N1*sg_gamma2_n1 + N2*sg_gamma2_n2 + N3*sg_gamma2_n3
print("\n\n\n\n\n\n")

# Shear strain field ----------------------------------------------------------
gx_u = sp.diff(sg_xi_sf,xi)*c/detJ + sp.diff(sg_xi_sf,eta)*-b/detJ
gy_u = sp.diff(sg_eta_sf,eta)*a/detJ + sp.diff(sg_eta_sf,xi)*-d/detJ

# vector of displacements------------------------------------------------------
UVector = sp.zeros(9,1) 
UVector[0] = W1
UVector[1] = W2
UVector[2] = W3
UVector[3] = PX1
UVector[4] = PX2
UVector[5] = PX3
UVector[6] = PY1
UVector[7] = PY2
UVector[8] = PY3

# Setup B matrix --------------------------------------------------------------
B = sp.zeros(2,9)
for col in range(9):
    B[0,col] = sp.diff(gx_u,UVector[col])
    B[1,col] = sp.diff(gy_u,UVector[col])

#Rearraging B-matrix to original DSG dofs for easier comparison ---------------
#    Here ----------------->    Original DSG
# [w1,w2,w3,phix1,...]'  -->    [w1,phix1,phiy1,w2,...]'
B_original_DSG_ordering = sp.zeros(2,9)
for gamma in range(2):
    for node in range(3):
        B_original_DSG_ordering[gamma,node*3] = B[gamma,node]
        B_original_DSG_ordering[gamma,node*3+1] = B[gamma,3+node]
        B_original_DSG_ordering[gamma,node*3+2] = B[gamma,6+node]

print("\n\n\nB Matrix (cartesian space, factor of 1/detJ taken out, ordered as per original DSG formulation) = \n")
sp.pprint(sp.factor(B_original_DSG_ordering)*detJ,wrap_line=False) #detJ taken out for clarity





print("\n\n\n\n\nPrinting individual entries of matrix above, just for easy copying into C++:\n")
Bsimp = sp.factor(B_original_DSG_ordering)*detJ #detJ taken out for clarity
for col in range(9):
    print("BSuper(0,",col,")=",Bsimp[0,col],";")
for col in range(9):
    print("BSuper(1,",col,")=",Bsimp[1,col],";")
    
    
