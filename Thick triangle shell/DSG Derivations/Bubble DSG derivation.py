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

loc1 = sp.Symbol('loc1')
loc2 = sp.Symbol('loc2')

W1 = sp.Symbol('W1')
W2 = sp.Symbol('W2')
W3 = sp.Symbol('W3')
PX1 = sp.Symbol('PX1')
PX2 = sp.Symbol('PX2')
PX3 = sp.Symbol('PX3')
PY1 = sp.Symbol('PY1')
PY2 = sp.Symbol('PY2')
PY3 = sp.Symbol('PY3')

# THESE ARE KRATOS NODAL ROTATIONS - DEFINED 
RX1 = sp.Symbol('RX1')
RX2 = sp.Symbol('RX2')
RX3 = sp.Symbol('RX3')
RY1 = sp.Symbol('RY1')
RY2 = sp.Symbol('RY2')
RY3 = sp.Symbol('RY3')

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
N4 = 27.0*xi*eta*(1.0-xi-eta)
N1 = (1.0-xi-eta-N4/3.0)
N2 = xi - N4/3.0
N3 = eta -N4/3.0


# Displacement field
W_test = N1*W1 + N2*W2 + N3*W3 + N4*(W1+W2+W3)/3.0

# Rotation field
phi_xi = N1*PX1 + N2*PX2 + N3*PX3 + N4*(PX1+PX2+PX3)/3.0
phi_eta = N1*PY1 + N2*PY2 + N3*PY3+ N4*(PY1+PY2+PY3)/3.0

# Shear gaps at each node
sg_xi_n1 = W_test.subs([(xi,0),(eta,0)])- W_test.subs([(xi,0),(eta,0)]) + sp.integrate(phi_xi*a + phi_eta*b,(xi,0,0)).subs([(xi,0),(eta,0)])

sg_xi_n2 = W_test.subs([(xi,1),(eta,0)])- W_test.subs([(xi,0),(eta,0)]) + sp.integrate(phi_xi*a + phi_eta*b,(xi,0,1)).subs([(xi,1),(eta,0)])

sg_xi_n3 = W_test.subs([(xi,0),(eta,1)])- W_test.subs([(xi,0),(eta,1)]) + sp.integrate(phi_xi*a + phi_eta*b,(xi,0,0)).subs([(xi,0),(eta,1)])

sg_xi_n4 = W_test.subs([(xi,1.0/3.0),(eta,1.0/3.0)])- W_test.subs([(xi,0),(eta,1.0/3.0)]) + sp.integrate(phi_xi*a + phi_eta*b,(xi,0,1.0/3.0)).subs([(xi,1.0/3.0),(eta,1.0/3.0)])

print("\n\n Shear gaps evaluated at each node =:")

print("\n\n sg_xi_n1 =:",sp.simplify(sg_xi_n1))

print("\n\n sg_xi_n2 =:",sp.simplify(sg_xi_n2))

print("\n\n sg_xi_n3 =:",sp.simplify(sg_xi_n3))

print("\n\n sg_xi_n4 =:",sp.simplify(sg_xi_n4))

sg_eta_n1 = W_test.subs([(xi,0),(eta,0)]) - W_test.subs([(xi,0),(eta,0)]) + sp.integrate(phi_xi*d + phi_eta*c,(eta,0,0)).subs([(xi,0),(eta,0)])

sg_eta_n2 = W_test.subs([(xi,1),(eta,0)])- W_test.subs([(xi,1),(eta,0)]) + sp.integrate(phi_xi*d + phi_eta*c,(eta,0,0)).subs([(xi,1),(eta,0)])

sg_eta_n3 = W_test.subs([(xi,0),(eta,1)])- W_test.subs([(xi,0),(eta,0)]) + sp.integrate(phi_xi*d + phi_eta*c,(eta,0,1)).subs([(xi,0),(eta,1)])

sg_eta_n4 = W_test.subs([(xi,1.0/3.0),(eta,1.0/3.0)])- W_test.subs([(xi,1.0/3.0),(eta,0)]) + sp.integrate(phi_xi*d + phi_eta*c,(eta,0,1.0/3.0)).subs([(xi,1.0/3.0),(eta,1.0/3.0)])

print("\n\n sg_eta_n1 =:",sg_eta_n1)

print("\n\n sg_eta_n2 =:",sg_eta_n2)

print("\n\n sg_eta_n3 =:",sg_eta_n3)

print("\n\n sg_eta_n4 =:",sg_eta_n4)

# Assemble shear gap field ----------------------------------------------------
sg_xi_sf = N1*sg_xi_n1 + N2*sg_xi_n2 + N3*sg_xi_n3 + N4*sg_xi_n4
sg_eta_sf = N1*sg_eta_n1 + N2*sg_eta_n2 + N3*sg_eta_n3 + N4*sg_eta_n4
print("\n\n\n\n\n\n")

# swap between plate theory rotations and FEM nodal rotations
sg_xi_sf = sg_xi_sf.subs([(PY1,-RX1),(PY2,-RX2),(PY3,-RX3),(PX1,RY1),(PX2,RY2),(PX3,RY3)])
sg_eta_sf = sg_eta_sf.subs([(PY1,-RX1),(PY2,-RX2),(PY3,-RX3),(PX1,RY1),(PX2,RY2),(PX3,RY3)])

# Shear strain field ----------------------------------------------------------
gx_u = sp.diff(sg_xi_sf,xi)*c/detJ + sp.diff(sg_eta_sf,eta)*-b/detJ
gy_u = sp.diff(sg_eta_sf,eta)*a/detJ + sp.diff(sg_xi_sf,xi)*-d/detJ



# vector of displacements------------------------------------------------------
UVector = sp.zeros(9,1) 
UVector[0] = W1
UVector[1] = W2
UVector[2] = W3
UVector[3] = RX1
UVector[4] = RX2
UVector[5] = RX3
UVector[6] = RY1
UVector[7] = RY2
UVector[8] = RY3

print("Vector of displacements (rotations are per FEM):\n",UVector)

# Setup B matrix --------------------------------------------------------------
B = sp.zeros(2,9)
for col in range(9):
    B[0,col] = sp.diff(gx_u,UVector[col])
    B[1,col] = sp.diff(gy_u,UVector[col])

    

print("\n\n\n\n\nPrinting individual entries of matrix above, just for easy copying into C++:\n")
Bsimp = sp.factor(B)*detJ #detJ taken out for clarity
Bsimp = Bsimp.subs([(xi,loc1),(eta,loc2)])
for col in range(9):
    print("BSuper(0,",col,")=",Bsimp[0,col],";")
for col in range(9):
    print("BSuper(1,",col,")=",Bsimp[1,col],";")
    
    
