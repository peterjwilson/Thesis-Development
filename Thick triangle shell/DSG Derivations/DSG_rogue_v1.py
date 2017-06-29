# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:26:16 2017

@author: Peter Wilson
"""

import numpy as np
import sympy as sp
print("Sympy version: ",sp.__version__)

x = sp.Symbol('x')
y = sp.Symbol('y')

loc1 = sp.Symbol('loc1')
loc2 = sp.Symbol('loc2')

alpha = sp.Symbol('alpha')
beta = sp.Symbol('beta')

a1 = sp.Symbol('a1')
a2 = sp.Symbol('a2')
a3 = sp.Symbol('a3')
a4 = sp.Symbol('a4')
a5 = sp.Symbol('a5')
a6 = sp.Symbol('a6')
a7 = sp.Symbol('a7')
a8 = sp.Symbol('a8')
a9 = sp.Symbol('a9')
a10 = sp.Symbol('a10')

PI = sp.Symbol('PI')
C = sp.Symbol('C')

Bub = sp.Symbol('Bub')

W1 = sp.Symbol('W1')
W2 = sp.Symbol('W2')
W3 = sp.Symbol('W3')
PX1 = sp.Symbol('PX1')
PX2 = sp.Symbol('PX2')
PX3 = sp.Symbol('PX3')
PY1 = sp.Symbol('PY1')
PY2 = sp.Symbol('PY2')
PY3 = sp.Symbol('PY3')

a = sp.Symbol('a') #x2 - x1
b = sp.Symbol('b') #y2 - y1
c = sp.Symbol('c') #y3 - y1
d = sp.Symbol('d') #x3 - x1
detJ = sp.Symbol('detJ')
A = sp.Symbol('A') #x2 - x1

#Nodal coords
x1 = sp.Symbol('x1')
y1 = sp.Symbol('y1')
x2 = sp.Symbol('x2')
y2 = sp.Symbol('y2')
x3 = sp.Symbol('x3')
y3 = sp.Symbol('y3')

N4 = 27*x*y*(1.0-x-y)
N1 = (1.0-x-y-N4/3.0)
N2 = (x-N4/3.0)
N3 = (y-N4/3.0)
x_sf = N1*x1 + N2*x2 + N3*x3 + (N4)*(x1 + x2 + x3)/3.0
y_sf = N1*y1 + N2*y2 + N3*y3 + (N4)*(y1 + y2 + y3)/3.0


#Node coordinates set to example values
x1 = 0.0
y1 = 0.0
x2 = 1.0
y2 = 0.0
x3 = 0.0
y3 = 1.0

#Displacement shear volume
W_test = (1.0-x-y)*W1 + x*W2 + y*W3
W_vol1 = sp.integrate(W_test,(y,0.0,1.0-x))
W_vol = sp.integrate(W_vol1,(x,0.0,1.0))



W_av = 2.0 * W_vol / detJ #average displacement gap across element
print("W_av=\n",W_av)
#beta_y_av = (W2-W1)/a #average rotations
#beta_x_av = (W3-W1)/c
#
##Gamma shear volume
gamma_x = -1*W1 + (1.0-x-y)*PX1 + W2 + x*PX2 + y*PX3
gamma_y = -1*W1 + (1.0-x-y)*PX1 + W3 + x*PY2 + y*PY3
gamma_x = -1*W1 + (1.0-x-y)*PX1 + W2 + x*PX2 + y*PX3
gamma_y = -1*W1 + (1.0-x-y)*PX1 + W3 + x*PY2 + y*PY3

sg_n1_gamma_xi = (0) + sp.integrate(gamma_x*a + gamma_y*b,(x,0.0,0.0)) #always 0
sg_n2_gamma_xi = (W2-W1) + sp.integrate(gamma_x*a + gamma_y*b,(x,0.0,1.0))
sg_n3_gamma_xi = (0) + sp.integrate(gamma_x*a + gamma_y*b,(x,0.0,0.0)) #always 0

print("sg_n1_gamma_xi\n",sg_n1_gamma_xi)
print("sg_n2_gamma_xi\n",sg_n2_gamma_xi)
print("sg_n3_gamma_xi\n",sg_n3_gamma_xi)

sg_n1_gamma_eta = (0) + sp.integrate(gamma_x*d + gamma_y*c,(y,0.0,0.0)) #always 0
sg_n2_gamma_eta = (0) + sp.integrate(gamma_x*d + gamma_y*c,(y,0.0,0.0)) #always 0
sg_n3_gamma_eta = (W3-W1) + sp.integrate(gamma_x*d + gamma_y*c,(y,0.0,1.0)) 

#print("sg_n1_gamma_eta\n",sg_n1_gamma_eta)
#print("sg_n2_gamma_eta\n",sg_n2_gamma_eta)
#print("sg_n3_gamma_eta\n",sg_n3_gamma_eta)

#Shear gap fields
sg_gamma_x = x*sg_n2_gamma_xi
sg_gamma_y = y*sg_n3_gamma_eta

# Shear volumes and average shear gaps
sg_gamma_x_vol1 = sp.integrate(sg_gamma_x,(y,0.0,1.0-x))
sg_gamma_x_vol = sp.integrate(sg_gamma_x_vol1,(x,0.0,1.0))
sg_gamma_x_av = 2.0*sg_gamma_x_vol/detJ #average shear gap across whole element

sg_gamma_y_vol1 = sp.integrate(sg_gamma_y,(y,0.0,1.0-x))
sg_gamma_y_vol = sp.integrate(sg_gamma_y_vol1,(x,0.0,1.0))
sg_gamma_y_av = 2.0*sg_gamma_y_vol/detJ #average shear gap across whole element

print("sg_gamma_x_av\n",sg_gamma_x_av)
print("sg_gamma_y_av\n",sg_gamma_y_av)


# Assign average displacement and shear gaps to all nodes with SFs
sg_n1_gamma_xi = sg_gamma_x_av
sg_n2_gamma_xi = sg_gamma_x_av
sg_n2_gamma_xi = sg_gamma_x_av

sg_n1_gamma_eta = sg_gamma_y_av
sg_n2_gamma_eta = sg_gamma_y_av
sg_n3_gamma_eta = sg_gamma_y_av

# Construct shear gap field  with shape functions
sg_xi_sf = (1.0-x-y)*sg_n1_gamma_xi + x*sg_n2_gamma_xi + y*sg_n3_gamma_xi
sg_eta_sf = (1.0-x-y)*sg_n1_gamma_eta + x*sg_n2_gamma_eta + y*sg_n3_gamma_eta

print("sg_xi_sf\n",sg_xi_sf)
print("sg_eta_sf\n",sg_eta_sf)

# Section 3 ---------------------------------------------------
#px_u = sp.diff(PX_u,x)
#py_u = sp.diff(PY_u,y)
gx_u = sp.diff(sg_xi_sf,x)
gy_u = sp.diff(sg_eta_sf,y)

#print("\nSymbolic GX_u:",GX_u)
#print("\nSymbolic GY_u:",GY_u)

UVector = sp.zeros(9,1) #vector of displacements
UVector[0] = W1
UVector[1] = W2
UVector[2] = W3
UVector[3] = PX1
UVector[4] = PX2
UVector[5] = PX3
UVector[6] = PY1
UVector[7] = PY2
UVector[8] = PY3

print("\n\n\nSymbolic B Matrix ( parametric space) = ")

# Assemble B Matrix ------------------------------------------
B = sp.zeros(2,9)
for col in range(9):
    B[0,col] = sp.diff(gx_u,UVector[col])
    B[1,col] = sp.diff(gy_u,UVector[col])

sp.pprint(B,wrap_line=False)



# conversion from parametric to cartesian space --------------------------------

# jacobian style taken from p8 and 9 of original DSG formulation
#Inverse jacobian =     1     [dAlpha/dX      dBeta/dX]     1   [ c      -b]
#                     -----   [                       ] = ----- [          ]
#                       2A    [dAlpha/dy      dBeta/dY]    2A   [-d       a]

#a = x2-x1, b=y2-y1, c=y3-y1, d=x3-x1

# d()/dx = d()/dalpha dalpha/dx + d()/dbeta dbeta/dx
# d()/dy = d()/dalpha dalpha/dy + d()/dbeta dbeta/dy

print("\n\n\nSymbolic B Matrix (cartesian space, factor of 1/detJ taken out) = ")
for col in range(9):
    B[0,col] = sp.diff(sp.diff(sg_xi_sf,x),UVector[col])*c/detJ + sp.diff(sp.diff(sg_xi_sf,y),UVector[col])*-b/detJ
    B[1,col] = sp.diff(sp.diff(sg_eta_sf,x),UVector[col])*-d/detJ + sp.diff(sp.diff(sg_eta_sf,y),UVector[col])*a/detJ

B = B.subs([(x,loc1),(y,loc2)])
sp.pprint(sp.factor(B)*detJ,wrap_line=False) #detJ taken out for clarity


#Rearraging to original DSG dofs for easier comparison
#    Here ---------------->    Original DSG
# [w1,w2,w3,phix1,...]'  -->    [w1,phix1,phiy1,w2,...]'
B_original_DSG_ordering = sp.zeros(2,9)
for gamma in range(2):
    for node in range(3):
        B_original_DSG_ordering[gamma,node*3] = B[gamma,node]
        B_original_DSG_ordering[gamma,node*3+1] = B[gamma,3+node]
        B_original_DSG_ordering[gamma,node*3+2] = B[gamma,6+node]
B_original_DSG_ordering = B_original_DSG_ordering.subs([(x,loc1),(y,loc2)])
print("\n\n\nSymbolic B Matrix (cartesian space, factor of 1/detJ taken out, ordered as per original DSG formulation) = ")
sp.pprint(sp.factor(B_original_DSG_ordering)*detJ,wrap_line=False) #detJ taken out for clarity




#print("\n\n\n\n\nPrinting individual entries of matrix above, just for easy copying into C++:")
#Bsimp = sp.factor(B)*detJ #detJ taken out for clarity
#for col in range(9):
#    print("BSuper(0,",col,")=",Bsimp[0,col],";")
#for col in range(9):
#    print("BSuper(1,",col,")=",Bsimp[1,col],";")
    
    
    
    
    
