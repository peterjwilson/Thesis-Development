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

xi = sp.Symbol('xi')
eta = sp.Symbol('eta')

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


# THESE ARE KRATOS NODAL ROTATIONS - DEFINED 
RX1 = sp.Symbol('RX1')
RX2 = sp.Symbol('RX2')
RX3 = sp.Symbol('RX3')
RY1 = sp.Symbol('RY1')
RY2 = sp.Symbol('RY2')
RY3 = sp.Symbol('RY3')

a = sp.Symbol('a') #x2 - x1
b = sp.Symbol('b') #y2 - y1
c = sp.Symbol('c') #y3 - y1
d = sp.Symbol('d') #x3 - x1
detJ = sp.Symbol('detJ')


#Nodal coords
x1 = sp.Symbol('x1')
y1 = sp.Symbol('y1')
x2 = sp.Symbol('x2')
y2 = sp.Symbol('y2')
x3 = sp.Symbol('x3')
y3 = sp.Symbol('y3')


#Node coordinates set to example values
x1 = 0.0
y1 = 0.0
x2 = 1.0
y2 = 0.0
x3 = 0.0
y3 = 1.0

#Shape functions
N1 = (1.0-xi-eta)
N2 = xi
N3 = eta

# Section 1 ---------------------------------------------------
PHI = a1 + a2*x +a3*y + a4*x**2 + 0.5*(a5+a6)*x*y + a7*y**2

#GAM = a8*x + a9*y - Bub*x*y #generic bubble mode
GAM = a8*x + a9*y - (a8+a9)*x*y #generic bubble mode
#GAM = a8*x + a9*y #no bubble

PD = 0.5*(a5-a6)*x*y
C = a5 - a6

# Section 2 ---------------------------------------------------
is_normal_DSG = False

PX = PHI + PD
PY = PHI - PD
if(is_normal_DSG):
    print("\nUsing basic DSG formulation!")
    GX = a8*x + a9*y #normal DSG
    GY = a8*x + a9*y #normal DSG
else:
    GX = GAM + PD
    GY = GAM - PD


# Section 3 ---------------------------------------------------
px = sp.diff(PX,x)
py = sp.diff(PY,y)
gx = sp.diff(GX,x)
gy = sp.diff(GY,y)

print("\nSymbolic gx:",gx)
print("\nSymbolic gy:",gy)

# Section 4 ---------------------------------------------------
# Skip internal energy derivation, just write bubble mode result
Bub = a8 + a9

# Alternative representation of the displacement field as the difference between thrust and Kirchhoff
W = GAM - PHI
WX = GX - PX
WY = GY - PY

#Identification of the Ansatz coefficients with the node values for shifts Wi and rotations
# (converted to zero value eqns)
eq1 = W1 - (W.subs(x,x1)).subs(y,y1)
eq2 = W2 - (W.subs(x,x2)).subs(y,y2)
eq3 = W3 - (W.subs(x,x3)).subs(y,y3)

eq4 = PX1 - (px.subs(x,x1)).subs(y,y1)
eq5 = PX2 - (px.subs(x,x2)).subs(y,y2)
eq6 = PX3 - (px.subs(x,x3)).subs(y,y3)

eq7 = PY1 - (py.subs(x,x1)).subs(y,y1)
eq8 = PY2 - (py.subs(x,x2)).subs(y,y2)
eq9 = PY3 - (py.subs(x,x3)).subs(y,y3)

# Setup system to solve [A][a] = [W] -------------------------
A = sp.zeros(9) #system matrix
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

results = list(sp.linsolve([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[a1,a2,a3,a4,a5,a6,a7,a8,a9]))
ansatzCoefficients = results[0]
print("\nAnsatz coefficients solved")
for i in range(9):
    print("a",i,"=\t",ansatzCoefficients[i])

#Go through and update everything "_u"
PHI_u = PHI.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

GAM_u = GAM.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

PD_u = PD.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

C_u = C.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
print("\nUpdated C",C_u)

#PX_u = PHI_u + PD_u
#PY_u = PHI_u - PD_u

if (is_normal_DSG):
    #normal DSG
    GX_u = GX.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
    GY_u = GY.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
else:
    # super DSG
    GX_u = GAM_u + PD_u
    GY_u = GAM_u - PD_u

print("\nGAM_u_n1:",sp.simplify(GAM_u.subs([(x,0),(y,0)])))
print("\nGAM_u_n2:",sp.simplify(GAM_u.subs([(x,1),(y,0)])))
print("\nGAM_u_n3:",sp.simplify(GAM_u.subs([(x,0),(y,1)])))

# Reconstruction of above field using nodal shear gaps and SFs ###########################
    
# 'Shear gaps' evaluated at nodes
print("\n\nShear gaps evaluated at nodes:")
print("\nGX_u_n1:",sp.simplify(GX_u.subs([(x,0),(y,0)])))
print("\nGX_u_n2:",sp.simplify(GX_u.subs([(x,1),(y,0)])))
print("\nGX_u_n3:",sp.simplify(GX_u.subs([(x,0),(y,1)])))

print("\n\nGY_u_n1:",sp.simplify(GY_u.subs([(x,0),(y,0)])))
print("\n\nGY_u_n2:",sp.simplify(GY_u.subs([(x,1),(y,0)])))
print("\n\nGY_u_n3:",sp.simplify(GY_u.subs([(x,0),(y,1)])))

sg_xi_n1 = GX_u.subs([(x,0),(y,0)])
sg_xi_n2 = GX_u.subs([(x,1),(y,0)])
sg_xi_n3 = GX_u.subs([(x,0),(y,1)])

sg_eta_n1 = GY_u.subs([(x,0),(y,0)])
sg_eta_n2 = GY_u.subs([(x,1),(y,0)])
sg_eta_n3 = GY_u.subs([(x,0),(y,1)])

# Manually override shear gaps, cartesian differences a thru d introduced
print("\n\nMANUALLY OVERRIDING SHEAR GAPS!!!!!!!!!!")
sg_xi_n1 = 0.0
sg_xi_n2 = 0.5*a*PX1 + 0.5*a*PX2 - 1.0*W1 + 1.0*W2
sg_xi_n3 = 0.5*b*PY1 + 0.5*b*PY3 - 1.0*W1 + 1.0*W3

sg_eta_n1 = 0.0
sg_eta_n2 = 0.5*d*PX1 + 0.5*d*PX2 - 1.0*W1 + 1.0*W2
sg_eta_n3 = 0.5*c*PY1 + 0.5*c*PY3 - 1.0*W1 + 1.0*W3

sg_xi_sf = (1.0-x-y)*sg_xi_n1 + x*sg_xi_n2 + y*sg_xi_n3
sg_eta_sf = (1.0-x-y)*sg_eta_n1 + x*sg_eta_n2 + y*sg_eta_n3

GX_u = sg_xi_sf
GY_u = sg_eta_sf

# Reconstruction of above field using nodal shear gaps and SFs ###########################

# DOF transformation from plate theory to FEM rotational dofs
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

GX_u = GX_u.subs([(PY1,-RX1),(PY2,-RX2),(PY3,-RX3),(PX1,RY1),(PX2,RY2),(PX3,RY3)])
GY_u = GY_u.subs([(PY1,-RX1),(PY2,-RX2),(PY3,-RX3),(PX1,RY1),(PX2,RY2),(PX3,RY3)])


# Section 3 ---------------------------------------------------
#px_u = sp.diff(PX_u,x)
#py_u = sp.diff(PY_u,y)
#gx_u = sp.diff(GX_u,x)
#gy_u = sp.diff(GY_u,y)

# Cartesian transformation
gx_u = sp.diff(GX_u,x)*c/detJ + sp.diff(GY_u,y)*-b/detJ
gy_u = sp.diff(GY_u,y)*a/detJ + sp.diff(GX_u,x)*-d/detJ

print("\n\n\nSymbolic B Matrix ( transformed to cartesian ) (1/detJ taken out) = ")

# Assemble B Matrix ------------------------------------------
B = sp.zeros(2,9)
for col in range(9):
    B[0,col] = sp.diff(gx_u,UVector[col])
    B[1,col] = sp.diff(gy_u,UVector[col])

sp.pprint(sp.simplify(B*detJ),wrap_line=False)


print("\n\n\n\n\nPrinting individual entries of original Bmat (1/detJ taken out), just for easy copying into C++:")
#Bsimp = sp.factor(B)*detJ #detJ taken out for clarity
#testing below
B = B.subs([(x,loc1),(y,loc2)])
Bsimp = sp.factor(B)*detJ
#testing end
for col in range(9):
    print("BSuper(0,",col,")=",Bsimp[0,col],";")
for col in range(9):
    print("BSuper(1,",col,")=",Bsimp[1,col],";")

    
    
    
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