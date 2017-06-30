# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:26:16 2017

@author: Peter Wilson
"""

import numpy as np
import sympy as sp
from sympy import nsolve
print("Sympy version: ",sp.__version__)

x = sp.Symbol('x')
y = sp.Symbol('y')

xc = sp.Symbol('xc')
yc = sp.Symbol('yc')

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

# Section 1 ---------------------------------------------------
# KT disp field
PHI = a1 + a2*x +a3*y + a4*x**2 + 0.5*(a5+a6)*x*y + a7*y**2
PHI = PHI.subs([(x,(x2*y3 - x3*y2 + xc*(y2-y3)+yc*(x3-x2))/detJ),(y,(x3*y1 - x1*y3 + xc*(y3-y1)+yc*(x1-x3))/detJ)])


# RM field
#GAM = a8*x + a9*y - Bub*x*y #generic bubble mode
#GAM = a8*x + a9*y - (a8+a9)*x*y #generic bubble mode
GAM = a8*x + a9*y #no bubble
GAM = GAM.subs([(x,(x2*y3 - x3*y2 + xc*(y2-y3)+yc*(x3-x2))/detJ),(y,(x3*y1 - x1*y3 + xc*(y3-y1)+yc*(x1-x3))/detJ)])

# moderator that reconciles the 2 together
PD = 0.5*(a5-a6)*x*y
PD = PD.subs([(x,(x2*y3 - x3*y2 + xc*(y2-y3)+yc*(x3-x2))/detJ),(y,(x3*y1 - x1*y3 + xc*(y3-y1)+yc*(x1-x3))/detJ)])
#C = a5 - a6

# Section 2 ---------------------------------------------------

PX = PHI + PD
PY = PHI - PD
GX = GAM + PD
GY = GAM - PD


# Section 3 ---------------------------------------------------
px = sp.diff(PX,xc)
py = sp.diff(PY,yc)
gx = sp.diff(GX,xc)
gy = sp.diff(GY,yc)

print("\nSymbolic gx:",gx)
print("\nSymbolic gy:",gy)

# Section 4 ---------------------------------------------------
#PI_eq1 = sp.integrate(gx**2 + gy**2,(y,0.0,1.0-x))
#PI_eq2 = sp.integrate(PI_eq1,(x,0.0,1.0))

#print("\nPI_eq2:",PI_eq2)
#Bub = a8 + a9

# Alternative representation of the displacement field as the difference between thrust and Kirchhoff
W = GAM - PHI
WX = GX - PX
WY = GY - PY

print("W=\t",W)

#Identification of the Ansatz coefficients with the node values for shifts Wi and rotations
# (converted to zero value eqns)
eq1 = W1 - (W.subs(xc,x1)).subs(yc,y1)
eq2 = W2 - (W.subs(xc,x2)).subs(yc,y2)
eq3 = W3 - (W.subs(xc,x3)).subs(yc,y3)

eq4 = PX1 - (px.subs(xc,x1)).subs(yc,y1)
eq5 = PX2 - (px.subs(xc,x2)).subs(yc,y2)
eq6 = PX3 - (px.subs(xc,x3)).subs(yc,y3)

eq7 = PY1 - (py.subs(xc,x1)).subs(yc,y1)
eq8 = PY2 - (py.subs(xc,x2)).subs(yc,y2)
eq9 = PY3 - (py.subs(xc,x3)).subs(yc,y3)

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

eqVec = sp.zeros(9,1)
eqVec[0] = eq1
eqVec[1] = eq2
eqVec[2] = eq3
eqVec[3] = eq4
eqVec[4] = eq5
eqVec[5] = eq6
eqVec[6] = eq7
eqVec[7] = eq8
eqVec[8] = eq9

#results = list(sp.nsolve([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[a1,a2,a3,a4,a5,a6,a7,a8,a9],(W1,1,1,1,1,1,1,1,1)))

results = (sp.solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[a1,a2,a3,a4,a5,a6,a7,a8,a9],warn=True,dict=True))
#ansatzCoefficients = results[0]


ansatzCoefficients = sp.zeros(9,1)
ansatzCoefficients[0] = results[0][a1]
ansatzCoefficients[1] = results[0][a2]
ansatzCoefficients[2] = results[0][a3]
ansatzCoefficients[3] = results[0][a4]
ansatzCoefficients[4] = results[0][a5]
ansatzCoefficients[5] = results[0][a6]
ansatzCoefficients[6] = results[0][a7]
ansatzCoefficients[7] = results[0][a8]
ansatzCoefficients[8] = results[0][a9]

print("\nAnsatz coefficients solved",ansatzCoefficients)
print("\na9=\t",results[0][a9])

#Go through and update everything "_u"
PHI_u = PHI.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

GAM_u = GAM.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

PD_u = PD.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

C_u = C.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
#print("\nUpdated C",C_u)

#PX_u = PHI_u + PD_u
#PY_u = PHI_u - PD_u

# super DSG
GX_u = GAM_u + PD_u
GY_u = GAM_u - PD_u

# EXTRA BIT NOW ==================================================================================
#print("\nGX_u =\t",sp.simplify(GX_u))
#GX_1 = GX_u.subs([(x,x1),(y,y1)])
#print("\nGX_1 =\t",GX_1)
#GX_2 = GX_u.subs([(x,x2),(y,y2)])
#print("\nGX_2 =\t",GX_2)
#GX_3 = GX_u.subs([(x,x3),(y,y3)])
#print("\nGX_3 =\t",GX_3)
#
#GX_sf = (1.0-x-y)*GX_1 + x*GX_2 + y*GX_3
#
#print("\nGX_u - GX_sf =\t",sp.simplify(GX_u-GX_sf))
#
#GY_1 = GY_u.subs([(x,x1),(y,y1)])
#print("\nGY_1 =\t",GY_1)
#GY_2 = GY_u.subs([(x,x2),(y,y2)])
#print("\nGY_2 =\t",GY_2)
#GY_3 = GY_u.subs([(x,x3),(y,y3)])
#print("\nGY_3 =\t",GY_3)

# EXTRA BIT NOW ==================================================================================

# Section 3 ---------------------------------------------------
#px_u = sp.diff(PX_u,x)
#py_u = sp.diff(PY_u,y)
gx_u = sp.diff(GX_u,xc)
gy_u = sp.diff(GY_u,yc)

#print("\nSymbolic GX_u:",GX_u)
#print("\nSymbolic GY_u:",GY_u)

print("\n\n\nSymbolic B Matrix ( parametric space) = ")

# Assemble B Matrix ------------------------------------------
B = sp.zeros(2,9)
for col in range(9):
    B[0,col] = sp.diff(gx_u,UVector[col])
    B[1,col] = sp.diff(gy_u,UVector[col])

sp.pprint(B,wrap_line=False)

disps = sp.zeros(9,1)
disps[0] = 0.5
disps[1] = 1.0
disps[2] = 0.5
disps[3] = 1.0 
disps[4] = -2.0 #solved for
disps[5] = 1.5
disps[6] = 1.0
disps[7] = 1.5 #solved for
disps[8] = -1.0 #solved for
B_test = B.subs([(x,1),(y,0)])
print("\n\n\nShear strains for example problem @ 0,0 = ")
sp.pprint(B_test*disps,wrap_line=False)


# conversion from parametric to cartesian space --------------------------------

# jacobian style taken from p8 and 9 of original DSG formulation
#Inverse jacobian =     1     [dAlpha/dX      dBeta/dX]     1   [ c      -b]
#                     -----   [                       ] = ----- [          ]
#                       2A    [dAlpha/dy      dBeta/dY]    2A   [-d       a]

#a = x2-x1, b=y2-y1, c=y3-y1, d=x3-x1

# d()/dx = d()/dalpha dalpha/dx + d()/dbeta dbeta/dx
# d()/dy = d()/dalpha dalpha/dy + d()/dbeta dbeta/dy

#print("\n\n\nSymbolic B Matrix (cartesian space, factor of 1/detJ taken out) = ")
#for col in range(9):
#    B[0,col] = sp.diff(sp.diff(GX_u,x),UVector[col])*c/detJ + sp.diff(sp.diff(GX_u,y),UVector[col])*-b/detJ
#    B[1,col] = sp.diff(sp.diff(GY_u,x),UVector[col])*-d/detJ + sp.diff(sp.diff(GY_u,y),UVector[col])*a/detJ
#
#B = B.subs([(x,loc1),(y,loc2)])
#sp.pprint(sp.factor(B)*detJ,wrap_line=False) #detJ taken out for clarity


##Rearraging to original DSG dofs for easier comparison
##    Here ---------------->    Original DSG
## [w1,w2,w3,phix1,...]'  -->    [w1,phix1,phiy1,w2,...]'
#B_original_DSG_ordering = sp.zeros(2,9)
#for gamma in range(2):
#    for node in range(3):
#        B_original_DSG_ordering[gamma,node*3] = B[gamma,node]
#        B_original_DSG_ordering[gamma,node*3+1] = B[gamma,3+node]
#        B_original_DSG_ordering[gamma,node*3+2] = B[gamma,6+node]
#B_original_DSG_ordering = B_original_DSG_ordering.subs([(x,loc1),(y,loc2)])
#print("\n\n\nSymbolic B Matrix (cartesian space, factor of 1/detJ taken out, ordered as per original DSG formulation) = ")
#sp.pprint(sp.factor(B_original_DSG_ordering)*detJ,wrap_line=False) #detJ taken out for clarity


x1 = sp.Symbol('x1')
y1 = sp.Symbol('y1')
x2 = sp.Symbol('x2')
y2 = sp.Symbol('y2')
x3 = sp.Symbol('x3')
y3 = sp.Symbol('y3')

print("\n\n\n\n\nPrinting individual entries of matrix above, just for easy copying into C++:")
Bsimp = sp.factor(B) #detJ taken out for clarity
for col in range(9):
    print("BSuper(0,",col,")=",Bsimp[0,col],";")
for col in range(9):
    print("BSuper(1,",col,")=",Bsimp[1,col],";")
    
print("\n\n\n\n\nPrinting individual entries of matrix above, just for easy copying into C++:")
Bsimp = sp.simplify(B.subs([(xc,(1.0-loc1-loc2)*x1+loc1*x2+loc2*x3),(yc,(1.0-loc1-loc2)*y1+loc1*y2+loc2*y3)])) #detJ taken out for clarity
for col in range(9):
    print("BSuper(0,",col,")=",Bsimp[0,col],";")
for col in range(9):
    print("BSuper(1,",col,")=",Bsimp[1,col],";")
    
    