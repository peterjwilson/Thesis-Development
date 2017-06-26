# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:26:16 2017

@author: Peter Wilson
"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
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

N4 = 27*x*y*(1.0-x-y)
N1 = (1.0-x-y-N4/3.0)
N2 = (x-N4/3.0)
N3 = (y-N4/3.0)
x_sf = N1*x1 + N2*x2 + N3*x3 + (N4)*(x1 + x2 + x3)/3.0
y_sf = N1*y1 + N2*y2 + N3*y3 + (N4)*(y1 + y2 + y3)/3.0

Jac = sp.zeros(2,2)

Jac[0,0] = sp.simplify(sp.diff(x_sf,x))
Jac[1,0] = sp.simplify(sp.diff(x_sf,y))
Jac[0,1] = sp.simplify(sp.diff(y_sf,x))
Jac[1,1] = sp.simplify(sp.diff(y_sf,y))

sp.pprint(Jac)

print("\ndx_sf/dxi = \t",sp.simplify(sp.diff(x_sf,x)))
print("\ndx_sf/deta = \t",sp.simplify(sp.diff(x_sf,y)))

#print("\ndy_sf/dxi = \t",sp.simplify(sp.diff(y_sf,x)))
#print("\ndy_sf/deta = \t",sp.simplify(sp.diff(y_sf,y)))




print("\nx_sf = \t",x_sf)
print("\nx_sf @ n1= \t",x_sf.subs([(x,0),(y,0)]))
print("\nx_sf @ n2= \t",x_sf.subs([(x,1),(y,0)]))
print("\nx_sf @ n3= \t",x_sf.subs([(x,0),(y,1)]))
print("\nx_sf @ n4= \t",x_sf.subs([(x,1.0/3.0),(y,1.0/3.0)]))


#Node coordinates set to example values
x1 = 0.0
y1 = 0.0
x2 = 1.0
y2 = 0.0
x3 = 0.0
y3 = 1.0

# Section 1 ---------------------------------------------------
# KT disp field
#PHI = a1 + a2*x +a3*y + a4*x**2 + 0.5*(a5+a6)*27*x*y*(1-x-y) + a7*y**2
PHI = a1 + a2*x +a3*y + a4*x**2 + 0.5*(a5+a6)*x*y + a7*y**2

# RM field
#GAM = a8*x + a9*y - Bub*x*y #generic bubble mode
#GAM = a8*x + a9*y - (a8+a9)*x*y #generic bubble mode
GAM = a8*x + a9*y #no bubble
#GAM = a8*x + a9*y - 27*x*y*(1-x-y)*(a8+a9)

# moderator that reconciles the 2 together
PD = 0.5*(a5-a6)*x*y
#PD = 0.5*(a5-a6)*27*x*y*(1-x-y) - 0.5*x*a5 + 0.5*y*a6
#C = a5 - a6

# Section 2 ---------------------------------------------------

PX = PHI + PD
PY = PHI - PD
GX = GAM + PD
GY = GAM - PD


# Section 3 ---------------------------------------------------
px = sp.diff(PX,x)
py = sp.diff(PY,y)
gx = sp.diff(GX,x)
gy = sp.diff(GY,y)

#print("\nSymbolic gx:",gx)
#print("\nSymbolic gy:",gy)

# Section 4 ---------------------------------------------------
PI_eq1 = sp.integrate(gx**2 + gy**2,(y,0.0,1.0-x))
PI_eq2 = sp.integrate(PI_eq1,(x,0.0,1.0))

#print("\nPI_eq2:",PI_eq2)
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
eq5 = PX2 - (px.subs(x,x2)).subs(y,y2) #problem
eq6 = PX3 - (px.subs(x,x3)).subs(y,y3) 

eq7 = PY1 - (py.subs(x,x1)).subs(y,y1)
eq8 = PY2 - (py.subs(x,x2)).subs(y,y2) #problem
eq9 = PY3 - (py.subs(x,x3)).subs(y,y3) 

print("eq1 = \n",eq1)
print("eq2 = \n",eq2)
print("eq3 = \n",eq3)
print("eq4 = \n",eq4)
print("eq5 = \n",eq5)
print("eq6 = \n",eq6)
print("eq7 = \n",eq7)
print("eq8 = \n",eq8)
print("eq9 = \n",eq9)

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

results = (sp.solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[a1,a2,a3,a4,a5,a6,a7,a8,a9],warn=True,dict=True))
#ansatzCoefficients = results[0]

print("results = \n",results)
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
#print("\nAnsatz coefficients solved",ansatzCoefficients)

#Go through and update everything "_u"
PHI_u = PHI.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

GAM_u = GAM.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

PD_u = PD.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

C_u = C.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
#print("\nUpdated C",C_u)

PX_u = PHI_u + PD_u
PY_u = PHI_u - PD_u

# super DSG
GX_u = GAM_u + PD_u
GY_u = GAM_u - PD_u

# EXTRA BIT NOW ==================================================================================
print("\nGX_u =\t",sp.simplify(GX_u))
GX_1 = GX_u.subs([(x,x1),(y,y1)])
print("\nGX_1 =\t",GX_1)
GX_2 = GX_u.subs([(x,x2),(y,y2)])
print("\nGX_2 =\t",GX_2)
GX_3 = GX_u.subs([(x,x3),(y,y3)])
print("\nGX_3 =\t",GX_3)
GX_4 = GX_u.subs([(x,(x1+x2+x3)/3.0),(y,(y1+y2+y3)/3.0)])
print("\nGX_4 =\t",GX_4)

GX_sf = N1*GX_1 + N2*GX_2 + N3*GX_3 + N4*GX_4

print("\nGX_u - GX_sf =\t",sp.simplify(GX_u-GX_sf))

GY_1 = GY_u.subs([(x,x1),(y,y1)])
print("\nGY_1 =\t",GY_1)
GY_2 = GY_u.subs([(x,x2),(y,y2)])
print("\nGY_2 =\t",GY_2)
GY_3 = GY_u.subs([(x,x3),(y,y3)])
print("\nGY_3 =\t",GY_3)
GY_4 = GY_u.subs([(x,(x1+x2+x3)/3.0),(y,(y1+y2+y3)/3.0)])
print("\nGY_4 =\t",GY_4)

GY_sf = N1*GY_1 + N2*GY_2 + N3*GY_3 + N4*GY_4

print("\nGY_u - GY_sf =\t",sp.simplify(GY_u-GY_sf))

GX_u = GX_sf
GY_u = GY_sf




# EXTRA BIT NOW ==================================================================================

# Section 3 ---------------------------------------------------
#px_u = sp.diff(PX_u,x)
#py_u = sp.diff(PY_u,y)
gx_u = sp.diff(GX_u,x)
gy_u = sp.diff(GY_u,y)

print("\nSymbolic GX_u:",sp.simplify(GX_u))
print("\nSymbolic GY_u:",sp.simplify(GY_u))




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
#B_test = B.subs([(x,1),(y,0)])
print("\n\n\nShear strains for example problem @ 0.5,0 = ")
sp.pprint((B.subs([(x,0.5),(y,0)]))*disps,wrap_line=False)

print("\n\n\nShear strains for example problem @ 0,1 = ")
sp.pprint((B.subs([(x,0),(y,1)]))*disps,wrap_line=False)

print("\n\n\nShear strains for example problem @ 0,0 = ")
sp.pprint((B.subs([(x,0),(y,0)]))*disps,wrap_line=False)


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
#sp.pprint(sp.simplify(B)*detJ,wrap_line=False) #detJ taken out for clarity





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




#print("\n\n\n\n\nPrinting individual entries of matrix above, just for easy copying into C++:")
#Bsimp = sp.factor(B)*detJ #detJ taken out for clarity
#for col in range(9):
#    print("BSuper(0,",col,")=",Bsimp[0,col],";")
#for col in range(9):
#    print("BSuper(1,",col,")=",Bsimp[1,col],";")
    
    
    
    
    
C = sp.zeros(2,2)
E = 1000
G = E/2.0 #poisson = 0.0
h = 0.05
C[0,0] = 1
C[1,1] = 1
C *= G*5.0/6.0*h

K = B.T * C * B
K1 = sp.integrate(K,(y,0.0,1.0-x))
Kint = sp.integrate(K1,(x,0.0,1.0))
print("\n\nCartesian stiffness matrix integrated symbolically:")
sp.pprint(sp.N(Kint,3),wrap_line=False)

#setup reduced system and get rid of prescribed displacements
redK = sp.zeros(3,3)
redF = sp.zeros(3,1)
displacedDOFs = [1,2,3,4,6,7]
freeDOFs = [5,8,9]
for DOFrowcheck in range(len(freeDOFs)):
    tempF = 0.0
    row = freeDOFs[DOFrowcheck] - 1
    for DOFcolcheck in range(len(displacedDOFs)):
        col = displacedDOFs[DOFcolcheck] - 1
        tempF += Kint[row,col] * disps[col]
    redF[DOFrowcheck] = 0.0 - tempF

for row in range(len(freeDOFs)):
    for col in range(len(freeDOFs)):
        redK[row,col] = Kint[freeDOFs[row]-1,freeDOFs[col]-1]

solvedDisps = (redK**-1)*redF
print("\nSolved displacements:")
sp.pprint(solvedDisps)

BletzyResults = sp.zeros(3,1)
BletzyResults[0] = -2.0
BletzyResults[1] = 1.5
BletzyResults[2] = -1.0
print("\nBletzinger's displacements:")
sp.pprint(BletzyResults)
#
#print("\n\n\n=====================================================================\n")
#print("\tNOW RE-RUNNING PROBLEM INCLUDING BENDING FORMULATION")
#print("\n=====================================================================\n")
#
##Bletzinger maple worksheet DOF order maintained here
#
#B_bend_and_shear = sp.zeros(5,9)
#
#a = x2 - x1
#b = y2 - y1
#c = y3 - y1
#d = x3 - x1
#
##Add in bending B matrix from DSG, re-arranged for maple DOF ordering
#
##phix's
#B_bend_and_shear[0,3] = b-c 
#B_bend_and_shear[1,3] = 0.0
#B_bend_and_shear[2,3] = d-a
#
#B_bend_and_shear[0,4] = c
#B_bend_and_shear[1,4] = 0.0
#B_bend_and_shear[2,4] = -d
#
#B_bend_and_shear[0,5] = -b
#B_bend_and_shear[1,5] = 0.0
#B_bend_and_shear[2,5] = a
#
##phiy's
#B_bend_and_shear[0,6] = 0.0
#B_bend_and_shear[1,6] = d-a
#B_bend_and_shear[2,6] = b-c
#
#B_bend_and_shear[0,7] = 0.0
#B_bend_and_shear[1,7] = -d
#B_bend_and_shear[2,7] = c
#
#B_bend_and_shear[0,8] = 0.0
#B_bend_and_shear[1,8] = a
#B_bend_and_shear[2,8] = -b
#
#area = 1.0*1.0/2.0
#B_bend_and_shear /= (2.0*area) #divide bending part of Bmat by detJ
#
##Copy in shear B matrix, already calculated above, no scaling needed
#for row in range(2):
#    for col in range(9):
#        B_bend_and_shear[3+row,col] = B[row,col]
#
## Setup combined material matrix
#C_bend_and_shear = sp.zeros(5,5)
#C_bend_and_shear[0,0] = 1.0
#C_bend_and_shear[1,1] = 1.0
#C_bend_and_shear[0,1] = 0.0
#C_bend_and_shear[1,0] = 0.0
#C_bend_and_shear[2,2] = (1.0-0.0)/2.0
#C_bend_and_shear *= (E*h**3/(12))
#for row in range(2):
#    for col in range(2):
#        C_bend_and_shear[3+row,3+col] = C[row,col]
#
##Calc combined K
#K_temp3 = B_bend_and_shear.T * C_bend_and_shear * B_bend_and_shear
#Ktemp1 = sp.integrate(K_temp3,(y,0.0,1.0-x))
#K_bend_and_shear = sp.integrate(Ktemp1,(x,0.0,1.0))
#print("\n\nCombined stiffness matrix :")
#sp.pprint(sp.N(K_bend_and_shear,3),wrap_line=False)
#
##setup reduced system and get rid of prescribed displacements
#for DOFrowcheck in range(len(freeDOFs)):
#    tempF = 0.0
#    row = freeDOFs[DOFrowcheck] - 1
#    for DOFcolcheck in range(len(displacedDOFs)):
#        col = displacedDOFs[DOFcolcheck] - 1
#        tempF += K_bend_and_shear[row,col] * disps[col]
#    redF[DOFrowcheck] = 0.0 - tempF
#
#for row in range(len(freeDOFs)):
#    for col in range(len(freeDOFs)):
#        redK[row,col] = K_bend_and_shear[freeDOFs[row]-1,freeDOFs[col]-1]
#
#solvedDisps = (redK**-1)*redF
#print("\nSolved displacements:")
#sp.pprint(solvedDisps)