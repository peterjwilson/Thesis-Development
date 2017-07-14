# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:26:16 2017

@author: Peter Wilson
"""

import numpy as np
import sympy as sp
print("sympy version: ",sp.__version__)

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
kirchhoffDispField = a1 + a2*xi +a3*eta + a4*xi**2 + 0.5*(a5+a6)*xi*eta + a7*eta**2

#rmDispField = a8*xi + a9*eta - Bub*xi*eta #generic bubble mode
#rmDispField = a8*xi + a9*eta - (a8+a9)*xi*eta #generic bubble mode
rmDispField = a8*xi + a9*eta #no bubble

moderatorDispField = 0.5*(a5-a6)*xi*eta
C = a5 - a6

# Section 2 ---------------------------------------------------
is_normal_DSG = False

PX = kirchhoffDispField + moderatorDispField
PY = kirchhoffDispField - moderatorDispField
if(is_normal_DSG):
    print("\nUsing basic DSG formulation!")
    GX = a8*xi + a9*eta #normal DSG
    GY = a8*xi + a9*eta #normal DSG
else:
    GX = rmDispField + moderatorDispField
    GY = rmDispField - moderatorDispField


# Section 3 ---------------------------------------------------
kirchhoffShearField_X = sp.diff(PX,xi)
kirchhoffShearField_Y = sp.diff(PY,eta)
rmShearField_X = sp.diff(GX,xi)
rmShearField_Y = sp.diff(GY,eta)

print("\nSymbolic rmShearField_X:",rmShearField_X)
print("\nSymbolic rmShearField_Y:",rmShearField_Y)

# Section 4 ---------------------------------------------------
# Skip internal enerrmShearField_Y derivation, just write bubble mode result
Bub = a8 + a9

# Alternative representation of the shear gap field
shearGapField = rmDispField - kirchhoffDispField
shearGapField_X = GX - PX
shearGapField_Y = GY - PY

#Identification of the Ansatz coefficients with the node values for shifts Wi and rotations
# (converted to zero value eqns)
eq1 = W1 - (shearGapField.subs(xi,x1)).subs(eta,y1)
eq2 = W2 - (shearGapField.subs(xi,x2)).subs(eta,y2)
eq3 = W3 - (shearGapField.subs(xi,x3)).subs(eta,y3)

eq4 = PX1 - (kirchhoffShearField_X.subs(xi,x1)).subs(eta,y1)
eq5 = PX2 - (kirchhoffShearField_X.subs(xi,x2)).subs(eta,y2)
eq6 = PX3 - (kirchhoffShearField_X.subs(xi,x3)).subs(eta,y3)

eq7 = PY1 - (kirchhoffShearField_Y.subs(xi,x1)).subs(eta,y1)
eq8 = PY2 - (kirchhoffShearField_Y.subs(xi,x2)).subs(eta,y2)
eq9 = PY3 - (kirchhoffShearField_Y.subs(xi,x3)).subs(eta,y3)

print("eq1 = ",eq1)
print("eq2 = ",eq2)
print("eq3 = ",eq3)
print("eq4 = ",eq4)
print("eq5 = ",eq5)
print("eq6 = ",eq6)
print("eq7 = ",eq7)
print("eq8 = ",eq8)
print("eq9 = ",eq9)

# Setup system to solve [A][a] = [shearGapField] -------------------------
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
kirchhoffDispField_u = kirchhoffDispField.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

rmDispField_u = rmDispField.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

print("\nrmDispField_u = ",rmDispField_u)

moderatorDispField_u = moderatorDispField.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

print("\nmoderatorDispField_u = ",moderatorDispField_u)

C_u = C.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

#PX_u = kirchhoffDispField_u + moderatorDispField_u
#PY_u = kirchhoffDispField_u - moderatorDispField_u


# Sub in Kratos rotations
rmDispField_u = rmDispField_u.subs([(PY1,-RX1),(PY2,-RX2),(PY3,-RX3),(PX1,RY1),(PX2,RY2),(PX3,RY3)])
moderatorDispField_u = moderatorDispField_u.subs([(PY1,-RX1),(PY2,-RX2),(PY3,-RX3),(PX1,RY1),(PX2,RY2),(PX3,RY3)])



if (is_normal_DSG):
    #normal DSG
    GX_u = GX.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
    GY_u = GY.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
else:
    # super DSG
    GX_u = rmDispField_u + moderatorDispField_u
    GY_u = rmDispField_u - moderatorDispField_u

    
    

print("\nprelim rmShearField_X_u:",sp.diff(GX_u,xi))
print("\nprelim rmShearField_Y_u:",sp.diff(GY_u,eta))




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
#kirchhoffShearField_X_u = sp.diff(PX_u,xi)
#kirchhoffShearField_Y_u = sp.diff(PY_u,eta)
#rmShearField_X_u = sp.diff(GX_u,xi)
#rmShearField_Y_u = sp.diff(GY_u,eta)

# Cartesian transformation
rmShearField_X_u = sp.diff(GX_u,xi)*c/detJ + sp.diff(GY_u,eta)*-b/detJ
rmShearField_Y_u = sp.diff(GY_u,eta)*a/detJ + sp.diff(GX_u,xi)*-d/detJ

print("\n\n\nSymbolic B Matrix ( transformed to cartesian ) (1/detJ taken out) = ")

# Assemble B Matrix ------------------------------------------
B = sp.zeros(2,9)
for col in range(9):
    B[0,col] = sp.diff(rmShearField_X_u,UVector[col])
    B[1,col] = sp.diff(rmShearField_Y_u,UVector[col])

sp.pprint(sp.simplify(B*detJ),wrap_line=False)


print("\n\n\n\n\nPrinting individual entries of original Bmat (1/detJ taken out), just for easy cokirchhoffShearField_Ying into C++:")
#Bsimp = sp.factor(B)*detJ #detJ taken out for clarity
#testing below
B = B.subs([(xi,loc1),(eta,loc2)])
Bsimp = sp.factor(B)*detJ
#testing end
for col in range(9):
    print("BSuper(0,",col,")=",Bsimp[0,col],";")
for col in range(9):
    print("BSuper(1,",col,")=",Bsimp[1,col],";")

    
    
    
#Rearraging B-matrix to original DSG dofs for easier comparison ---------------
#    Here ----------------->    Original DSG
# [w1,w2,w3,kirchhoffDispFieldx1,...]'  -->    [w1,kirchhoffDispFieldx1,kirchhoffDispFieldy1,w2,...]'
B_original_DSG_ordering = sp.zeros(2,9)
for rmDispFieldma in range(2):
    for node in range(3):
        B_original_DSG_ordering[rmDispFieldma,node*3] = B[rmDispFieldma,node]
        B_original_DSG_ordering[rmDispFieldma,node*3+1] = B[rmDispFieldma,3+node]
        B_original_DSG_ordering[rmDispFieldma,node*3+2] = B[rmDispFieldma,6+node]

print("\n\n\nB Matrix (cartesian space, factor of 1/detJ taken out, ordered as per original DSG formulation) = \n")
sp.pprint(sp.factor(B_original_DSG_ordering)*detJ,wrap_line=False) #detJ taken out for clarity