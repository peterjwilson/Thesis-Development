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

a = sp.Symbol('a') #xi2 - xi1
b = sp.Symbol('b') #eta2 - eta1
c = sp.Symbol('c') #eta3 - eta1
d = sp.Symbol('d') #xi3 - xi1
detJ = sp.Symbol('detJ')


#Nodal coords
xi1 = sp.Symbol('xi1')
eta1 = sp.Symbol('eta1')
xi2 = sp.Symbol('xi2')
eta2 = sp.Symbol('eta2')
xi3 = sp.Symbol('xi3')
eta3 = sp.Symbol('eta3')


#Node coordinates set to example values
xi1 = 0.0
eta1 = 0.0
xi2 = 1.0
eta2 = 0.0
xi3 = 0.0
eta3 = 1.0


#Shape functions
N1 = (1.0-xi-eta)
N2 = xi
N3 = eta



# -----------------------------------------------------------------------------
# Section 1 : Define diplacement fields 
# -----------------------------------------------------------------------------

kirchhoffDispField = a1 + a2*xi +a3*eta + a4*xi**2 + 0.5*(a5+a6)*xi*eta + a7*eta**2 # original Bletzinger


#rmDispField = a8*xi + a9*eta - Bub*xi*eta #generic bubble mode # original Bletzinger
#rmDispField = a8*xi + a9*eta - (a8+a9)*xi*eta #generic bubble mode # original Bletzinger
#rmDispField = a8*xi + a9*eta #no bubble # original Bletzinger no bubbl field
rmDispField = N1*W1 + xi*W2 + eta*W3 + a8*xi + a9*eta # modified by explicitly including displacements, this means the constants a8 and a9 are only rotations

moderatorDispField = 0.5*(a5-a6)*N2*N3 # original Bletzinger
C = a5 - a6


# -----------------------------------------------------------------------------
# Section 2 : Introduce moderator field into the Kirchhoff and RM disp fields 
# -----------------------------------------------------------------------------

PX = kirchhoffDispField + moderatorDispField
PY = kirchhoffDispField - moderatorDispField
GX = rmDispField + moderatorDispField
GY = rmDispField - moderatorDispField



# -----------------------------------------------------------------------------
# Section 3 : Differentiate displacement fields to recover shear deformation fields
# -----------------------------------------------------------------------------

kirchhoffShearField_xi = sp.diff(PX,xi)
kirchhoffShearField_eta = sp.diff(PY,eta)
rmShearField_xi = sp.diff(GX,xi)
rmShearField_eta = sp.diff(GY,eta)


# -----------------------------------------------------------------------------
# Section 4 : Setup equations
# -----------------------------------------------------------------------------

# Alternative representation of the shear gap field
shearGapField = rmDispField - kirchhoffDispField
shearGapField_xi = GX - PX
shearGapField_eta = GY - PY

#Identification of the Ansatz coefficients with the node values for shifts Wi and rotations
# (converted to zero value eqns)

# Shear gaps at nodes should be zero
eq1 = W1 - (shearGapField.subs(xi,xi1)).subs(eta,eta1)
eq2 = W2 - (shearGapField.subs(xi,xi2)).subs(eta,eta2) #xi2, eta2
eq3 = W3 - (shearGapField.subs(xi,xi3)).subs(eta,eta3) #xi3, eta3

# Nodal shear deformations should be zero
eq4 = PX1 - (kirchhoffShearField_xi.subs(xi,xi1)).subs(eta,eta1)
eq5 = PX2 - (kirchhoffShearField_xi.subs(xi,xi2)).subs(eta,eta2)
eq6 = PX3 - (kirchhoffShearField_xi.subs(xi,xi3)).subs(eta,eta3)

eq7 = PY1 - (kirchhoffShearField_eta.subs(xi,xi1)).subs(eta,eta1)
eq8 = PY2 - (kirchhoffShearField_eta.subs(xi,xi2)).subs(eta,eta2)
eq9 = PY3 - (kirchhoffShearField_eta.subs(xi,xi3)).subs(eta,eta3)

print("eq1 = ",eq1)
print("eq2 = ",eq2)
print("eq3 = ",eq3)
print("eq4 = ",eq4)
print("eq5 = ",eq5)
print("eq6 = ",eq6)
print("eq7 = ",eq7)
print("eq8 = ",eq8)
print("eq9 = ",eq9)



# -----------------------------------------------------------------------------
# Section 5 : Solve equations and recover Ansatz coefficients
# -----------------------------------------------------------------------------

UVector = sp.zeros(9,1) #vector of displacements, plate theory rotations still used here
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
    print("a",i+1,"=\t",ansatzCoefficients[i])
    

    
# -----------------------------------------------------------------------------
# Section 6 : Update previously defined fields with recovered Ansatz coefficients
# -----------------------------------------------------------------------------
# updated fields denoted by appended "_u"

kirchhoffDispField_u = kirchhoffDispField.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
#print("\n kirchhoffDispField_u = ",kirchhoffDispField_u)

rmDispField_u = rmDispField.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
#print("\nrmDispField_u = ",rmDispField_u)

moderatorDispField_u = moderatorDispField.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
print("\nmoderatorDispField_u = ",moderatorDispField_u)

C_u = C.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])



# -----------------------------------------------------------------------------
# Section 7 : Construct the updated RM (+/-) moderator displacement fields
# -----------------------------------------------------------------------------

#Currently the RM disp field is overwritten with a field based on that obtained
#form the original DSG formulation. The only difference is instead of representing
#the displacement gap it represents the displacement (ie. +W1 to all nodes)

#The rotation entries are all multiplied with lengths (a,b,c,d), resulting
# in transverse displacement entries. Correct, since we're defining a displacement field
rmDispField_w = N1*W1 + N2*W2 + N3*W3
rmDispField_xi = N2*0.5*(a*(PX1+PX2) + b*(PY1 + PY2))
rmDispField_eta = N3*0.5*(d*(PX1+PX3) + c*(PY1 + PY3))
rmDispField_u = rmDispField_w + rmDispField_xi + rmDispField_eta
GX_rm = rmDispField_u

#The rotation entries in this field are NOT multiplied with lengths (a,b,c,d), resulting
#in incorrect rotation entries to the displacement field. In fact, they only become 
#correct displacement entries in the right angle unit triangle special case (because side lengths are 1.0).
#If the side lengths are 1.0 but the triangle is generally skewed, entries will be missing.
GX_mod = moderatorDispField_u
GX_u = GX_rm + GX_mod
GX_u = GX_u.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])


# rmDispField_u defined above
GY_rm = rmDispField_u

GY_mod = moderatorDispField_u 
GY_u = GY_rm - GY_mod
GY_u = GY_u.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])



# -----------------------------------------------------------------------------
# Section 8 : Switch from plate theory rotational DOFs to geometric rotations
# -----------------------------------------------------------------------------

# vector of displacements is updated
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

# Sub in Kratos rotations
GX_u = GX_u.subs([(PY1,-RX1),(PY2,-RX2),(PY3,-RX3),(PX1,RY1),(PX2,RY2),(PX3,RY3)])
GY_u = GY_u.subs([(PY1,-RX1),(PY2,-RX2),(PY3,-RX3),(PX1,RY1),(PX2,RY2),(PX3,RY3)])


# Jacobian entries
dxi_dx = c/detJ
deta_dx = -b/detJ
dxi_dy = -d/detJ
deta_dy = a/detJ



# -----------------------------------------------------------------------------
# Section 9 : Cartesian transformation
# -----------------------------------------------------------------------------

# Big Cartesian transformation - I think this is incorrect
#rmShearField_xi_u = sp.diff(GX_u,xi)*dxi_dx + sp.diff(GX_u,eta)*deta_dx + sp.diff(GY_u,xi)*dxi_dx + sp.diff(GY_u,eta)*deta_dx
#rmShearField_eta_u = sp.diff(GY_u,xi)*dxi_dy + sp.diff(GY_u,eta)*deta_dy + sp.diff(GX_u,xi)*dxi_dy + sp.diff(GX_u,eta)*deta_dy

# Cartesian transformation used in Original DSG derivation
rmShearField_xi_u = sp.diff(GX_u,xi)*c/detJ + sp.diff(GY_u,eta)*-b/detJ
rmShearField_eta_u = sp.diff(GY_u,eta)*a/detJ + sp.diff(GX_u,xi)*-d/detJ



# -----------------------------------------------------------------------------
# Section 10 : Print out results
# -----------------------------------------------------------------------------

print("\n\n\nSymbolic B Matrix ( transformed to cartesian ) (1/detJ taken out) = ")

# Assemble B Matrix ------------------------------------------
B = sp.zeros(2,9)
for col in range(9):
    B[0,col] = sp.diff(rmShearField_xi_u,UVector[col])
    B[1,col] = sp.diff(rmShearField_eta_u,UVector[col])

sp.pprint(sp.simplify(B*detJ),wrap_line=False)


print("\n\n\n\n\nPrinting individual entries of original Bmat (1/detJ taken out), just for easy copying into C++:")
B = B.subs([(xi,loc1),(eta,loc2)])
Bsimp = sp.factor(B)*detJ
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