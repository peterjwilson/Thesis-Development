# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:26:16 2017

@author: Peter Wilson
"""

import numpy as np
import sympy as sp
print("sympy version: ",sp.__version__)
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

myfontsize = 20
labelfontsize = 14
#===============================================================================
## Set style of plots
font = {'family' : 'sans-serif',
        'family' : 'Arial',
        'style'  : 'normal',
        'weight' : 'normal',
        'size'   : labelfontsize}

plt.rc('font', **font)
plt.rc('font', serif='Helvetica Neue') 
#plt.rcParams['font.family'] = 'Arial'
#plt.ticklabel_format(style='sci', axis='y', scilimits=(1e-4,10000))
#===============================================================================

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
rmDispField = a8*xi + a9*eta #no bubble # original Bletzinger no bubbl field
#rmDispField = N1*W1 + xi*W2 + eta*W3 + a8*xi + a9*eta # modified by explicitly including displacements, this means the constants a8 and a9 are only rotations

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

print("shearGapField = ",sp.simplify(shearGapField))

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

shearGapField_u = shearGapField.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

kirchhoffDispField_u = kirchhoffDispField.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
#print("\n kirchhoffDispField_u = ",kirchhoffDispField_u)

rmDispField_u = rmDispField.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
print("\nrmDispField_u = ",rmDispField_u)

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
#rmDispField_u = rmDispField_w + rmDispField_xi + rmDispField_eta
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

## vector of displacements is updated
#UVector = sp.zeros(9,1) 
#UVector[0] = W1
#UVector[1] = W2
#UVector[2] = W3
#UVector[3] = RX1
#UVector[4] = RX2
#UVector[5] = RX3
#UVector[6] = RY1
#UVector[7] = RY2
#UVector[8] = RY3
#print("Vector of displacements (rotations are per FEM):\n",UVector)
#
## Sub in Kratos rotations
#GX_u = GX_u.subs([(PY1,-RX1),(PY2,-RX2),(PY3,-RX3),(PX1,RY1),(PX2,RY2),(PX3,RY3)])
#GY_u = GY_u.subs([(PY1,-RX1),(PY2,-RX2),(PY3,-RX3),(PX1,RY1),(PX2,RY2),(PX3,RY3)])


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
#rmShearField_xi_u = sp.diff(GX_u,xi)*c/detJ + sp.diff(GY_u,eta)*-b/detJ
#rmShearField_eta_u = sp.diff(GY_u,eta)*a/detJ + sp.diff(GX_u,xi)*-d/detJ

# Parametric only differentiation
rmShearField_xi_u = sp.diff(GX_u,xi)
rmShearField_eta_u = sp.diff(GY_u,eta)

print("\ parametric only rmShearField_xi_u = ",rmShearField_xi_u)
print("\ parametric only rmShearField_eta_u = ",rmShearField_eta_u)



# -----------------------------------------------------------------------------
# Section 10 : Print out results
# -----------------------------------------------------------------------------

print("\n\n\nSymbolic B Matrix ( transformed to cartesian ) (1/detJ taken out) = ")

# Assemble B Matrix ------------------------------------------
B = sp.zeros(2,9)
for col in range(9):
    B[0,col] = sp.diff(rmShearField_xi_u,UVector[col])
    B[1,col] = sp.diff(rmShearField_eta_u,UVector[col])

sp.pprint(sp.simplify(B),wrap_line=False)


# -----------------------------------------------------------------------------
# Section 11 : Solve example problem
# -----------------------------------------------------------------------------


disps = sp.zeros(9,1)
disps[0] = 0.5 #
disps[1] = 1.0 #
disps[2] = 0.5 #
disps[3] = 1.0 #
disps[4] = -2.0 #solved for
disps[5] = 1.5
disps[6] = 1.0 #
disps[7] = 1.5 #solved for
disps[8] = -1.0 #solved for

C = sp.zeros(2,2)
E = 1000
G = E/2.0 #poisson = 0.0
h = 0.05
C[0,0] = 1
C[1,1] = 1
C *= G*5.0/6.0*h

K = B.T * C * B
K1 = sp.integrate(K,(eta,0.0,1.0-xi))
Kint = sp.integrate(K1,(xi,0.0,1.0))
print("\n\nCartesian stiffness matrix integrated symbolically:")
sp.pprint(sp.N(Kint,3),wrap_line=False)

#setup reduced system and get rid of prescribed displacements
redK = sp.zeros(3,3)
redF = sp.zeros(3,1)
displacedDOFs = [1,2,3,4,6,7] #dof number, not index
freeDOFs = [5,8,9] #dof number, not index
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

print("\n\n\n=====================================================================\n")
print("\tNOW RE-RUNNING PROBLEM INCLUDING BENDING FORMULATION")
print("\n=====================================================================\n")

#Bletzinger maple worksheet DOF order maintained here

B_bend_and_shear = sp.zeros(5,9)

a = xi2 - xi1
b = eta2 - eta1
c = eta3 - eta1
d = xi3 - xi1

#Add in bending B matrix from DSG, re-arranged for maple DOF ordering

#phix's
B_bend_and_shear[0,3] = b-c 
B_bend_and_shear[1,3] = 0.0
B_bend_and_shear[2,3] = d-a

B_bend_and_shear[0,4] = c
B_bend_and_shear[1,4] = 0.0
B_bend_and_shear[2,4] = -d

B_bend_and_shear[0,5] = -b
B_bend_and_shear[1,5] = 0.0
B_bend_and_shear[2,5] = a

#phiy's
B_bend_and_shear[0,6] = 0.0
B_bend_and_shear[1,6] = d-a
B_bend_and_shear[2,6] = b-c

B_bend_and_shear[0,7] = 0.0
B_bend_and_shear[1,7] = -d
B_bend_and_shear[2,7] = c

B_bend_and_shear[0,8] = 0.0
B_bend_and_shear[1,8] = a
B_bend_and_shear[2,8] = -b

area = 1.0*1.0/2.0
B_bend_and_shear /= (2.0*area) #divide bending part of Bmat by detJ

#Copy in shear B matrix, already calculated above, no scaling needed
for row in range(2):
    for col in range(9):
        B_bend_and_shear[3+row,col] = B[row,col]

# Setup combined material matrix
C_bend_and_shear = sp.zeros(5,5)
C_bend_and_shear[0,0] = 1.0
C_bend_and_shear[1,1] = 1.0
C_bend_and_shear[0,1] = 0.0
C_bend_and_shear[1,0] = 0.0
C_bend_and_shear[2,2] = (1.0-0.0)/2.0
C_bend_and_shear *= (E*h**3/(12))
for row in range(2):
    for col in range(2):
        C_bend_and_shear[3+row,3+col] = C[row,col]

#Calc combined K
K_temp3 = B_bend_and_shear.T * C_bend_and_shear * B_bend_and_shear
Ktemp1 = sp.integrate(K_temp3,(eta,0.0,1.0-xi))
K_bend_and_shear = sp.integrate(Ktemp1,(xi,0.0,1.0))
print("\n\nCombined stiffness matrix :")
sp.pprint(sp.N(K_bend_and_shear,3),wrap_line=False)

#setup reduced system and get rid of prescribed displacements
for DOFrowcheck in range(len(freeDOFs)):
    tempF = 0.0
    row = freeDOFs[DOFrowcheck] - 1
    for DOFcolcheck in range(len(displacedDOFs)):
        col = displacedDOFs[DOFcolcheck] - 1
        tempF += K_bend_and_shear[row,col] * disps[col]
    redF[DOFrowcheck] = 0.0 - tempF

for row in range(len(freeDOFs)):
    for col in range(len(freeDOFs)):
        redK[row,col] = K_bend_and_shear[freeDOFs[row]-1,freeDOFs[col]-1]

solvedDisps = (redK**-1)*redF
print("\nSolved displacements:")
sp.pprint(solvedDisps)




# -----------------------------------------------------------------------------
# Section 12 : Plot results
# -----------------------------------------------------------------------------

print('\n\n\n\ntotal displacement vector  = ')
for i in range(9):
    print(disps[i])

shearGapField_u = shearGapField_u.subs([(W1,disps[0]),(W2,disps[1]),(W3,disps[2]),(PX1,disps[3]),(PX2,disps[4]),(PX3,disps[5]),(PY1,disps[6]),(PY2,disps[7]),(PY3,disps[8])])
print('\n\nshearGapField_u = ',shearGapField_u)

kirchhoffDispField_u = kirchhoffDispField_u.subs([(W1,disps[0]),(W2,disps[1]),(W3,disps[2]),(PX1,disps[3]),(PX2,disps[4]),(PX3,disps[5]),(PY1,disps[6]),(PY2,disps[7]),(PY3,disps[8])])
print('\n\nkirchhoffDispField_u = ',kirchhoffDispField_u)

rmDispField_u = rmDispField_u.subs([(W1,disps[0]),(W2,disps[1]),(W3,disps[2]),(PX1,disps[3]),(PX2,disps[4]),(PX3,disps[5]),(PY1,disps[6]),(PY2,disps[7]),(PY3,disps[8])])
print('\n\nrmDispField_u = ',rmDispField_u)

moderatorDispField_u=moderatorDispField_u.subs([(W1,disps[0]),(W2,disps[1]),(W3,disps[2]),(PX1,disps[3]),(PX2,disps[4]),(PX3,disps[5]),(PY1,disps[6]),(PY2,disps[7]),(PY3,disps[8])])
print('\n\nmoderatorDispField_u = ',moderatorDispField_u)

rmShearField_xi_u = rmShearField_xi_u.subs([(W1,disps[0]),(W2,disps[1]),(W3,disps[2]),(PX1,disps[3]),(PX2,disps[4]),(PX3,disps[5]),(PY1,disps[6]),(PY2,disps[7]),(PY3,disps[8])])
print('\n\nrmShearField_xi_u = ',rmShearField_xi_u)


kirchhoffShearField_xi_u = kirchhoffShearField_xi.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
kirchhoffShearField_xi_u = kirchhoffShearField_xi_u.subs([(W1,disps[0]),(W2,disps[1]),(W3,disps[2]),(PX1,disps[3]),(PX2,disps[4]),(PX3,disps[5]),(PY1,disps[6]),(PY2,disps[7]),(PY3,disps[8])])
print('\n\nkirchhoffShearField_xi_u = ',kirchhoffShearField_xi_u)


kirchhoffShearField_eta_u = kirchhoffShearField_eta.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])
kirchhoffShearField_eta_u = kirchhoffShearField_eta_u.subs([(W1,disps[0]),(W2,disps[1]),(W3,disps[2]),(PX1,disps[3]),(PX2,disps[4]),(PX3,disps[5]),(PY1,disps[6]),(PY2,disps[7]),(PY3,disps[8])])
print('\n\nkirchhoffShearField_eta_u = ',kirchhoffShearField_eta_u)


divisions = 8
u = np.linspace(0, 1, divisions)
v = np.linspace(1, 0, divisions)

x_mesh = np.outer(np.ones_like(u), v)
y_mesh = np.outer(v, u)

divisions = 10
u = np.linspace(0, 1, divisions)
v = np.linspace(1, 0, divisions)

x = np.outer(np.ones_like(u), v)
y = np.outer(v, u)


#cm.coolwarm

fig = plt.figure(1)

# shearGapField_u
z = 1.0*y**2 - 0.5*y*x - 1.0*y + 1.5*x**2 - 1.0*x + 0.5
ax = fig.gca(projection='3d')
ax.plot_wireframe(x_mesh,y_mesh,0.0)
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet,
    linewidth=0, antialiased=True)

#ax.set_zlabel('shearGapField_u',fontsize=myfontsize)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_zticks([0, 0.25, 0.5, 0.75, 1])
plt.xlabel('xi')
plt.ylabel('eta')
ax.set_zlim(0, 1)
surf.set_clim(vmin=0, vmax=1)
fig.colorbar(surf, shrink=0.75, aspect=20)
plt.title('Displacement gap field')

plt.savefig('displacement_gap_field.pdf',bbox_inches="tight")
ax.elev = 89
ax.azim = -90
ax.w_zaxis.line.set_lw(0.)
ax.set_zticks([])
plt.savefig('displacement_gap_field_90.pdf',bbox_inches="tight")





fig = plt.figure(2)

# kirchhoffDispField_u
z = -1.0*y**2 + 0.5*y*x + 1.0*y - 1.5*x**2 + 1.0*x - 0.5
ax = fig.gca(projection='3d')
ax.plot_wireframe(x_mesh,y_mesh,0.0)
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet,
    linewidth=0, antialiased=True)
#ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_zticks([-1, -0.75, -0.5, -0.25, 0])
plt.xlabel('xi')
plt.ylabel('eta')
ax.set_zlim(-1, 0)
surf.set_clim(vmin=-1, vmax=0)
fig.colorbar(surf, shrink=0.75, aspect=20)
plt.title('KT displacement field')

plt.savefig('kt_displacement_field.pdf',bbox_inches="tight")
ax.elev = 89
ax.azim = -90
ax.w_zaxis.line.set_lw(0.)
ax.set_zticks([])
plt.savefig('kt_displacement_field_90.pdf',bbox_inches="tight")


fig = plt.figure(3)

# rmDispField_u
z = 0.0
ax = fig.gca(projection='3d')
ax.plot_wireframe(x_mesh,y_mesh,0.0)
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet,
    linewidth=0, antialiased=True)
#ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_zticks([-1, 0,  1])
plt.xlabel('xi')
plt.ylabel('eta')
ax.set_zlim(-1, 1)
surf.set_clim(vmin=-1, vmax=1)
fig.colorbar(surf, shrink=0.75, aspect=20)
plt.title('RM shear-displacement field')

plt.savefig('rm_shear_disp_field.pdf',bbox_inches="tight")
#ax.elev = 89
#ax.azim = -90
#ax.w_zaxis.line.set_lw(0.)
#ax.set_zticks([])
#plt.savefig('rm_shear_disp_field_90.pdf',bbox_inches="tight")


fig = plt.figure(4)

# addition of fields
z = 1.0*y**2 - 0.5*y*x - 1.0*y + 1.5*x**2 - 1.0*x + 0.5
ax = fig.gca(projection='3d')
ax.plot_wireframe(x_mesh,y_mesh,0.0)
surf1 = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet,
    linewidth=0, antialiased=True)
z = -1.0*y**2 + 0.5*y*x + 1.0*y - 1.5*x**2 + 1.0*x - 0.5
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet,
    linewidth=0, antialiased=True)
z = 0.0
ax = fig.gca(projection='3d')
surf3 = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet,
    linewidth=0, antialiased=True)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

ax.set_zticks([-1, -0.5, 0, 0.5, 1])
plt.xlabel('xi')
plt.ylabel('eta')
ax.set_zlim(-1, 1)
#fig.colorbar(surf1, shrink=1, aspect=10)
surf1.set_clim(vmin=-1, vmax=1)
surf2.set_clim(vmin=-1, vmax=1)
surf3.set_clim(vmin=-1, vmax=1)
fig.colorbar(surf, shrink=0.75, aspect=20)
plt.title('Displacement gap + KT = RM shear-displacement')

plt.savefig('sum.pdf',bbox_inches="tight")
#ax.elev = 89
#ax.azim = -90
#ax.w_zaxis.line.set_lw(0.)
#ax.set_zticks([])
#plt.savefig('sum_90.pdf',bbox_inches="tight")




fig = plt.figure(5)

# kirchhoffShearField_xi_u
z = 0.5*y - 3.0*x + 1.0
ax = fig.gca(projection='3d')
ax.plot_wireframe(x_mesh,y_mesh,0.0)
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet,
    linewidth=0, antialiased=True)
#ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_zticks([-2, -1.5, -1, -0.5, 0, 0.5, 1 ,1.5])
plt.xlabel('xi')
plt.ylabel('eta')
ax.set_zlim(-2, 1.5)
surf.set_clim(vmin=-2, vmax=1.5)
fig.colorbar(surf, shrink=0.75, aspect=20)
ax.set_zlim(-2, 1.5)
surf.set_clim(vmin=-2, vmax=1.5)
plt.title('Rotation xi field')

plt.savefig('rotation_xi.pdf',bbox_inches="tight")
ax.elev = 89
ax.azim = -90
ax.w_zaxis.line.set_lw(0.)
ax.set_zticks([])
plt.savefig('rotation_xi_90.pdf',bbox_inches="tight")



fig = plt.figure(6)

# kirchhoffShearField_eta_u
z = -2.0*y + 0.5*x + 1.0
ax = fig.gca(projection='3d')
ax.plot_wireframe(x_mesh,y_mesh,0.0)
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet,
    linewidth=0, antialiased=True)
#ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_zlim(-1, 1.5)
surf.set_clim(vmin=-1, vmax=1.5)
fig.colorbar(surf, shrink=0.75, aspect=20)
ax.set_zticks([-1, -0.5, 0, 0.5, 1 ,1.5])
plt.xlabel('xi')
plt.ylabel('eta')
ax.set_zlim(-1, 1.5)
surf.set_clim(vmin=-1, vmax=1.5)
plt.title('Rotation eta field')

plt.savefig('rotation_eta.pdf',bbox_inches="tight")
ax.elev = 89
ax.azim = -90
ax.w_zaxis.line.set_lw(0.)
ax.set_zticks([])
plt.savefig('rotation_eta_90.pdf',bbox_inches="tight")






plt.show()


