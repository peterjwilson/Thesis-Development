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

xi = sp.Symbol('xi')
eta = sp.Symbol('eta')

xi1 = sp.Symbol('xi1')
xi2 = sp.Symbol('xi2')
xi3 = sp.Symbol('xi3')
eta1 = sp.Symbol('eta1')
eta2 = sp.Symbol('eta2')
eta3 = sp.Symbol('eta3')


#N4 = 27*xi*eta*(1.0-xi-eta)
#N1 = (1.0-xi-eta-N4/3.0)
#N2 = (x-N4/3.0)
#N3 = (y-N4/3.0)
#x_sf = N1*x1 + N2*x2 + N3*x3 + (N4)*(x1 + x2 + x3)/3.0
#y_sf = N1*y1 + N2*y2 + N3*y3 + (N4)*(y1 + y2 + y3)/3.0

N1 = (1.0-xi-eta)
N2 = xi
N3 = eta


#Node coordinates set to example values
#x1 = 0.0
#y1 = 0.0
#x2 = 1.0
#y2 = 0.0
#x3 = 0.0
#y3 = 1.0

#Displacement field
#W_test = N1*W1 + N2*W2 + N3*W3 + N4*(W1 + W2 + W3)/3.0
W_test = N1*W1 + N2*W2 + N3*W3
#dW_dX = sp.diff(W_test,x)
#dW_dY = sp.diff(W_test,y)



#
##Gamma field
#gamma_x = N1*PX1 + N2*PX2 + N3*PX3 + N4*(PX1 + PX2 + PX3)/3.0
#gamma_y = N1*PY1 + N2*PY2 + N3*PY3 + N4*(PY1 + PY2 + PY3)/3.0
gamma_xi = N1*PX1 + N2*PX2 + N3*PX3 
gamma_eta = N1*PY1 + N2*PY2 + N3*PY3 

print("\n\n Shear gaps evaluated at each node =:")

sg_xi_n1 = W_test.subs([(xi,xi1),(eta,eta1)])- W_test.subs([(xi,xi1),(eta,eta1)]) + (sp.integrate(gamma_xi*a + gamma_eta*b,(xi,xi1,xi1))).subs([(xi,xi1),(eta1,eta1)])

print("\n\n sg_xi_n1 =:",sg_xi_n1.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)]))

sg_xi_n2 = W_test.subs([(xi,xi2),(eta,eta2)])- W_test.subs([(xi,xi1),(eta,eta2)]) + sp.integrate(gamma_xi*a + gamma_eta*b,(xi,xi1,xi2)).subs([(xi,xi2),(eta,eta2)])

print("\n\n sg_xi_n2 =:",sg_xi_n2.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)]))

sg_xi_n3 = W_test.subs([(xi,xi3),(eta,eta3)])- W_test.subs([(xi,xi1),(eta,eta3)]) + sp.integrate(gamma_xi*a + gamma_eta*b,(xi,xi1,xi3)).subs([(xi,xi3),(eta,eta3)])

print("\n\n sg_xi_n3 =:",sg_xi_n3.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)]))

#sg_x_n4 = W_test.subs([(x,1.0/3.0),(y,1.0/3.0)])- W_test.subs([(x,0.0),(y,1.0/3.0)]) - sp.integrate(gamma_x*a + gamma_y*b,(x,0.0,1.0/3.0)).subs([(x,1.0/3.0),(y,1.0/3.0)])
#print("\n\n sg_x_n4 =:",sg_x_n4)





sg_eta_n1 = W_test.subs([(xi,xi1),(eta,eta1)]) - W_test.subs([(xi,xi1),(eta,eta1)]) + (sp.integrate(gamma_xi*d + gamma_eta*c,(eta,eta1,eta1))).subs([(xi,xi1),(eta,eta1)])
print("\n\n sg_eta_n1 =:",sg_eta_n1.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)]))

sg_eta_n2 = W_test.subs([(xi,xi2),(eta,eta2)])- W_test.subs([(xi,xi2),(eta,eta1)]) + sp.integrate(gamma_xi*d + gamma_eta*c,(eta,eta1,eta2)).subs([(xi,xi2),(eta,eta2)])
print("\n\n sg_eta_n2 =:",sg_eta_n2.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)]))

sg_eta_n3 = W_test.subs([(xi,xi3),(eta,eta3)])- W_test.subs([(xi,xi3),(eta,eta1)]) + sp.integrate(gamma_xi*d + gamma_eta*c,(eta,eta1,eta3)).subs([(xi,xi3),(eta,eta3)])

print("\n\n sg_eta_n3 =:",sg_eta_n3.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)]))

##sg_y_n4 = W_test.subs([(x,1.0/3.0),(y,1.0/3.0)])- W_test.subs([(x,1.0/3.0),(y,0.0)]) - sp.integrate(gamma_x*d + gamma_y*c,(y,0.0,1.0/3.0)).subs([(x,1.0/3.0),(y,1.0/3.0)])
##print("\n\n sg_y_n4 =:",sg_y_n4)


print("\n\nShear gap field:")
#sg_x_sf = N1*sg_x_n1 + N2*sg_x_n2 + N3*sg_x_n3 + N4*sg_x_n4
sg_xi_sf = N1*sg_xi_n1 + N2*sg_xi_n2 + N3*sg_xi_n3
#sg_xi_sf = sg_xi_sf.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)])
print("\n\n sg_xi_sf=",sg_xi_sf.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)]))
#
##sg_y_sf = N1*sg_y_n1 + N2*sg_y_n2 + N3*sg_y_n3 + N4*sg_y_n4
sg_eta_sf = N1*sg_eta_n1 + N2*sg_eta_n2 + N3*sg_eta_n3
#sg_eta_sf = sg_eta_sf.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)])
print("\n\n sg_y_sf=",sg_eta_sf.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)]))

print("\n\n\n\n\n\n")

# Section 3 ---------------------------------------------------
gx_u = sp.diff(sg_xi_sf,xi)
gy_u = sp.diff(sg_eta_sf,eta)

#gx_u = gx_u.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)])


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

sp.pprint(B.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)]),wrap_line=False)



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
    B[0,col] = sp.diff(sp.diff(sg_xi_sf,xi),UVector[col])*c/detJ + sp.diff(sp.diff(sg_xi_sf,eta),UVector[col])*-b/detJ
    B[1,col] = sp.diff(sp.diff(sg_eta_sf,xi),UVector[col])*-d/detJ + sp.diff(sp.diff(sg_eta_sf,eta),UVector[col])*a/detJ

B = B.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)])
sp.pprint(sp.simplify(B)*detJ,wrap_line=False) #detJ taken out for clarity

#print("\n\n\nSymbolic B Matrix ( parametric space, from DSG formulation) = ")
#gx_u = sp.diff(N2,x)*sg_x_n2
#gy_u = sp.diff(N3,y)*sg_y_n3
#B = sp.zeros(2,9)
#for col in range(9):
#    B[0,col] = sp.diff(gx_u,UVector[col])
#    B[1,col] = sp.diff(gy_u,UVector[col])
#
#sp.pprint(B,wrap_line=False)
#
#
#print("\n\n\nSymbolic B Matrix ( cartesian space, from DSG formulation) = ")
#gx_u = sp.diff(N2,xi)*sg_xi_n2*c/detJ + sp.diff(N3,eta)*sg_xi_n3*-1.0*b/detJ
##gx_u = gx_u.subs([(xi1,0),(xi2,1.0),(xi3,0),(eta1,0),(eta2,0),(eta3,1.0)])
#gy_u = sp.diff(N2,xi)*sg_eta_n2*-1*d/detJ + sp.diff(N3,eta)*sg_eta_n3*a/detJ
#
#print("\n\n gx_u = ",gx_u)
#B = sp.zeros(2,9)
#for col in range(9):
#    B[0,col] = sp.diff(gx_u,UVector[col])
#    B[1,col] = sp.diff(gy_u,UVector[col])
#
#sp.pprint(sp.simplify(B*detJ),wrap_line=False)


#Rearraging to original DSG dofs for easier comparison
#    Here ---------------->    Original DSG
# [w1,w2,w3,phix1,...]'  -->    [w1,phix1,phiy1,w2,...]'
B_original_DSG_ordering = sp.zeros(2,9)
for gamma in range(2):
    for node in range(3):
        B_original_DSG_ordering[gamma,node*3] = B[gamma,node]
        B_original_DSG_ordering[gamma,node*3+1] = B[gamma,3+node]
        B_original_DSG_ordering[gamma,node*3+2] = B[gamma,6+node]
B_original_DSG_ordering = B_original_DSG_ordering.subs([(xi1,x1),(xi2,x2),(xi3,x3),(eta1,y1),(eta2,y2),(eta3,y3)])
print("\n\n\nSymbolic B Matrix (cartesian space, factor of 1/detJ taken out, ordered as per original DSG formulation) = ")
sp.pprint(sp.factor(B_original_DSG_ordering)*detJ,wrap_line=False) #detJ taken out for clarity



#B_original_DSG_ordering = B_original_DSG_ordering.subs([(a,x2-x1),(b,y2-y1),(c,y3-x1),(d,x3-x1)])
B_original_DSG_ordering = B_original_DSG_ordering.subs([(x2-x1,a),(y2-y1,b),(y3-x1,c),(x3-x1,d)])

print("\n\n\n\n\nPrinting individual entries of matrix above, just for easy copying into C++:")
Bsimp = sp.factor(B_original_DSG_ordering) #detJ taken out for clarity
for col in range(9):
    print("BSuper(0,",col,")=",Bsimp[0,col],";")
for col in range(9):
    print("BSuper(1,",col,")=",Bsimp[1,col],";")
    
    
    
    
    
