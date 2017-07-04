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


a = sp.Symbol('a') #x2 - x1
b = sp.Symbol('b') #y2 - y1
c = sp.Symbol('c') #y3 - y1
d = sp.Symbol('d') #x3 - x1
detJ = sp.Symbol('detJ')
A = sp.Symbol('A')

#Shape functions
N1 = 1.0/detJ*((x2*y3-x3*y2)+x*(y2-y3)+y*(x3-x2))
N2 = 1.0/detJ*((x3*y1-x1*y3)+x*(y3-y1)+y*(x1-x3))
N3 = 1.0/detJ*((x1*y2-x2*y1)+x*(y1-y2)+y*(x2-x1))

# Displacement field
W_test = N1*W1 + N2*W2 + N3*W3

# Rotation field
phi_x = N1*PX1 + N2*PX2 + N3*PX3 
phi_y = N1*PY1 + N2*PY2 + N3*PY3 

# Shear gaps at each node
sg_x_n1 = W_test.subs([(x,0)])- W_test.subs([(x,0)]) + sp.integrate(phi_x*a + phi_y*b,(x,x1,x1)).subs([(x,0),(y,y1)])

sg_x_n2 = W_test.subs([(x,x2)])- W_test.subs([(x,0)]) + sp.integrate(phi_x*a + phi_y*b,(x,0,x2)).subs([(x,x2),(y,y2)])

sg_x_n3 = W_test.subs([(x,x3)])- W_test.subs([(x,0)]) + sp.integrate(phi_x*a + phi_y*b,(x,0,x3)).subs([(x,x3),(y,y3)])


sg_y_n1 = W_test.subs([(x,0),(y,y1)]) - W_test.subs([(x,0),(y,y1)]) + sp.integrate(phi_x*d + phi_y*c,(y,y1,y1)).subs([(x,0),(y,y1)])

sg_y_n2 = W_test.subs([(x,x2),(y,y2)])- W_test.subs([(x,x2),(y,y1)]) + sp.integrate(phi_x*d + phi_y*c,(y,y1,y2)).subs([(x,x2),(y,y2)])

sg_y_n3 = W_test.subs([(x,x3),(y,y3)])- W_test.subs([(x,x3),(y,y1)]) + sp.integrate(phi_x*d + phi_y*c,(y,y1,y3)).subs([(x,x3),(y,y3)])

print("\n\n Shear gaps evaluated at each node =:")

print("\n\n sg_x_n1 =:",sp.simplify(sg_x_n1))

print("\n\n sg_x_n2 =:",sp.simplify(sg_x_n2))

print("\n\n sg_x_n3 =:",sp.simplify(sg_x_n3))

print("\n\n sg_y_n1 =:",sg_y_n1)

print("\n\n sg_y_n2 =:",sg_y_n2)

print("\n\n sg_y_n3 =:",sg_y_n3)



# Assemble shear gap field ----------------------------------------------------
sg_x_sf = N1*sg_x_n1 + N2*sg_x_n2 + N3*sg_x_n3
sg_y_sf = N1*sg_y_n1 + N2*sg_y_n2 + N3*sg_y_n3
print("\n\n\n\n\n\n")



# Shear strain field ----------------------------------------------------------
gx_u = sp.diff(sg_x_sf,x)*c/detJ + sp.diff(sg_y_sf,y)*-b/detJ
gy_u = sp.diff(sg_y_sf,y)*a/detJ + sp.diff(sg_x_sf,x)*-d/detJ

gx_u = sp.diff(sg_x_sf,x)
gy_u = sp.diff(sg_y_sf,y)

#gx_u = sg_x_n2*c/detJ + sg_y_n3*-b/detJ
#gy_u = sg_y_n3*a/detJ + sg_x_n2*-d/detJ

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


print("\n\n\n\n\nPrinting individual entries of matrix above (1/detJ taken out), just for easy copying into C++:")
#Bsimp = sp.factor(B)*detJ #detJ taken out for clarity
#testing below
Bsimp = sp.simplify(B_original_DSG_ordering)*detJ
#testing end
for col in range(9):
    print("BSuper(0,",col,")=",Bsimp[0,col],";")
for col in range(9):
    print("BSuper(1,",col,")=",Bsimp[1,col],";")


#print("\n\n\n\n\nPrinting individual entries of matrix above, just for easy copying into C++:\n")
#Bsimp = sp.factor(B)*detJ #detJ taken out for clarity
#B_C_plus_plus = sp.zeros(2,18)
#for gamma in range(2):
#    for node in range(3):
#        B_C_plus_plus[gamma,2+node*6] = B[gamma,node]
#        B_C_plus_plus[gamma,3+node*6] = B[gamma,node+3]
#        B_C_plus_plus[gamma,4+node*6] = B[gamma,node+6]
#
#B_C_plus_plus *= detJ
#B_C_plus_plus = sp.simplify(B_C_plus_plus)
#for col in range(18):
#    print("data.B(6,",col,")=",B_C_plus_plus[0,col],";")
#for col in range(18):
#    print("data.B(7,",col,")=",B_C_plus_plus[1,col],";")
#    
#    