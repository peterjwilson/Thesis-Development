# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:26:16 2017

@author: Peter Wilson
"""

import numpy as np
import sympy as sp

x = sp.Symbol('x')
y = sp.Symbol('y')

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



def assembleEquationIntoMatrix( equation, matrix, UVector, row):
    row -= 1
    for col in range(9):
        filterVector = sp.zeros(9,1)
        filterVector[col] = 1
        test = equation.subs([(a1,filterVector[0]),(a2,filterVector[1]),(a3,filterVector[2]),(a4,filterVector[3]),(a5,filterVector[4]),(a6,filterVector[5]),(a7,filterVector[6]),(a8,filterVector[7]),(a9,filterVector[8])])
        matrix[row,col] =sp.solve(test,UVector[row])

    

#Node coordinates
x1 = 0.0
y1 = 0.0
x2 = 1.0
y2 = 0.0
x3 = 0.0
y3 = 1.0

# Section 1 ---------------------------------------------------
PHI = a1 + a2*x +a3*y + a4*x**2 + 0.5*(a5+a6)*x*y + a7*y**2
#GAM = a8*x + a9*y - Bub*x*y #generic bubble mode
GAM = a8*x + a9*y - (a8+a9)*x*y #generic bubble mode
PD = 0.5*(a5-a6)*x*y
C = a5 - a6

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
UVector = sp.zeros(9,1) #vector of equations
UVector[0] = W1
UVector[1] = W2
UVector[2] = W3
UVector[3] = PX1
UVector[4] = PX2
UVector[5] = PX3
UVector[6] = PY1
UVector[7] = PY2
UVector[8] = PY3

#assemble coefficients into matrix
assembleEquationIntoMatrix(eq1,A,UVector,1)
assembleEquationIntoMatrix(eq2,A,UVector,2)
assembleEquationIntoMatrix(eq3,A,UVector,3)
assembleEquationIntoMatrix(eq4,A,UVector,4)
assembleEquationIntoMatrix(eq5,A,UVector,5)
assembleEquationIntoMatrix(eq6,A,UVector,6)
assembleEquationIntoMatrix(eq7,A,UVector,7)
assembleEquationIntoMatrix(eq8,A,UVector,8)
assembleEquationIntoMatrix(eq9,A,UVector,9)
#sp.pprint(A)

#solve
ansatzCoefficients = (A**-1)*UVector
#.pprint(ansatzCoefficients)

#assign ansatz coefficients into a1 -> a9 and sub into px, py, gx, gy
pxSub = px.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

pySub = py.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

gxSub = gx.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

gySub = gy.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])


#Go through and update everything "_u"
PHI_u = PHI.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

GAM_u = GAM.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

PD_u = PD.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

C_u = C.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

PX_u = PHI_u + PD_u
PY_u = PHI_u - PD_u
GX_u = GAM_u + PD_u
GY_u = GAM_u - PD_u

# Section 3 ---------------------------------------------------
px_u = sp.diff(PX_u,x)
py_u = sp.diff(PY_u,y)
gx_u = sp.diff(GX_u,x)
gy_u = sp.diff(GY_u,y)


# Assemble B Matrix ------------------------------------------
B = sp.zeros(2,9)
for col in range(9):
    B[0,col] = sp.diff(gx_u,UVector[col])
    B[1,col] = sp.diff(gy_u,UVector[col])

sp.pprint(B,wrap_line=False)
