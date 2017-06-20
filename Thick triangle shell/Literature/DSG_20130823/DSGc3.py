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


#Nodal coords
x1 = sp.Symbol('x1')
y1 = sp.Symbol('y1')
x2 = sp.Symbol('x2')
y2 = sp.Symbol('y2')
x3 = sp.Symbol('x3')
y3 = sp.Symbol('y3')


def assembleEquationIntoMatrix( equation, matrix, UVector, row):
    print("Assembling row ",row)
    row -= 1
    for col in range(9):
        filterVector = sp.zeros(9,1)
        filterVector[col] = 1
        test = equation.subs([(a1,filterVector[0]),(a2,filterVector[1]),(a3,filterVector[2]),(a4,filterVector[3]),(a5,filterVector[4]),(a6,filterVector[5]),(a7,filterVector[6]),(a8,filterVector[7]),(a9,filterVector[8])])
        matrix[row,col] =sp.solve(test,UVector[row])

    

#Node coordinates set to example values
x1 = 0.0
y1 = 0.0
x2 = 1.0
y2 = 0.0
x3 = 0.0
y3 = 1.0

# Section 1 ---------------------------------------------------
PHI = a1 + a2*x +a3*y + a4*x**2 + 0.5*(a5+a6)*x*y + a7*y**2

#GAM = a8*x + a9*y - Bub*x*y #generic bubble mode
#GAM = a8*x + a9*y - (a8+a9)*x*y #generic bubble mode
GAM = a8*x + a9*y #no bubble

PD = 0.5*(a5-a6)*x*y
C = a5 - a6

# Section 2 ---------------------------------------------------
PX = PHI + PD
PY = PHI - PD
GX = GAM + PD
GY = GAM - PD

#GX = a8*x + a9*y #normal DSG
#GY = a8*x + a9*y #normal DSG

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

print("Prelims complete")

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

print("Equations setup")
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

#assemble coefficients into matrix
#assembleEquationIntoMatrix(eq1,A,UVector,1)
#assembleEquationIntoMatrix(eq2,A,UVector,2)
#assembleEquationIntoMatrix(eq3,A,UVector,3)
#assembleEquationIntoMatrix(eq4,A,UVector,4)
#assembleEquationIntoMatrix(eq5,A,UVector,5)
#assembleEquationIntoMatrix(eq6,A,UVector,6)
#assembleEquationIntoMatrix(eq7,A,UVector,7)
#assembleEquationIntoMatrix(eq8,A,UVector,8)
#assembleEquationIntoMatrix(eq9,A,UVector,9)
#print("Printing inverse Ansatz coefficients")
#sp.pprint(A,wrap_line=False)

#solve

results = list(sp.linsolve([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[a1,a2,a3,a4,a5,a6,a7,a8,a9]))
ansatzCoefficients = results[0]
print("Ansatz coefficients solved",ansatzCoefficients)

#ansatzCoefficients = sp.zeros(9,1)
#ansatzCoefficients[0] = -(2*W1*x2**2*y3**2 + 2*W1*x3**2*y2**2 + 2*W2*x1**2*y3**2 + 2*W2*x3**2*y1**2 + 2*W3*x1**2*y2**2 + 2*W3*x2**2*y1**2 - 2*W1*x1*x2*y3**2 - 2*W1*x1*x3*y2**2 - 2*W2*x1*x2*y3**2 - 2*W2*x2*x3*y1**2 - 2*W3*x1*x3*y2**2 - 2*W3*x2*x3*y1**2 - 2*W1*x2**2*y1*y3 - 2*W1*x3**2*y1*y2 - 2*W2*x1**2*y2*y3 - 2*W2*x3**2*y1*y2 - 2*W3*x1**2*y2*y3 - 2*W3*x2**2*y1*y3 + PX1*x1*x2**2*y3**2 + PX1*x1*x3**2*y2**2 - PX1*x1**2*x2*y3**2 - PX1*x1**2*x3*y2**2 - PX2*x1*x2**2*y3**2 + PX2*x2*x3**2*y1**2 + PX2*x1**2*x2*y3**2 - PX2*x2**2*x3*y1**2 - PX3*x1*x3**2*y2**2 - PX3*x2*x3**2*y1**2 + PX3*x1**2*x3*y2**2 + PX3*x2**2*x3*y1**2 + PY1*x2**2*y1*y3**2 - PY1*x2**2*y1**2*y3 + PY1*x3**2*y1*y2**2 - PY1*x3**2*y1**2*y2 + PY2*x1**2*y2*y3**2 - PY2*x1**2*y2**2*y3 - PY2*x3**2*y1*y2**2 + PY2*x3**2*y1**2*y2 - PY3*x1**2*y2*y3**2 + PY3*x1**2*y2**2*y3 - PY3*x2**2*y1*y3**2 + PY3*x2**2*y1**2*y3 - 2*W2*x1*x3**2*y1**2 - 2*W3*x1*x2**2*y1**2 - 2*W1*x2*x3**2*y2**2 - 2*W3*x1**2*x2*y2**2 - 2*W1*x2**2*x3*y3**2 - 2*W2*x1**2*x3*y3**2 - 2*W2*x1**2*y1*y3**2 - 2*W3*x1**2*y1*y2**2 - 2*W1*x2**2*y2*y3**2 - 2*W3*x2**2*y1**2*y2 - 2*W1*x3**2*y2**2*y3 - 2*W2*x3**2*y1**2*y3 - PX1*x1*x2*x3**2*y2**2 + PX1*x1**2*x2*x3*y2**2 - PX2*x1*x2*x3**2*y1**2 + PX2*x1*x2**2*x3*y1**2 - PX1*x1*x2**2*x3*y3**2 + PX1*x1**2*x2*x3*y3**2 + PX3*x1*x2*x3**2*y1**2 - PX3*x1*x2**2*x3*y1**2 + PX2*x1*x2**2*x3*y3**2 - PX2*x1**2*x2*x3*y3**2 + PX3*x1*x2*x3**2*y2**2 - PX3*x1**2*x2*x3*y2**2 + PX1*x1*x2**2*y1*y3**2 + PX1*x1*x3**2*y1*y2**2 - 2*PX1*x1*x2**2*y2*y3**2 - PX1*x2*x3**2*y1*y2**2 + PX1*x1**2*x2*y2*y3**2 - PX1*x1**2*x2*y2**2*y3 + PX2*x1*x2**2*y1*y3**2 - PX2*x1*x2**2*y1**2*y3 - PX2*x1*x3**2*y1**2*y2 - 2*PX2*x1**2*x2*y1*y3**2 - 2*PX1*x1*x3**2*y2**2*y3 - PX1*x1**2*x3*y2*y3**2 + PX1*x1**2*x3*y2**2*y3 - PX1*x2**2*x3*y1*y3**2 + PX2*x2*x3**2*y1**2*y2 + PX2*x1**2*x2*y2*y3**2 - PX3*x1*x2**2*y1**2*y3 + PX3*x1*x3**2*y1*y2**2 - PX3*x1*x3**2*y1**2*y2 - 2*PX3*x1**2*x3*y1*y2**2 - 2*PX2*x2*x3**2*y1**2*y3 - PX2*x1**2*x3*y2*y3**2 - PX2*x2**2*x3*y1*y3**2 + PX2*x2**2*x3*y1**2*y3 - PX3*x2*x3**2*y1*y2**2 + PX3*x2*x3**2*y1**2*y2 - PX3*x1**2*x2*y2**2*y3 - 2*PX3*x2**2*x3*y1**2*y2 + PX3*x1**2*x3*y2**2*y3 + PX3*x2**2*x3*y1**2*y3 + PY1*x1*x2**2*y1*y3**2 + PY1*x1*x3**2*y1*y2**2 - PY1*x1*x2**2*y2*y3**2 - 2*PY1*x2*x3**2*y1*y2**2 + PY1*x2*x3**2*y1**2*y2 - PY1*x2**2*x3*y1**2*y2 + PY2*x1*x3**2*y1*y2**2 - 2*PY2*x1*x3**2*y1**2*y2 - PY2*x1**2*x2*y1*y3**2 - PY2*x1**2*x3*y1*y2**2 - PY1*x1*x3**2*y2**2*y3 - PY1*x2*x3**2*y1**2*y3 - 2*PY1*x2**2*x3*y1*y3**2 + PY1*x2**2*x3*y1**2*y3 + PY2*x2*x3**2*y1**2*y2 + PY2*x1**2*x2*y2*y3**2 + PY3*x1*x2**2*y1*y3**2 - 2*PY3*x1*x2**2*y1**2*y3 - PY3*x1**2*x2*y1*y3**2 - PY3*x1**2*x3*y1*y2**2 - PY2*x1*x3**2*y2**2*y3 - PY2*x2*x3**2*y1**2*y3 - 2*PY2*x1**2*x3*y2*y3**2 + PY2*x1**2*x3*y2**2*y3 - PY3*x1*x2**2*y2*y3**2 + PY3*x1**2*x2*y2*y3**2 - 2*PY3*x1**2*x2*y2**2*y3 - PY3*x2**2*x3*y1**2*y2 + PY3*x1**2*x3*y2**2*y3 + PY3*x2**2*x3*y1**2*y3 - PY1*x2**2*y1*y2*y3**2 + PY1*x2**2*y1**2*y2*y3 - PY2*x1**2*y1*y2*y3**2 + PY2*x1**2*y1*y2**2*y3 - PY1*x3**2*y1*y2**2*y3 + PY1*x3**2*y1**2*y2*y3 + PY3*x1**2*y1*y2*y3**2 - PY3*x1**2*y1*y2**2*y3 + PY2*x3**2*y1*y2**2*y3 - PY2*x3**2*y1**2*y2*y3 + PY3*x2**2*y1*y2*y3**2 - PY3*x2**2*y1**2*y2*y3 + 2*W1*x1*x2*y2*y3 + 2*W1*x2*x3*y1*y2 + 2*W2*x1*x2*y1*y3 + 2*W2*x1*x3*y1*y2 - 4*W3*x1*x2*y1*y2 + 2*W1*x1*x3*y2*y3 + 2*W1*x2*x3*y1*y3 - 4*W2*x1*x3*y1*y3 + 2*W3*x1*x2*y1*y3 + 2*W3*x1*x3*y1*y2 - 4*W1*x2*x3*y2*y3 + 2*W2*x1*x3*y2*y3 + 2*W2*x2*x3*y1*y3 + 2*W3*x1*x2*y2*y3 + 2*W3*x2*x3*y1*y2 - PX1*x1*x2**2*y1*y3 - PX1*x1*x3**2*y1*y2 + PX1*x1**2*x2*y2*y3 + PX2*x1*x2**2*y1*y3 + PX1*x1**2*x3*y2*y3 - PX2*x2*x3**2*y1*y2 - PX2*x1**2*x2*y2*y3 + PX3*x1*x3**2*y1*y2 + PX2*x2**2*x3*y1*y3 + PX3*x2*x3**2*y1*y2 - PX3*x1**2*x3*y2*y3 - PX3*x2**2*x3*y1*y3 - PY1*x1*x2*y1*y3**2 - PY1*x1*x3*y1*y2**2 + PY1*x2*x3*y1**2*y2 + PY2*x1*x3*y1*y2**2 + PY1*x2*x3*y1**2*y3 - PY2*x1*x2*y2*y3**2 - PY2*x2*x3*y1**2*y2 + PY3*x1*x2*y1*y3**2 + PY2*x1*x3*y2**2*y3 + PY3*x1*x2*y2*y3**2 - PY3*x1*x3*y2**2*y3 - PY3*x2*x3*y1**2*y3 + 2*W1*x1*x2*x3*y2**2 + 2*W2*x1*x2*x3*y1**2 + 2*W1*x1*x2*x3*y3**2 + 2*W3*x1*x2*x3*y1**2 + 2*W2*x1*x2*x3*y3**2 + 2*W3*x1*x2*x3*y2**2 + 2*W1*x1*x2*y2*y3**2 - 2*W1*x1*x2*y2**2*y3 + 2*W1*x2*x3**2*y1*y2 - 2*W1*x2**2*x3*y1*y2 + 2*W2*x1*x2*y1*y3**2 - 2*W2*x1*x2*y1**2*y3 + 2*W2*x1*x3**2*y1*y2 - 2*W2*x1**2*x3*y1*y2 + 2*W3*x1*x2*y1*y2**2 + 2*W3*x1*x2*y1**2*y2 + 2*W3*x1*x2**2*y1*y2 + 2*W3*x1**2*x2*y1*y2 - 2*W1*x1*x3*y2*y3**2 + 2*W1*x1*x3*y2**2*y3 - 2*W1*x2*x3**2*y1*y3 + 2*W1*x2**2*x3*y1*y3 + 2*W2*x1*x3*y1*y3**2 + 2*W2*x1*x3*y1**2*y3 + 2*W2*x1*x3**2*y1*y3 + 2*W2*x1**2*x3*y1*y3 + 2*W3*x1*x3*y1*y2**2 - 2*W3*x1*x3*y1**2*y2 + 2*W3*x1*x2**2*y1*y3 - 2*W3*x1**2*x2*y1*y3 + 2*W1*x2*x3*y2*y3**2 + 2*W1*x2*x3*y2**2*y3 + 2*W1*x2*x3**2*y2*y3 + 2*W1*x2**2*x3*y2*y3 - 2*W2*x1*x3**2*y2*y3 - 2*W2*x2*x3*y1*y3**2 + 2*W2*x2*x3*y1**2*y3 + 2*W2*x1**2*x3*y2*y3 - 2*W3*x1*x2**2*y2*y3 - 2*W3*x2*x3*y1*y2**2 + 2*W3*x2*x3*y1**2*y2 + 2*W3*x1**2*x2*y2*y3 + 2*W1*x2**2*y1*y2*y3 + 2*W2*x1**2*y1*y2*y3 + 2*W1*x3**2*y1*y2*y3 + 2*W3*x1**2*y1*y2*y3 + 2*W2*x3**2*y1*y2*y3 + 2*W3*x2**2*y1*y2*y3 + PX1*x1*x2*x3*y1*y2 + PX1*x1*x2*x3*y1*y3 + PX2*x1*x2*x3*y1*y2 - 2*PX1*x1*x2*x3*y2*y3 - 2*PX2*x1*x2*x3*y1*y3 - 2*PX3*x1*x2*x3*y1*y2 + PX2*x1*x2*x3*y2*y3 + PX3*x1*x2*x3*y1*y3 + PX3*x1*x2*x3*y2*y3 + PY1*x1*x2*y1*y2*y3 + PY1*x1*x3*y1*y2*y3 + PY2*x1*x2*y1*y2*y3 - 2*PY1*x2*x3*y1*y2*y3 - 2*PY2*x1*x3*y1*y2*y3 - 2*PY3*x1*x2*y1*y2*y3 + PY2*x2*x3*y1*y2*y3 + PY3*x1*x3*y1*y2*y3 + PY3*x2*x3*y1*y2*y3 - 4*W1*x1*x2*x3*y2*y3 - 4*W2*x1*x2*x3*y1*y3 - 4*W3*x1*x2*x3*y1*y2 - 4*W1*x2*x3*y1*y2*y3 - 4*W2*x1*x3*y1*y2*y3 - 4*W3*x1*x2*y1*y2*y3 + PX1*x1*x2*x3**2*y1*y2 - PX1*x1*x2**2*x3*y1*y2 - PX1*x1*x2*x3**2*y1*y3 + PX1*x1*x2**2*x3*y1*y3 + PX2*x1*x2*x3**2*y1*y2 - PX2*x1**2*x2*x3*y1*y2 + 2*PX1*x1*x2*x3*y2*y3**2 + 2*PX1*x1*x2*x3*y2**2*y3 + PX1*x1*x2*x3**2*y2*y3 + PX1*x1*x2**2*x3*y2*y3 - 2*PX1*x1**2*x2*x3*y2*y3 + 2*PX2*x1*x2*x3*y1*y3**2 + 2*PX2*x1*x2*x3*y1**2*y3 + PX2*x1*x2*x3**2*y1*y3 - 2*PX2*x1*x2**2*x3*y1*y3 + PX2*x1**2*x2*x3*y1*y3 + 2*PX3*x1*x2*x3*y1*y2**2 + 2*PX3*x1*x2*x3*y1**2*y2 - 2*PX3*x1*x2*x3**2*y1*y2 + PX3*x1*x2**2*x3*y1*y2 + PX3*x1**2*x2*x3*y1*y2 - PX2*x1*x2*x3**2*y2*y3 + PX2*x1**2*x2*x3*y2*y3 + PX3*x1*x2**2*x3*y1*y3 - PX3*x1**2*x2*x3*y1*y3 - PX3*x1*x2**2*x3*y2*y3 + PX3*x1**2*x2*x3*y2*y3 + PY1*x1*x2*x3*y1*y2**2 + PY1*x1*x2*x3*y1*y3**2 + PY2*x1*x2*x3*y1**2*y2 + PY1*x1*x2*x3*y2*y3**2 + PY1*x1*x2*x3*y2**2*y3 + PY2*x1*x2*x3*y1*y3**2 + PY2*x1*x2*x3*y1**2*y3 + PY3*x1*x2*x3*y1*y2**2 + PY3*x1*x2*x3*y1**2*y2 + PY2*x1*x2*x3*y2*y3**2 + PY3*x1*x2*x3*y1**2*y3 + PY3*x1*x2*x3*y2**2*y3 + PX1*x1*x2**2*y1*y2*y3 + PX1*x1*x3**2*y1*y2*y3 + PX2*x1**2*x2*y1*y2*y3 + PX1*x2*x3**2*y1*y2*y3 + PX1*x2**2*x3*y1*y2*y3 + PX2*x1*x3**2*y1*y2*y3 + PX2*x1**2*x3*y1*y2*y3 + PX3*x1*x2**2*y1*y2*y3 + PX3*x1**2*x2*y1*y2*y3 + PX2*x2*x3**2*y1*y2*y3 + PX3*x1**2*x3*y1*y2*y3 + PX3*x2**2*x3*y1*y2*y3 + PY1*x1*x2*y1*y2*y3**2 - PY1*x1*x2*y1*y2**2*y3 - PY1*x1*x3*y1*y2*y3**2 + PY1*x1*x3*y1*y2**2*y3 + PY2*x1*x2*y1*y2*y3**2 - PY2*x1*x2*y1**2*y2*y3 + PY1*x2*x3*y1*y2*y3**2 + PY1*x2*x3*y1*y2**2*y3 - 2*PY1*x2*x3*y1**2*y2*y3 + 2*PY1*x2*x3**2*y1*y2*y3 + 2*PY1*x2**2*x3*y1*y2*y3 + PY2*x1*x3*y1*y2*y3**2 - 2*PY2*x1*x3*y1*y2**2*y3 + PY2*x1*x3*y1**2*y2*y3 + 2*PY2*x1*x3**2*y1*y2*y3 + 2*PY2*x1**2*x3*y1*y2*y3 - 2*PY3*x1*x2*y1*y2*y3**2 + PY3*x1*x2*y1*y2**2*y3 + PY3*x1*x2*y1**2*y2*y3 + 2*PY3*x1*x2**2*y1*y2*y3 + 2*PY3*x1**2*x2*y1*y2*y3 - PY2*x2*x3*y1*y2*y3**2 + PY2*x2*x3*y1**2*y2*y3 + PY3*x1*x3*y1*y2**2*y3 - PY3*x1*x3*y1**2*y2*y3 - PY3*x2*x3*y1*y2**2*y3 + PY3*x2*x3*y1**2*y2*y3 - 4*PX1*x1*x2*x3*y1*y2*y3 - 4*PX2*x1*x2*x3*y1*y2*y3 - 4*PX3*x1*x2*x3*y1*y2*y3 - 4*PY1*x1*x2*x3*y1*y2*y3 - 4*PY2*x1*x2*x3*y1*y2*y3 - 4*PY3*x1*x2*x3*y1*y2*y3)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2 + x1*x2*y1 - x1*x2*y2 - x1*x3*y1 + x1*x3*y3 + x2*x3*y2 - x2*x3*y3 - x1*y1*y2 + x1*y1*y3 + x2*y1*y2 - x2*y2*y3 - x3*y1*y3 + x3*y2*y3))
#
#ansatzCoefficients[1] = (PX1*x2*y3 - PX1*x3*y2 - PX2*x1*y3 + PX2*x3*y1 + PX3*x1*y2 - PX3*x2*y1)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)
#
#ansatzCoefficients[2] = (PY1*x2*y3 - PY1*x3*y2 - PY2*x1*y3 + PY2*x3*y1 + PY3*x1*y2 - PY3*x2*y1)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)
#
#ansatzCoefficients[3] = (PX1*y2 - PX2*y1 - PX1*y3 + PX3*y1 + PX2*y3 - PX3*y2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2))
#
#ansatzCoefficients[4] = -(PX1*x2 - PX2*x1 - PX1*x3 + PX3*x1 + PX2*x3 - PX3*x2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)
#
#ansatzCoefficients[5] = (PY1*y2 - PY2*y1 - PY1*y3 + PY3*y1 + PY2*y3 - PY3*y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)
#
#ansatzCoefficients[6] = -(PY1*x2 - PY2*x1 - PY1*x3 + PY3*x1 + PY2*x3 - PY3*x2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2))
#
#ansatzCoefficients[7] = (2*W1*x1*y2**2 + 2*W1*x1*y3**2 + 2*W2*x2*y1**2 - 2*W1*x2*y3**2 - 2*W1*x3*y2**2 - 2*W2*x1*y3**2 - 2*W2*x3*y1**2 - 2*W3*x1*y2**2 - 2*W3*x2*y1**2 + 2*W2*x2*y3**2 + 2*W3*x3*y1**2 + 2*W3*x3*y2**2 + PX1*x1**2*y2**2 + PX1*x1**2*y3**2 + PX2*x2**2*y1**2 + PX1*x2**2*y3**2 + PX1*x3**2*y2**2 + PX2*x1**2*y3**2 + PX2*x3**2*y1**2 + PX3*x1**2*y2**2 + PX3*x2**2*y1**2 + PX2*x2**2*y3**2 + PX3*x3**2*y1**2 + PX3*x3**2*y2**2 - 2*W1*x2*y1*y2 - 2*W2*x1*y1*y2 - 4*W1*x1*y2*y3 + 2*W1*x2*y1*y3 + 2*W1*x3*y1*y2 + 2*W2*x1*y1*y3 + 2*W3*x1*y1*y2 + 2*W1*x2*y2*y3 - 2*W1*x3*y1*y3 + 2*W2*x1*y2*y3 - 4*W2*x2*y1*y3 + 2*W2*x3*y1*y2 - 2*W3*x1*y1*y3 + 2*W3*x2*y1*y2 + 2*W1*x3*y2*y3 + 2*W2*x3*y1*y3 + 2*W3*x1*y2*y3 + 2*W3*x2*y1*y3 - 4*W3*x3*y1*y2 - 2*W2*x3*y2*y3 - 2*W3*x2*y2*y3 - 2*PX1*x1*x2*y3**2 - 2*PX1*x1*x3*y2**2 - 2*PX2*x1*x2*y3**2 - 2*PX2*x2*x3*y1**2 - 2*PX3*x1*x3*y2**2 - 2*PX3*x2*x3*y1**2 - 2*PX1*x1**2*y2*y3 - PX1*x2**2*y1*y3 - PX1*x3**2*y1*y2 - PX2*x1**2*y2*y3 - 2*PX2*x2**2*y1*y3 - PX2*x3**2*y1*y2 - PX3*x1**2*y2*y3 - PX3*x2**2*y1*y3 - 2*PX3*x3**2*y1*y2 + PY1*x1*y1*y2**2 + PY1*x1*y1*y3**2 - PY1*x2*y1**2*y2 - PY2*x1*y1*y2**2 - PY1*x2*y1*y3**2 + PY1*x2*y1**2*y3 - PY1*x3*y1*y2**2 + PY1*x3*y1**2*y2 + PY2*x2*y1**2*y2 - PY1*x3*y1**2*y3 - PY2*x1*y2*y3**2 + PY2*x1*y2**2*y3 + PY2*x3*y1*y2**2 - PY2*x3*y1**2*y2 - PY3*x1*y1*y3**2 + PY2*x2*y2*y3**2 + PY3*x1*y2*y3**2 - PY3*x1*y2**2*y3 + PY3*x2*y1*y3**2 - PY3*x2*y1**2*y3 - PY2*x3*y2**2*y3 - PY3*x2*y2*y3**2 + PY3*x3*y1**2*y3 + PY3*x3*y2**2*y3 - 2*W1*x1*x2*y2**2 - 2*W2*x1*x2*y1**2 + 2*W2*x1*x3*y1**2 + 2*W3*x1*x2*y1**2 - 2*W1*x1*x3*y3**2 + 2*W1*x2*x3*y2**2 + 2*W3*x1*x2*y2**2 - 2*W3*x1*x3*y1**2 + 2*W1*x2*x3*y3**2 + 2*W2*x1*x3*y3**2 - 2*W2*x2*x3*y3**2 - 2*W3*x2*x3*y2**2 + 2*W1*x2**2*y1*y2 + 2*W2*x1**2*y1*y2 - 2*W2*x1**2*y1*y3 - 2*W3*x1**2*y1*y2 - 2*W1*x2**2*y2*y3 + 2*W1*x3**2*y1*y3 + 2*W3*x1**2*y1*y3 - 2*W3*x2**2*y1*y2 - 2*W1*x3**2*y2*y3 - 2*W2*x3**2*y1*y3 + 2*W2*x3**2*y2*y3 + 2*W3*x2**2*y2*y3 - PX1*x1**2*x2*y2**2 - PX2*x1*x2**2*y1**2 - PX2*x1*x3**2*y1**2 - PX3*x1*x2**2*y1**2 - PX1*x2*x3**2*y2**2 - PX1*x1**2*x3*y3**2 - PX3*x1*x3**2*y1**2 - PX3*x1**2*x2*y2**2 - PX1*x2**2*x3*y3**2 - PX2*x1**2*x3*y3**2 - PX2*x2**2*x3*y3**2 - PX3*x2*x3**2*y2**2 + PY1*x2**2*y1**2*y2 + PY2*x1**2*y1*y2**2 + PY2*x1**2*y1*y3**2 + PY3*x1**2*y1*y2**2 + PY1*x2**2*y2*y3**2 + PY1*x3**2*y1**2*y3 + PY3*x1**2*y1*y3**2 + PY3*x2**2*y1**2*y2 + PY1*x3**2*y2**2*y3 + PY2*x3**2*y1**2*y3 + PY2*x3**2*y2**2*y3 + PY3*x2**2*y2*y3**2 - PX1*x1*x2*y1*y2 + PX1*x1*x2*y1*y3 + PX1*x1*x3*y1*y2 - PX2*x1*x2*y1*y2 + 2*PX1*x1*x2*y2*y3 - PX1*x1*x3*y1*y3 + PX1*x2*x3*y1*y2 + 2*PX2*x1*x2*y1*y3 + PX2*x1*x3*y1*y2 - 2*PX3*x1*x2*y1*y2 + 2*PX1*x1*x3*y2*y3 + PX1*x2*x3*y1*y3 + PX2*x1*x2*y2*y3 - 2*PX2*x1*x3*y1*y3 + PX2*x2*x3*y1*y2 + PX3*x1*x2*y1*y3 + 2*PX3*x1*x3*y1*y2 - 2*PX1*x2*x3*y2*y3 + PX2*x1*x3*y2*y3 + 2*PX2*x2*x3*y1*y3 + PX3*x1*x2*y2*y3 - PX3*x1*x3*y1*y3 + 2*PX3*x2*x3*y1*y2 - PX2*x2*x3*y2*y3 + PX3*x1*x3*y2*y3 + PX3*x2*x3*y1*y3 - PX3*x2*x3*y2*y3 - 2*PY1*x1*y1*y2*y3 + PY1*x2*y1*y2*y3 + PY2*x1*y1*y2*y3 + PY1*x3*y1*y2*y3 - 2*PY2*x2*y1*y2*y3 + PY3*x1*y1*y2*y3 + PY2*x3*y1*y2*y3 + PY3*x2*y1*y2*y3 - 2*PY3*x3*y1*y2*y3 + 2*W1*x1*x2*y2*y3 - 2*W1*x2*x3*y1*y2 + 2*W2*x1*x2*y1*y3 - 2*W2*x1*x3*y1*y2 + 2*W1*x1*x3*y2*y3 - 2*W1*x2*x3*y1*y3 - 2*W3*x1*x2*y1*y3 + 2*W3*x1*x3*y1*y2 - 2*W2*x1*x3*y2*y3 + 2*W2*x2*x3*y1*y3 - 2*W3*x1*x2*y2*y3 + 2*W3*x2*x3*y1*y2 + 2*PX1*x1*x2*x3*y2**2 + 2*PX2*x1*x2*x3*y1**2 + 2*PX1*x1*x2*x3*y3**2 + 2*PX3*x1*x2*x3*y1**2 + 2*PX2*x1*x2*x3*y3**2 + 2*PX3*x1*x2*x3*y2**2 + PX1*x1*x2**2*y1*y2 + PX1*x1*x2**2*y1*y3 + PX1*x1*x3**2*y1*y2 + PX2*x1**2*x2*y1*y2 - 2*PX1*x1*x2**2*y2*y3 + PX1*x1*x3**2*y1*y3 + PX1*x1**2*x2*y2*y3 + PX2*x1*x2**2*y1*y3 - 2*PX2*x1**2*x2*y1*y3 + PX3*x1*x2**2*y1*y2 + PX3*x1**2*x2*y1*y2 - 2*PX1*x1*x3**2*y2*y3 + PX1*x1**2*x3*y2*y3 + PX2*x1*x3**2*y1*y3 + PX2*x2*x3**2*y1*y2 + PX2*x1**2*x2*y2*y3 + PX2*x1**2*x3*y1*y3 + PX3*x1*x3**2*y1*y2 - 2*PX3*x1**2*x3*y1*y2 + PX1*x2*x3**2*y2*y3 + PX1*x2**2*x3*y2*y3 - 2*PX2*x2*x3**2*y1*y3 + PX2*x2**2*x3*y1*y3 + PX3*x2*x3**2*y1*y2 + PX3*x1**2*x3*y1*y3 - 2*PX3*x2**2*x3*y1*y2 + PX2*x2*x3**2*y2*y3 + PX3*x1**2*x3*y2*y3 + PX3*x2**2*x3*y1*y3 + PX3*x2**2*x3*y2*y3 - PY1*x1*x2*y1*y2**2 - PY1*x1*x2*y1*y3**2 - PY1*x1*x3*y1*y2**2 - PY2*x1*x2*y1**2*y2 - PY1*x1*x3*y1*y3**2 + 2*PY1*x2*x3*y1*y2**2 - PY1*x2*x3*y1**2*y2 - PY2*x1*x3*y1*y2**2 + 2*PY2*x1*x3*y1**2*y2 - PY3*x1*x2*y1*y2**2 - PY3*x1*x2*y1**2*y2 + 2*PY1*x2*x3*y1*y3**2 - PY1*x2*x3*y1**2*y3 - PY2*x1*x2*y2*y3**2 - PY2*x1*x3*y1*y3**2 - PY2*x1*x3*y1**2*y3 - PY2*x2*x3*y1**2*y2 - PY3*x1*x2*y1*y3**2 + 2*PY3*x1*x2*y1**2*y3 - PY1*x2*x3*y2*y3**2 - PY1*x2*x3*y2**2*y3 + 2*PY2*x1*x3*y2*y3**2 - PY2*x1*x3*y2**2*y3 - PY3*x1*x2*y2*y3**2 + 2*PY3*x1*x2*y2**2*y3 - PY3*x1*x3*y1**2*y3 - PY2*x2*x3*y2*y3**2 - PY3*x1*x3*y2**2*y3 - PY3*x2*x3*y1**2*y3 - PY3*x2*x3*y2**2*y3 - 2*PY1*x2**2*y1*y2*y3 - 2*PY2*x1**2*y1*y2*y3 - 2*PY1*x3**2*y1*y2*y3 - 2*PY3*x1**2*y1*y2*y3 - 2*PY2*x3**2*y1*y2*y3 - 2*PY3*x2**2*y1*y2*y3 - 2*PX1*x1*x2*x3*y1*y2 - 2*PX1*x1*x2*x3*y1*y3 - 2*PX2*x1*x2*x3*y1*y2 - 2*PX2*x1*x2*x3*y2*y3 - 2*PX3*x1*x2*x3*y1*y3 - 2*PX3*x1*x2*x3*y2*y3 + 2*PY1*x1*x2*y1*y2*y3 + 2*PY1*x1*x3*y1*y2*y3 + 2*PY2*x1*x2*y1*y2*y3 + 2*PY2*x2*x3*y1*y2*y3 + 2*PY3*x1*x3*y1*y2*y3 + 2*PY3*x2*x3*y1*y2*y3)/(2*(x1**2*x2*y1*y2 - x1**2*x2*y1*y3 - x1**2*x2*y2**2 + x1**2*x2*y2*y3 - x1**2*x3*y1*y2 + x1**2*x3*y1*y3 + x1**2*x3*y2*y3 - x1**2*x3*y3**2 - x1**2*y1*y2**2 + 2*x1**2*y1*y2*y3 - x1**2*y1*y3**2 + x1**2*y2**2 - 2*x1**2*y2*y3 + x1**2*y3**2 - x1*x2**2*y1**2 + x1*x2**2*y1*y2 + x1*x2**2*y1*y3 - x1*x2**2*y2*y3 + 2*x1*x2*x3*y1**2 - 2*x1*x2*x3*y1*y2 - 2*x1*x2*x3*y1*y3 + 2*x1*x2*x3*y2**2 - 2*x1*x2*x3*y2*y3 + 2*x1*x2*x3*y3**2 + x1*x2*y1**2*y2 - x1*x2*y1**2*y3 + x1*x2*y1*y2**2 - 2*x1*x2*y1*y2*y3 - 2*x1*x2*y1*y2 + x1*x2*y1*y3**2 + 2*x1*x2*y1*y3 - x1*x2*y2**2*y3 + x1*x2*y2*y3**2 + 2*x1*x2*y2*y3 - 2*x1*x2*y3**2 - x1*x3**2*y1**2 + x1*x3**2*y1*y2 + x1*x3**2*y1*y3 - x1*x3**2*y2*y3 - x1*x3*y1**2*y2 + x1*x3*y1**2*y3 + x1*x3*y1*y2**2 - 2*x1*x3*y1*y2*y3 + 2*x1*x3*y1*y2 + x1*x3*y1*y3**2 - 2*x1*x3*y1*y3 + x1*x3*y2**2*y3 - 2*x1*x3*y2**2 - x1*x3*y2*y3**2 + 2*x1*x3*y2*y3 - x2**2*x3*y1*y2 + x2**2*x3*y1*y3 + x2**2*x3*y2*y3 - x2**2*x3*y3**2 - x2**2*y1**2*y2 + x2**2*y1**2 + 2*x2**2*y1*y2*y3 - 2*x2**2*y1*y3 - x2**2*y2*y3**2 + x2**2*y3**2 + x2*x3**2*y1*y2 - x2*x3**2*y1*y3 - x2*x3**2*y2**2 + x2*x3**2*y2*y3 + x2*x3*y1**2*y2 + x2*x3*y1**2*y3 - 2*x2*x3*y1**2 - x2*x3*y1*y2**2 - 2*x2*x3*y1*y2*y3 + 2*x2*x3*y1*y2 - x2*x3*y1*y3**2 + 2*x2*x3*y1*y3 + x2*x3*y2**2*y3 + x2*x3*y2*y3**2 - 2*x2*x3*y2*y3 - x3**2*y1**2*y3 + x3**2*y1**2 + 2*x3**2*y1*y2*y3 - 2*x3**2*y1*y2 - x3**2*y2**2*y3 + x3**2*y2**2))
#
#ansatzCoefficients[8] = (2*W1*x2**2*y1 + 2*W1*x3**2*y1 + 2*W2*x1**2*y2 - 2*W1*x2**2*y3 - 2*W1*x3**2*y2 - 2*W2*x1**2*y3 - 2*W2*x3**2*y1 - 2*W3*x1**2*y2 - 2*W3*x2**2*y1 + 2*W2*x3**2*y2 + 2*W3*x1**2*y3 + 2*W3*x2**2*y3 + PY1*x2**2*y1**2 + PY1*x3**2*y1**2 + PY2*x1**2*y2**2 + PY1*x2**2*y3**2 + PY1*x3**2*y2**2 + PY2*x1**2*y3**2 + PY2*x3**2*y1**2 + PY3*x1**2*y2**2 + PY3*x2**2*y1**2 + PY2*x3**2*y2**2 + PY3*x1**2*y3**2 + PY3*x2**2*y3**2 - 2*W1*x1*x2*y2 - 2*W2*x1*x2*y1 + 2*W1*x1*x2*y3 + 2*W1*x1*x3*y2 - 4*W1*x2*x3*y1 + 2*W2*x1*x3*y1 + 2*W3*x1*x2*y1 - 2*W1*x1*x3*y3 + 2*W1*x2*x3*y2 + 2*W2*x1*x2*y3 - 4*W2*x1*x3*y2 + 2*W2*x2*x3*y1 + 2*W3*x1*x2*y2 - 2*W3*x1*x3*y1 + 2*W1*x2*x3*y3 + 2*W2*x1*x3*y3 - 4*W3*x1*x2*y3 + 2*W3*x1*x3*y2 + 2*W3*x2*x3*y1 - 2*W2*x2*x3*y3 - 2*W3*x2*x3*y2 + PX1*x1*x2**2*y1 + PX1*x1*x3**2*y1 - PX1*x1**2*x2*y2 - PX2*x1*x2**2*y1 - PX1*x1*x2**2*y3 - PX1*x1*x3**2*y2 + PX1*x1**2*x2*y3 + PX1*x1**2*x3*y2 + PX2*x1**2*x2*y2 - PX1*x1**2*x3*y3 + PX2*x1*x2**2*y3 - PX2*x2*x3**2*y1 - PX2*x1**2*x2*y3 + PX2*x2**2*x3*y1 - PX3*x1*x3**2*y1 + PX2*x2*x3**2*y2 + PX3*x1*x3**2*y2 + PX3*x2*x3**2*y1 - PX3*x1**2*x3*y2 - PX3*x2**2*x3*y1 - PX2*x2**2*x3*y3 - PX3*x2*x3**2*y2 + PX3*x1**2*x3*y3 + PX3*x2**2*x3*y3 - PY1*x1*x2*y3**2 - PY1*x1*x3*y2**2 - 2*PY1*x2*x3*y1**2 - PY2*x1*x2*y3**2 - 2*PY2*x1*x3*y2**2 - PY2*x2*x3*y1**2 - 2*PY3*x1*x2*y3**2 - PY3*x1*x3*y2**2 - PY3*x2*x3*y1**2 - 2*PY1*x2**2*y1*y3 - 2*PY1*x3**2*y1*y2 - 2*PY2*x1**2*y2*y3 - 2*PY2*x3**2*y1*y2 - 2*PY3*x1**2*y2*y3 - 2*PY3*x2**2*y1*y3 + 2*W1*x1*x2*y2**2 + 2*W2*x1*x2*y1**2 - 2*W2*x1*x3*y1**2 - 2*W3*x1*x2*y1**2 + 2*W1*x1*x3*y3**2 - 2*W1*x2*x3*y2**2 - 2*W3*x1*x2*y2**2 + 2*W3*x1*x3*y1**2 - 2*W1*x2*x3*y3**2 - 2*W2*x1*x3*y3**2 + 2*W2*x2*x3*y3**2 + 2*W3*x2*x3*y2**2 - 2*W1*x2**2*y1*y2 - 2*W2*x1**2*y1*y2 + 2*W2*x1**2*y1*y3 + 2*W3*x1**2*y1*y2 + 2*W1*x2**2*y2*y3 - 2*W1*x3**2*y1*y3 - 2*W3*x1**2*y1*y3 + 2*W3*x2**2*y1*y2 + 2*W1*x3**2*y2*y3 + 2*W2*x3**2*y1*y3 - 2*W2*x3**2*y2*y3 - 2*W3*x2**2*y2*y3 + PX1*x1**2*x2*y2**2 + PX2*x1*x2**2*y1**2 + PX2*x1*x3**2*y1**2 + PX3*x1*x2**2*y1**2 + PX1*x2*x3**2*y2**2 + PX1*x1**2*x3*y3**2 + PX3*x1*x3**2*y1**2 + PX3*x1**2*x2*y2**2 + PX1*x2**2*x3*y3**2 + PX2*x1**2*x3*y3**2 + PX2*x2**2*x3*y3**2 + PX3*x2*x3**2*y2**2 - PY1*x2**2*y1**2*y2 - PY2*x1**2*y1*y2**2 - PY2*x1**2*y1*y3**2 - PY3*x1**2*y1*y2**2 - PY1*x2**2*y2*y3**2 - PY1*x3**2*y1**2*y3 - PY3*x1**2*y1*y3**2 - PY3*x2**2*y1**2*y2 - PY1*x3**2*y2**2*y3 - PY2*x3**2*y1**2*y3 - PY2*x3**2*y2**2*y3 - PY3*x2**2*y2*y3**2 - 2*PX1*x1*x2*x3*y1 + PX1*x1*x2*x3*y2 + PX2*x1*x2*x3*y1 + PX1*x1*x2*x3*y3 - 2*PX2*x1*x2*x3*y2 + PX3*x1*x2*x3*y1 + PX2*x1*x2*x3*y3 + PX3*x1*x2*x3*y2 - 2*PX3*x1*x2*x3*y3 - PY1*x1*x2*y1*y2 + PY1*x1*x2*y1*y3 + PY1*x1*x3*y1*y2 - PY2*x1*x2*y1*y2 + PY1*x1*x2*y2*y3 - PY1*x1*x3*y1*y3 + 2*PY1*x2*x3*y1*y2 + PY2*x1*x2*y1*y3 + 2*PY2*x1*x3*y1*y2 - 2*PY3*x1*x2*y1*y2 + PY1*x1*x3*y2*y3 + 2*PY1*x2*x3*y1*y3 + PY2*x1*x2*y2*y3 - 2*PY2*x1*x3*y1*y3 + PY2*x2*x3*y1*y2 + 2*PY3*x1*x2*y1*y3 + PY3*x1*x3*y1*y2 - 2*PY1*x2*x3*y2*y3 + 2*PY2*x1*x3*y2*y3 + PY2*x2*x3*y1*y3 + 2*PY3*x1*x2*y2*y3 - PY3*x1*x3*y1*y3 + PY3*x2*x3*y1*y2 - PY2*x2*x3*y2*y3 + PY3*x1*x3*y2*y3 + PY3*x2*x3*y1*y3 - PY3*x2*x3*y2*y3 - 2*W1*x1*x2*y2*y3 + 2*W1*x2*x3*y1*y2 - 2*W2*x1*x2*y1*y3 + 2*W2*x1*x3*y1*y2 - 2*W1*x1*x3*y2*y3 + 2*W1*x2*x3*y1*y3 + 2*W3*x1*x2*y1*y3 - 2*W3*x1*x3*y1*y2 + 2*W2*x1*x3*y2*y3 - 2*W2*x2*x3*y1*y3 + 2*W3*x1*x2*y2*y3 - 2*W3*x2*x3*y1*y2 - 2*PX1*x1*x2*x3*y2**2 - 2*PX2*x1*x2*x3*y1**2 - 2*PX1*x1*x2*x3*y3**2 - 2*PX3*x1*x2*x3*y1**2 - 2*PX2*x1*x2*x3*y3**2 - 2*PX3*x1*x2*x3*y2**2 - PX1*x1*x2**2*y1*y2 - PX1*x1*x2**2*y1*y3 - PX1*x1*x3**2*y1*y2 - PX2*x1**2*x2*y1*y2 + 2*PX1*x1*x2**2*y2*y3 - PX1*x1*x3**2*y1*y3 - PX1*x1**2*x2*y2*y3 - PX2*x1*x2**2*y1*y3 + 2*PX2*x1**2*x2*y1*y3 - PX3*x1*x2**2*y1*y2 - PX3*x1**2*x2*y1*y2 + 2*PX1*x1*x3**2*y2*y3 - PX1*x1**2*x3*y2*y3 - PX2*x1*x3**2*y1*y3 - PX2*x2*x3**2*y1*y2 - PX2*x1**2*x2*y2*y3 - PX2*x1**2*x3*y1*y3 - PX3*x1*x3**2*y1*y2 + 2*PX3*x1**2*x3*y1*y2 - PX1*x2*x3**2*y2*y3 - PX1*x2**2*x3*y2*y3 + 2*PX2*x2*x3**2*y1*y3 - PX2*x2**2*x3*y1*y3 - PX3*x2*x3**2*y1*y2 - PX3*x1**2*x3*y1*y3 + 2*PX3*x2**2*x3*y1*y2 - PX2*x2*x3**2*y2*y3 - PX3*x1**2*x3*y2*y3 - PX3*x2**2*x3*y1*y3 - PX3*x2**2*x3*y2*y3 + PY1*x1*x2*y1*y2**2 + PY1*x1*x2*y1*y3**2 + PY1*x1*x3*y1*y2**2 + PY2*x1*x2*y1**2*y2 + PY1*x1*x3*y1*y3**2 - 2*PY1*x2*x3*y1*y2**2 + PY1*x2*x3*y1**2*y2 + PY2*x1*x3*y1*y2**2 - 2*PY2*x1*x3*y1**2*y2 + PY3*x1*x2*y1*y2**2 + PY3*x1*x2*y1**2*y2 - 2*PY1*x2*x3*y1*y3**2 + PY1*x2*x3*y1**2*y3 + PY2*x1*x2*y2*y3**2 + PY2*x1*x3*y1*y3**2 + PY2*x1*x3*y1**2*y3 + PY2*x2*x3*y1**2*y2 + PY3*x1*x2*y1*y3**2 - 2*PY3*x1*x2*y1**2*y3 + PY1*x2*x3*y2*y3**2 + PY1*x2*x3*y2**2*y3 - 2*PY2*x1*x3*y2*y3**2 + PY2*x1*x3*y2**2*y3 + PY3*x1*x2*y2*y3**2 - 2*PY3*x1*x2*y2**2*y3 + PY3*x1*x3*y1**2*y3 + PY2*x2*x3*y2*y3**2 + PY3*x1*x3*y2**2*y3 + PY3*x2*x3*y1**2*y3 + PY3*x2*x3*y2**2*y3 + 2*PY1*x2**2*y1*y2*y3 + 2*PY2*x1**2*y1*y2*y3 + 2*PY1*x3**2*y1*y2*y3 + 2*PY3*x1**2*y1*y2*y3 + 2*PY2*x3**2*y1*y2*y3 + 2*PY3*x2**2*y1*y2*y3 + 2*PX1*x1*x2*x3*y1*y2 + 2*PX1*x1*x2*x3*y1*y3 + 2*PX2*x1*x2*x3*y1*y2 + 2*PX2*x1*x2*x3*y2*y3 + 2*PX3*x1*x2*x3*y1*y3 + 2*PX3*x1*x2*x3*y2*y3 - 2*PY1*x1*x2*y1*y2*y3 - 2*PY1*x1*x3*y1*y2*y3 - 2*PY2*x1*x2*y1*y2*y3 - 2*PY2*x2*x3*y1*y2*y3 - 2*PY3*x1*x3*y1*y2*y3 - 2*PY3*x2*x3*y1*y2*y3)/(2*(x1**2*x2*y1*y2 - x1**2*x2*y1*y3 - x1**2*x2*y2**2 + x1**2*x2*y2*y3 - x1**2*x3*y1*y2 + x1**2*x3*y1*y3 + x1**2*x3*y2*y3 - x1**2*x3*y3**2 - x1**2*y1*y2**2 + 2*x1**2*y1*y2*y3 - x1**2*y1*y3**2 + x1**2*y2**2 - 2*x1**2*y2*y3 + x1**2*y3**2 - x1*x2**2*y1**2 + x1*x2**2*y1*y2 + x1*x2**2*y1*y3 - x1*x2**2*y2*y3 + 2*x1*x2*x3*y1**2 - 2*x1*x2*x3*y1*y2 - 2*x1*x2*x3*y1*y3 + 2*x1*x2*x3*y2**2 - 2*x1*x2*x3*y2*y3 + 2*x1*x2*x3*y3**2 + x1*x2*y1**2*y2 - x1*x2*y1**2*y3 + x1*x2*y1*y2**2 - 2*x1*x2*y1*y2*y3 - 2*x1*x2*y1*y2 + x1*x2*y1*y3**2 + 2*x1*x2*y1*y3 - x1*x2*y2**2*y3 + x1*x2*y2*y3**2 + 2*x1*x2*y2*y3 - 2*x1*x2*y3**2 - x1*x3**2*y1**2 + x1*x3**2*y1*y2 + x1*x3**2*y1*y3 - x1*x3**2*y2*y3 - x1*x3*y1**2*y2 + x1*x3*y1**2*y3 + x1*x3*y1*y2**2 - 2*x1*x3*y1*y2*y3 + 2*x1*x3*y1*y2 + x1*x3*y1*y3**2 - 2*x1*x3*y1*y3 + x1*x3*y2**2*y3 - 2*x1*x3*y2**2 - x1*x3*y2*y3**2 + 2*x1*x3*y2*y3 - x2**2*x3*y1*y2 + x2**2*x3*y1*y3 + x2**2*x3*y2*y3 - x2**2*x3*y3**2 - x2**2*y1**2*y2 + x2**2*y1**2 + 2*x2**2*y1*y2*y3 - 2*x2**2*y1*y3 - x2**2*y2*y3**2 + x2**2*y3**2 + x2*x3**2*y1*y2 - x2*x3**2*y1*y3 - x2*x3**2*y2**2 + x2*x3**2*y2*y3 + x2*x3*y1**2*y2 + x2*x3*y1**2*y3 - 2*x2*x3*y1**2 - x2*x3*y1*y2**2 - 2*x2*x3*y1*y2*y3 + 2*x2*x3*y1*y2 - x2*x3*y1*y3**2 + 2*x2*x3*y1*y3 + x2*x3*y2**2*y3 + x2*x3*y2*y3**2 - 2*x2*x3*y2*y3 - x3**2*y1**2*y3 + x3**2*y1**2 + 2*x3**2*y1*y2*y3 - 2*x3**2*y1*y2 - x3**2*y2**2*y3 + x3**2*y2**2))




#Go through and update everything "_u"
PHI_u = PHI.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

GAM_u = GAM.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

PD_u = PD.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

#C_u = C.subs([(a1,ansatzCoefficients[0]),(a2,ansatzCoefficients[1]),(a3,ansatzCoefficients[2]),(a4,ansatzCoefficients[3]),(a5,ansatzCoefficients[4]),(a6,ansatzCoefficients[5]),(a7,ansatzCoefficients[6]),(a8,ansatzCoefficients[7]),(a9,ansatzCoefficients[8])])

print("Update half complete")

#PX_u = PHI_u + PD_u
#PY_u = PHI_u - PD_u
GX_u = GAM_u + PD_u
GY_u = GAM_u - PD_u

# Section 3 ---------------------------------------------------
#px_u = sp.diff(PX_u,x)
#py_u = sp.diff(PY_u,y)
gx_u = sp.diff(GX_u,x)
gy_u = sp.diff(GY_u,y)

print("Update complete\n\nSymbolic B Matrix = ")

# Assemble B Matrix ------------------------------------------
#Bsym = sp.zeros(2,9)
#for col in range(9):
#    Bsym[0,col] = sp.diff(gx,UVector[col])
#    Bsym[1,col] = sp.diff(gy,UVector[col])
#
#sp.pprint(Bsym,wrap_line=False)


B = sp.zeros(2,9)
for col in range(9):
    B[0,col] = sp.diff(gx_u,UVector[col])
    B[1,col] = sp.diff(gy_u,UVector[col])
    
print("\nB Matrix for [0,0] [1,0] [0,1] = ")
sp.pprint(B,wrap_line=False)
#sp.pprint(B.subs([(x1,0),(y1,0),(x2,1.0),(y2,0.0),(x3,0.0),(y3,1.0)]),wrap_line=False)


# TESTING FOR SIMPLIFIED FORM

# gx and gy feed into B
# only a5,a6,a8,a9 are in gx and gy
#a5simp = sp.zeros(9,1)
#a6simp = sp.zeros(9,1)
#a8simp = sp.zeros(9,1)
#a9simp = sp.zeros(9,1)
#
#for col in range(9):
#    a5simp[col] = sp.diff(ansatzCoefficients[4],UVector[col])
#    a6simp[col] = sp.diff(ansatzCoefficients[5],UVector[col])
#    a8simp[col] = sp.diff(ansatzCoefficients[7],UVector[col])
#    a9simp[col] = sp.diff(ansatzCoefficients[8],UVector[col])
#
#Btest = sp.zeros(2,9)
#for col in range(9):
#    Btest[0,col] = a8simp[col] + y*(0.5*a5simp[col] - 0.5*a6simp[col]) - y*(a8simp[col] + a9simp[col])
#    Btest[1,col] = a9simp[col] - x*(0.5*a5simp[col] - 0.5*a6simp[col]) - x*(a8simp[col] + a9simp[col])
#    
#print("\nTest B Matrix for [0,0] [1,0] [0,1] = ")
#sp.pprint(Btest.subs([(x1,0),(y1,0),(x2,1.0),(y2,0.0),(x3,0.0),(y3,1.0)]),wrap_line=False)
#
#print("\n\n\nPrinting simplified forms of ansatz coefficients!\n")
#for col in range(9):
#    print("a5[",col,"] = ", sp.factor(a5simp[col]))
#    
#print("\n\n")    
#for col in range(9):
#    print("a6[",col,"] = ", sp.factor(a6simp[col]))
#    
#print("\n\n")       
#for col in range(9):
#    print("a8[",col,"] = ", sp.factor(a8simp[col]))
#    
#print("\n\n")       
#for col in range(9):
#    print("a9[",col,"] = ", sp.factor(a9simp[col]))