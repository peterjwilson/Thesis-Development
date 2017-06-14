# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 11:00:16 2017

@author: Peter
"""

import numpy as np
import sympy as sp

lambda_s = sp.Symbol('lambda')
l_s = sp.Symbol('l')

A = sp.zeros(4)

A[0,1] = 1
A[0,3] = 1

A[1,0] = sp.sin(lambda_s*l_s)
A[1,1] = sp.cos(lambda_s*l_s)
A[1,2] = l_s
A[1,3] = 1

A[2,0] = lambda_s
A[2,2] = 1

A[3,0] = lambda_s*sp.cos(lambda_s*l_s)
A[3,1] = -1*lambda_s*sp.sin(lambda_s*l_s)
A[3,2] = 1

sp.pprint(A)

detA = A.det()

print(detA)

#sol = sp.solve(detA,lambda_s)
#
#print(sol)