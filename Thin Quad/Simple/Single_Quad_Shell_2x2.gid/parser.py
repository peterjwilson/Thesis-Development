# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 19:58:51 2016

@author: Peter
"""

#===============================================================================
## Display which python version is used
print('Python version:') 
#import sys
#print(sys.version)

## Display which modules are imported
print('Imported:') 

import sympy as sp
print(sp) 
from sympy import Matrix

import numpy


def enforceRestraints(restraintVec,Kin,Fin):
    Kdim = Kin.shape[0]
    #numbers in restraintvec start from 1....
    for i in range(len(restraintVec)):
        currentDOF = restraintVec[i] - 1
        for j in range(Kdim):
            Kin[j,currentDOF] = 0
            Kin[currentDOF,j] = 0
        Kin[currentDOF,currentDOF] = 1
        Fin[currentDOF] = 0




result = numpy.loadtxt(open("Kthin.csv","rb"),delimiter=",",skiprows=0)
#print(result)
#print(result[0,0])

Kthin =  sp.zeros(24)
for row in range (24):
    for col in range (24):
        Kthin[row,col] = result[row,col]

sp.pprint(Kthin,wrap_line=False)

myRestrainedDOFs = [1,2,3,6,7,8,9,12,15,18,21,24]
force = sp.zeros(24,1)
force[13] = 0.5
force[19] = 0.5

sp.pprint(force,wrap_line=False)

enforceRestraints(myRestrainedDOFs,Kthin,force)

sp.pprint(Kthin,wrap_line=False)

#sp.pprint((Kthin**-1),wrap_line=False)

displacements = (Kthin**-1)*force

sp.pprint(displacements,wrap_line=False)
