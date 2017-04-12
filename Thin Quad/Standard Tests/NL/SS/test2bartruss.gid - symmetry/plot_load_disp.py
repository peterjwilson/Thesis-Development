# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 21:30:56 2016

@author: sautt
"""

import matplotlib.pyplot as plt
import numpy as np

load = []
disp = []

with open("load.txt", "r") as fol:
  for line in fol:
    load.append(abs(float(line)))      
fol.close()

with open("displacement.txt", "r") as fod:
  for line in fod:
    disp.append(abs(float(line)))         
fod.close()


#custom_edit
E = 210000000000
A =0.01
EA =E*A
x =2
y =1
L =np.sqrt(x**2+y**2)
L3 = L*L*L

v =0.00
dv =0.01
loadAN =[]
dispAN =[]

for i in range(220):
    loadAN.append((EA/(2*L3))*(v*v -2*y*v)*(v-y))
    dispAN.append(v)
    v = v + dv
    
    

plt.plot(disp,load,color = '#FFA500', linewidth=3.0, label = 'TRUSS3D2N',antialiased=True)
plt.plot(dispAN,loadAN,'--', color = '#696969', label = 'analytical',antialiased=True)
plt.legend()
plt.xlabel('displacement [m]')
plt.ylabel('load [N]')
plt.title('load-displacement')
plt.grid()
#plt.savefig('disp_v_truss.pdf')
plt.show()

