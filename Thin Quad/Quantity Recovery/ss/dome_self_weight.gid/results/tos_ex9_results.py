# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 21:30:56 2016

@author: wilson
"""


#===============================================
#
#       tos ex9 results
#
#===============================================



import matplotlib.pyplot as plt
import numpy as np

myfontsize = 20
labelfontsize = 14
#===============================================================================
## Set style of plots
font = {'family' : 'sans-serif',
        'family' : 'Arial',
        'style'  : 'normal',
        'weight' : 'normal',
        'size'   : myfontsize}

plt.rc('font', **font)
plt.rc('font', serif='Helvetica Neue') 
#plt.rcParams['font.family'] = 'Arial'
#===============================================================================

t = 0.01

y = []
force_theta = []
force_phi_1 = []
force_phi_2 = []
force_phi = []

with open("tos_ex9_n_theta.txt", "r") as fol:
  for line in fol:
    line = line.split(' ')
    y.append((float(line[0]))) 
    force_theta.append((float(line[1]))*0.01)
fol.close()
#
#with open("tos_ex9_force_xx.txt", "r") as fol:
#  for line in fol:
#    line = line.split(' ')
#    force_phi_1.append((float(line[1])))
#fol.close()
#
#with open("tos_ex9_force_yy.txt", "r") as fol:
#  for line in fol:
#    line = line.split(' ')
#    force_phi_2.append((float(line[1])))
#fol.close()

#for i in range(len(force_phi_1)):
#    force_phi.append((force_phi_1[i]**2 + force_phi_2[i]**2)**0.5)

phi = []

for i in range(len((y))):
    phi.append(180*y[i]/np.pi/5)


#custom_edit




force_theta_ref = []
phi_ref = []
with open("ref_n_theta.txt", "r") as fol:
  for line in fol:
    line = line.split('\t')
    phi_ref.append((float(line[0]))) 
    force_theta_ref.append((float(line[1])))
fol.close()


##77B5FE = nice blue

plt.plot(phi,force_theta, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi,force_phi, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(phi_ref,force_theta_ref,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=2,frameon=False,fontsize=labelfontsize+2)
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlabel('Phi (degrees)')
plt.ylabel('n_phi [N]')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('Load_displacement_curve_open_cylinder_pullout.pdf')
plt.show()

