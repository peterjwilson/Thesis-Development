# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 21:30:56 2016

@author: wilson
"""

#===============================================
#
#       hinged cylindrical roof
#
#===============================================






import matplotlib.pyplot as plt
import numpy as np

myfontsize = 17
labelfontsize = 12
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



load = [0]
dispz = [0]
dispx=[0]



with open("load.txt", "r") as fol:
  for line in fol:
    load.append((float(line)))      
fol.close()

with open("displacementz.txt", "r") as fod:
  for line in fod:
    dispz.append((float(line)))         
fod.close()

with open("displacementx.txt", "r") as fod:
  for line in fod:
    dispx.append((float(line)))         
fod.close()

#https://www.researchgate.net/publication/222423070_Popular_benchmark_problems_for_geometric_nonlinear_analysis_of_shells
#load_ref = [0]
#disp_ref = [0]
#max_load_ref = 3000
#with open("load_benchmark.txt", "r") as fol:
#  for line in fol:
#    load_ref.append(abs(float(line))*max_load_ref)      
#fol.close()
#
#with open("displacement_benchmark.txt", "r") as fod:
#  for line in fod:
#    disp_ref.append(abs(float(line))/disp_scale)         
#fod.close()


#custom_edit

##77B5FE = nice blue

fig = plt.figure(1)

plt.plot(dispz,load, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linestyle='None', markerfacecolor= 'None', marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=2,frameon=False,fontsize=labelfontsize)
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlabel('DisplacementZ [m]')
plt.ylabel('Load [N]')
#plt.title('Load-displacement curve of hinged cylindrical roof')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Load_displacement_curve_hinged_cylindrical_roof.pdf')



fig = plt.figure(2)

plt.plot(dispx,load, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linestyle='None', markerfacecolor= 'None', marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=2,frameon=False,fontsize=labelfontsize)
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlabel('DisplacementX [m]')
plt.ylabel('Load [N]')
#plt.title('Load-displacement curve of hinged cylindrical roof')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Load_displacement_curve_hinged_cylindrical_roof.pdf')





plt.show()

