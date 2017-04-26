# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:40:14 2017

@author: Peter
"""

#===============================================
#
#       simply supported plate analysis
#
#===============================================



import matplotlib.pyplot as plt
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
#===============================================================================

import numpy as np
import sympy as sp


def evaluateExpression(expression,x_in,y_in):
    return ((expression.subs(x_s,x_in)).subs(y_s,y_in))


# -----------------------------------------------------------------------------
#       READ KRATOS DISPLACEMENTS_Z
# -----------------------------------------------------------------------------
kratos_t = [0]
composite_thin_quad_vibrating_square = [0]
with open("composite_thin_quad_vibrating_square.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split(' ')
      kratos_t.append(float(line[0]))
      composite_thin_quad_vibrating_square.append(float(line[1]))
      

composite_thick_tri_vibrating_square = [0]
with open("composite_thick_tri_vibrating_square.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split(' ')
      #kratos_t.append(float(line[0]))
      composite_thick_tri_vibrating_square.append(float(line[1]))     
      
      
      
      
ref_clpt_omega = 0.03900058471
ref_clpt_w = []
for i in range(len(kratos_t)):
    ref_clpt_w.append(0.19/2*(-1+np.cos(ref_clpt_omega*kratos_t[i])))
    
    
time_strand = [0]
strand_composite_vibrating_lumped = [0]
counter = 0
with open("strand_composite_vibrating_lumped.txt", "r") as fol:
  for line in fol:
    counter += 1
    if counter%1 == 0:
        line = line.split('\t')
        time_strand.append((float(line[1])))    
        strand_composite_vibrating_lumped.append((float(line[2])))   
fol.close()

          
# -----------------------------------------------------------------------------
#       PLOT GRAPHS
## -----------------------------------------------------------------------------         
          
# DISPLACEMENTS ---------------------------------------------------------------
fig = plt.figure(1)
plt.plot(kratos_t,composite_thin_quad_vibrating_square, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ (Lumped)',antialiased=True)

plt.plot(kratos_t,composite_thick_tri_vibrating_square, color = '#FF91A4', linewidth=3.0, linestyle='--',label = 'DSG (Lumped)',antialiased=True)

plt.plot(time_strand,strand_composite_vibrating_lumped,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Strand7 (Lumped)',antialiased=True)

#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,150])
plt.xlabel('Time')
plt.ylabel('Transverse displacment')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('composite_vibrating_square.pdf',bbox_inches="tight")





plt.show()