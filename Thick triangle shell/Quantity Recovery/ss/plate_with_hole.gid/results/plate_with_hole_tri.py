# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 21:30:56 2016

@author: wilson
"""


#===============================================
#
#       plate with hole - tri
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

r_ref = 10
s_yy_applied = 1
R = 1

r = []
s_yy = []
s_yy_ref = []

with open("s_yy.txt", "r") as fol:
  for line in fol:
    line = line.split(' ')
    r.append(abs((float(line[0]))-r_ref)) 
    s_yy.append((float(line[1])))
fol.close()

for i in range(len(r)):
    s_yy_ref.append(s_yy_applied/2*(1+R*R/r[i]/r[i] + 1 + 3*R**4/(r[i])**4))


s_xx = []
s_xx_ref = []   
with open("s_xx.txt", "r") as fol:
  for line in fol:
    line = line.split(' ')
    #r.append(abs((float(line[0]))-r_ref)) 
    s_xx.append((float(line[1])))
fol.close()

for i in range(len(r)):
    s_xx_ref.append(s_yy_applied/2*(1-R*R/r[i]/r[i]) 
    - s_yy_applied/2*(1 - 4*R*R/r[i]/r[i] + 3*R**4/(r[i])**4 ))  
    
    
    
    
    


DSG = '#FF91A4'             #salmon
DSGBasic = '#F5A352'        #orange
KRATOSTRI = '#BB7365'       #brown

fig = plt.figure(1)
plt.plot(r,s_yy, color = DSG, linewidth=2.0, 
         markersize = 7.0, marker='None', markeredgewidth = 2.0, markeredgecolor = DSG, 
         markerfacecolor= 'None', label = 'DSG',antialiased=True)

#plt.plot(dofs,disp_dsg_basic, color = DSGBasic, linewidth=2.0, 
#         markersize = 7.0, marker='o', markeredgewidth = 2.0, markeredgecolor = DSGBasic, 
#         markerfacecolor= 'None', label = 'Basic-T3',antialiased=True)

plt.plot(r,s_yy_ref,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(r,s_yy,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.plot(r,s_yy_ref, color = 'grey', linewidth=2.0, 
#         markersize = 7.0, marker='o', markeredgewidth = 2.0, markeredgecolor = 'grey', 
#         markerfacecolor= 'None', label = 'KRATOS T3',antialiased=True)


plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False,fontsize=labelfontsize+2)

plt.xlabel('Distance from hole center [m]',fontsize=myfontsize)
plt.ylabel('Stress_yy [Pa]',fontsize=myfontsize)
plt.grid()
plt.ylim([0,3.5])
plt.xlim([0,10])
plt.tick_params(labelsize=labelfontsize)
plt.savefig('plate_with_hole_results.pdf',bbox_inches="tight")



plt.show()