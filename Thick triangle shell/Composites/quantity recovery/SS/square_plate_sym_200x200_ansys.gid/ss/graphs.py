# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 14:35:48 2017

@author: Peter
"""

import matplotlib.pyplot as plt
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
#===============================================================================

kratos_x = []    
kratos_s_xz_middle = []
with open("kratos_s_xz.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split(' ')
      kratos_x.append(line[0])
      kratos_s_xz_middle.append(float(line[1]))   
      
      
ansys_x = []    
ansys_s_xz_middle = []
with open("ansys_s_xz.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      ansys_x.append(line[0])
      ansys_s_xz_middle.append(float(line[1]))         
      
      
      
         
kratos_s_xx_bot = []
with open("kratos_s_xx_bot.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split(' ')
      #kratos_x.append(line[0])
      kratos_s_xx_bot.append(float(line[1]))   
      
ansys_s_xx_bot = []
with open("ansys_s_xx_bot.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      #ansys_x.append(line[0])
      ansys_s_xx_bot.append(float(line[1]))       
      
      
      
      
      
      
      
# Stress XX bot surface
fig = plt.figure(30)
#plt.plot(kratos_x,kratos_s_xz_middle, color = '#77B5FE', linewidth=3.0, label = 'PYTHON',antialiased=True)
plt.plot(kratos_x,kratos_s_xx_bot, color = '#FF91A4', linewidth=3.0, label = 'GID',antialiased=True)
plt.plot(ansys_x,ansys_s_xx_bot,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'ANSYS',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate')
plt.ylabel('S_xx bot surface')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")      
      
      
      
      
      
      
      
      
      
      
# Stress XZ middle surface
fig = plt.figure(40)
#plt.plot(kratos_x,kratos_s_xz_middle, color = '#77B5FE', linewidth=3.0, label = 'PYTHON',antialiased=True)
plt.plot(kratos_x,kratos_s_xz_middle, color = '#FF91A4', linewidth=3.0, label = 'GID',antialiased=True)
plt.plot(ansys_x,ansys_s_xz_middle,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'ANSYS',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate')
plt.ylabel('S_xz middle surface')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")

















plt.show()