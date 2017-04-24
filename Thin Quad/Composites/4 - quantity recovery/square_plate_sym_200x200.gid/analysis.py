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
kratos_w = []
kratos_w_x = []
with open("displacements.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      kratos_w_x.append(float(line[0]))
      kratos_w.append(float(line[1]))
      
# -----------------------------------------------------------------------------
#       READ KRATOS ROTATION_Y
# -----------------------------------------------------------------------------
kratos_phi_y = []
with open("rotations.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      #kratos_w_x.append(float(line[0]))
      kratos_phi_y.append(float(line[1]))

# -----------------------------------------------------------------------------
#       READ KRATOS CURVATURES
# -----------------------------------------------------------------------------
kratos_x =[]
kratos_y =[]
kratos_kappa_x = []
kratos_kappa_y = []
kratos_kappa_xy = []

with open("curvatures.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      kratos_x.append(float(line[0]))
      kratos_y.append(float(line[1]))
      kratos_kappa_x.append(float(line[2]))
      kratos_kappa_y.append(float(line[6]))
      kratos_kappa_xy.append(float(line[3]))
      
# -----------------------------------------------------------------------------
#       READ KRATOS MOMENTS
# -----------------------------------------------------------------------------      
kratos_M_xx = []
kratos_M_yy = []
kratos_M_xy = []
with open("moments.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      kratos_M_xx.append(float(line[2]))
      kratos_M_yy.append(float(line[6]))
      kratos_M_xy.append(float(line[3]))    
      
# -----------------------------------------------------------------------------
#       READ KRATOS TOP SURFACE STRESSES
# -----------------------------------------------------------------------------      
kratos_s_xx_top = []
kratos_s_yy_top = []
kratos_s_xy_top = []
with open("stress_top_surface.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      kratos_s_xx_top.append(float(line[2]))
      kratos_s_yy_top.append(float(line[6])) 
      kratos_s_xy_top.append(float(line[3]))

# -----------------------------------------------------------------------------
#       READ KRATOS BOTTOM SURFACE STRESSES
# -----------------------------------------------------------------------------      
kratos_s_xx_bottom = []
kratos_s_yy_bottom = []
kratos_s_xy_bottom = []
with open("stress_bottom_surface.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      kratos_s_xx_bottom.append(float(line[2]))
      kratos_s_yy_bottom.append(float(line[6]))    
      kratos_s_xy_bottom.append(float(line[3]))
     
kratos_s_vm_bottom = []
with open("stress_vm_bottom_surface_file.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      kratos_s_vm_bottom.append(float(line[2]))   
      
      
# -----------------------------------------------------------------------------
#       READ KRATOS THRU THICKNESS RESULTS
# -----------------------------------------------------------------------------  

kratos_thru_z = [0.5,0.25,0.25,0,0,-0.25,-0.25,-0.5]
ticks =[-0.5,-0.25,0,0.25,0.5]

kratos_thru_s_xx = []
with open("laminate_thru_s_xx.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      for i in range(8):
          kratos_thru_s_xx.append(float(line[2+i]))  

kratos_thru_s_yy = []
with open("laminate_thru_s_yy.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      for i in range(8):
          kratos_thru_s_yy.append(float(line[2+i])) 
 
kratos_thru_s_xy = []
with open("laminate_thru_s_xy.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      for i in range(8):
          kratos_thru_s_xy.append(float(line[2+i]))           


kratos_thru_e_xx = []
with open("laminate_thru_e_xx.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      for i in range(8):
          kratos_thru_e_xx.append(float(line[2+i]))  

kratos_thru_e_yy = []
with open("laminate_thru_e_yy.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      for i in range(8):
          kratos_thru_e_yy.append(float(line[2+i])) 
 
kratos_thru_e_xy = []
with open("laminate_thru_e_xy.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      for i in range(8):
          kratos_thru_e_xy.append(float(line[2+i])) 
          
# -----------------------------------------------------------------------------
#       PLOT GRAPHS
# -----------------------------------------------------------------------------         
          
# DISPLACEMENTS ---------------------------------------------------------------
fig = plt.figure(1)
plt.plot(kratos_w_x,kratos_w, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_w_x,analytical_w,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate')
plt.ylabel('Transverse displacment')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")




#### CURVATURES   ---------------------------------------------------------------
#
## Kappa X
#fig = plt.figure(11)
#plt.plot(kratos_x,kratos_kappa_x, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_x,analytical_kappa_x,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X Coordinate')
#plt.ylabel('Kappa X')
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")
#
#
## Kappa y
#fig = plt.figure(12)
#plt.plot(kratos_x,kratos_kappa_y, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_x,analytical_kappa_y,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X Coordinate')
#plt.ylabel('Kappa Y')
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")
#
#
## Kappa xy
#fig = plt.figure(13)
#plt.plot(kratos_x,kratos_kappa_xy, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_x,analytical_kappa_xy,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X Coordinate')
#plt.ylabel('Kappa XY')
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")
#
#
#
## MOMENTS      ---------------------------------------------------------------
#
## M xx
#fig = plt.figure(21)
#plt.plot(kratos_x,kratos_M_xx, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_x,analytical_M_xx,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X Coordinate')
#plt.ylabel('M_xx')
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")
#
#
## M YY
#fig = plt.figure(22)
#plt.plot(kratos_x,kratos_M_yy, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_x,analytical_M_yy,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X Coordinate')
#plt.ylabel('M_yy')
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")
#
#
## M XY
#fig = plt.figure(23)
#plt.plot(kratos_x,kratos_M_xy, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_x,analytical_M_xy,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X Coordinate')
#plt.ylabel('M_xy')
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")



# STRESSES      ---------------------------------------------------------------

## Stress XX top surface
#fig = plt.figure(31)
#plt.plot(kratos_x,kratos_s_xx_top, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
##plt.plot(kratos_x,analytical_s_x_top,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X coordinate',fontsize=myfontsize)
#plt.ylabel('Stress (XX) @ top surface',fontsize=myfontsize)
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('navier_plate_quad_s_xx_top.pdf',bbox_inches="tight")
#
#
#
## Stress YY top surface
#fig = plt.figure(32)
#plt.plot(kratos_x,kratos_s_yy_top, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
##plt.plot(kratos_x,analytical_s_y_top,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X Coordinate',fontsize=myfontsize)
#plt.ylabel('S_yy top surface',fontsize=myfontsize)
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")
#
#
## Stress XY top surface
#fig = plt.figure(33)
#plt.plot(kratos_x,kratos_s_xy_top, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
##plt.plot(kratos_x,analytical_s_xy_top,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X Coordinate',fontsize=myfontsize)
#plt.ylabel('S_xy top surface',fontsize=myfontsize)
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")





## Stress XX bottom surface
#fig = plt.figure(34)
#plt.plot(kratos_x,kratos_s_xx_bottom, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
##plt.plot(kratos_x,analytical_s_x_bottom,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X coordinate',fontsize=myfontsize)
#plt.ylabel('Stress_XX @ bottom surface',fontsize=myfontsize)
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('navier_plate_quad_s_xx_bot.pdf',bbox_inches="tight")
#
#
#
## Stress YY bottom surface
#fig = plt.figure(35)
#plt.plot(kratos_x,kratos_s_yy_bottom, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
##plt.plot(kratos_x,analytical_s_y_bottom,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X coordinate',fontsize=myfontsize)
#plt.ylabel('Stress_YY @ bottom surface',fontsize=myfontsize)
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('navier_plate_quad_s_yy_bot.pdf',bbox_inches="tight")
#
#
#
## Stress XY bottom surface
#fig = plt.figure(36)
#plt.plot(kratos_x,kratos_s_xy_bottom, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
##plt.plot(kratos_x,analytical_s_xy_bottom,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X coordinate',fontsize=myfontsize)
#plt.ylabel('Stress_XY @ bottom surface',fontsize=myfontsize)
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('navier_plate_quad_s_xy_bot.pdf',bbox_inches="tight")




## Stress von mises bottom surface
#fig = plt.figure(37)
#plt.plot(kratos_x,kratos_s_vm_bottom, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
##plt.plot(kratos_x,analytical_s_vm_bottom,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
#plt.xlim([0,100])
#plt.xlabel('X coordinate',fontsize=myfontsize)
#plt.ylabel('Von Mises stress @ bottom surface',fontsize=myfontsize)
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('navier_plate_quad_s_vm_bot.pdf',bbox_inches="tight")






# Stress XX thru thickness @ plate center
fig = plt.figure(41)

ref_z = [0.5,0,-0.5]
ref_val =[21548,0,-21548]
plt.plot(kratos_thru_s_xx,kratos_thru_z, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(ref_val,ref_z,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.ylim([-0.75,0.75])
plt.yticks(ticks)
plt.ylabel('Z coordinate = z/h',fontsize=myfontsize)
plt.xlabel('Stress_XX',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('navier_plate_composite_quad_s_xx.pdf',bbox_inches="tight")




# Stress YY thru thickness @ plate center
fig = plt.figure(42)
ref_z = [0.25,0,-0.25]
ref_val =[10668,0,-10668]
plt.plot(kratos_thru_s_yy,kratos_thru_z, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(ref_val,ref_z,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.ylim([-0.75,0.75])
plt.yticks(ticks)
plt.ylabel('Z coordinate = z/h',fontsize=myfontsize)
plt.xlabel('Stress_YY',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('navier_plate_composite_quad_s_yy.pdf',bbox_inches="tight")



# Stress XY thru thickness @ plate center
fig = plt.figure(43)
ref_z = [0.5,0,-0.5]
ref_val =[-852,0,852]
plt.plot(kratos_thru_s_xy,kratos_thru_z, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(ref_val,ref_z,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.ylim([-0.75,0.75])
plt.yticks(ticks)
plt.ylabel('Z coordinate = z/h',fontsize=myfontsize)
plt.xlabel('Stress_XY',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('navier_plate_composite_quad_s_xy.pdf',bbox_inches="tight")



# Strain XX thru thickness @ plate center
fig = plt.figure(51)
plt.plot(kratos_thru_e_xx,kratos_thru_z, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_x,analytical_s_x_top,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.ylim([-0.75,0.75])
plt.yticks(ticks)
plt.ylabel('Z coordinate = z/h',fontsize=myfontsize)
plt.xlabel('Strain_XX',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('navier_plate_composite_quad_e_xx.pdf',bbox_inches="tight")


# Strain YY thru thickness @ plate center
fig = plt.figure(52)
plt.plot(kratos_thru_e_yy,kratos_thru_z, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_x,analytical_s_x_top,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.ylim([-0.75,0.75])
plt.yticks(ticks)
plt.ylabel('Z coordinate = z/h',fontsize=myfontsize)
plt.xlabel('Strain_YY',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('navier_plate_composite_quad_e_yy.pdf',bbox_inches="tight")


# Strain XY thru thickness @ plate center
fig = plt.figure(53)
plt.plot(kratos_thru_e_xy,kratos_thru_z, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_x,analytical_s_x_top,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.ylim([-0.75,0.75])
plt.yticks(ticks)
plt.ylabel('Z coordinate = z/h',fontsize=myfontsize)
plt.xlabel('Strain_XY',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('navier_plate_composite_quad_e_xy.pdf',bbox_inches="tight")



plt.show()