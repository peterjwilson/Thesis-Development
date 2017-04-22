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
#       ESTABLISH ANALYTICAL SOLUTION - ONLY FOR SQUARE PLATES!
# -----------------------------------------------------------------------------

#setup symbolic stuff
x_s = sp.Symbol('x')
y_s = sp.Symbol('y')

#problem constants
a = b = 200.0
q0 = 1e7

E = 2e11
Poisson = 0.3
G = E/2.0/(1+Poisson)
h = 10

# integrated bending material matrix
D = sp.zeros(3,3)
D[0,0]= 1.0
D[1,1]= 1.0
D[0,1]= Poisson
D[1,0] = Poisson
D[2,2]= (1-Poisson)/2.0
D*= (E*h**3)/(12.0*(1-Poisson**2))

# pure plane stress material matrix
C = sp.zeros(3,3)
C[0,0]= 1.0
C[1,1]= 1.0
C[0,1]= Poisson
C[1,0] = Poisson
C[2,2] = (1-Poisson) #strain is already halved in Kratos 
C*= E/(1-Poisson**2)

# solution values
Qmn = q0
dmn = (np.pi**4)/(b**4) * (D[0,0] + 2.0*(D[0,1] + 2.0*D[2,2]) + D[1,1])
Wmn = Qmn/dmn
alpha = np.pi/a
beta = np.pi/a

w = -1.0*Wmn*(sp.sin(alpha*x_s))*(sp.sin(beta*y_s)) #transverse disp field
phi_x = (-1.0*sp.diff(w,x_s))                       #director angles         
phi_y = -1.0*sp.diff(w,y_s)
kappa_x = sp.diff(phi_x,x_s)                        #curvatures
kappa_y = sp.diff(phi_y,y_s)
kappa_xy = 0.5*(sp.diff(phi_x,y_s) + sp.diff(phi_y,x_s))

M_xx = (D[0,0]*alpha**2 + D[0,1]*beta**2) * Wmn * (sp.sin(np.pi*x_s/a)) * (sp.sin(np.pi*y_s/b))
M_yy = (D[0,1]*alpha**2 + D[1,1]*beta**2) * Wmn * (sp.sin(np.pi*x_s/a)) * (sp.sin(np.pi*y_s/b))
M_xy = -2.0*alpha*beta*D[2,2]* Wmn * (sp.cos(np.pi*x_s/a)) * (sp.cos(np.pi*y_s/b))


# max central disp
print("Analytical w_max = ",evaluateExpression(w,a/2.0,b/2.0))
print("\n\n")

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
#       CALCULATE DISCRETE ANALYTICAL RESULTS
# -----------------------------------------------------------------------------
analytical_w = []
for i in range(len(kratos_w_x)):
    # Displacements
    analytical_w.append(evaluateExpression(w,kratos_w_x[i],100.0))





analytical_kappa_y =[] 
analytical_kappa_x =[]
analytical_kappa_xy =[]
analytical_phi_x = []
analytical_phi_y = []
analytical_M_xx = []
analytical_M_yy = []
analytical_M_xy = []
e_x =[]
e_y =[]
e_xy =[]
analytical_s_x_top = []
analytical_s_y_top = []
analytical_s_xy_top = []
analytical_s_x_bottom = []
analytical_s_y_bottom = []
analytical_s_xy_bottom = []
analytical_s_vm_bottom =[]

for i in range(len(kratos_x)):   
    # Curvatures
    analytical_kappa_x.append(evaluateExpression(kappa_x,kratos_x[i],kratos_y[i]))
    analytical_kappa_y.append(evaluateExpression(kappa_y,kratos_x[i],kratos_y[i]))
    analytical_kappa_xy.append(evaluateExpression(kappa_xy,kratos_x[i],kratos_y[i]))
    
    # Rotations
    analytical_phi_x.append(evaluateExpression(phi_x,kratos_x[i],kratos_y[i]))
    analytical_phi_y.append(evaluateExpression(phi_y,kratos_x[i],kratos_y[i]))
    
    #Moments
    analytical_M_xx.append(evaluateExpression(M_xx,kratos_x[i],kratos_y[i]))
    analytical_M_yy.append(evaluateExpression(M_yy,kratos_x[i],kratos_y[i]))
    analytical_M_xy.append(evaluateExpression(M_xy,kratos_x[i],kratos_y[i]))
    
    # Stress and strain - top surface
    e_x.append(h/2.0*analytical_kappa_x[i])
    e_y.append(h/2.0*analytical_kappa_y[i])
    e_xy.append(h/2.0*analytical_kappa_xy[i]) #already halved in analytic formula
    Mat_e = sp.zeros(3,1)
    Mat_e[0,0] = e_x[i]
    Mat_e[1,0] = e_y[i]
    Mat_e[2,0] = e_xy[i]
    Mat_s = C*Mat_e
    analytical_s_x_top.append(Mat_s[0,0])
    analytical_s_y_top.append(Mat_s[1,0])
    analytical_s_xy_top.append(Mat_s[2,0])
    
    # Stress and strain - bottom surface
    analytical_s_x_bottom.append(-1.0*analytical_s_x_top[i])
    analytical_s_y_bottom.append(-1.0*analytical_s_y_top[i])
    analytical_s_xy_bottom.append(-1.0*analytical_s_xy_top[i])
    
    #Von mises - bottom surface
    vm_temp = analytical_s_x_bottom[i]**2 - analytical_s_x_bottom[i]*analytical_s_y_bottom[i] + analytical_s_y_bottom[i]**2 + 3.0*(analytical_s_xy_bottom[i]**2)
    analytical_s_vm_bottom.append(vm_temp**0.5)  



# -----------------------------------------------------------------------------
#       PLOT GRAPHS
# -----------------------------------------------------------------------------         
          
# DISPLACEMENTS ---------------------------------------------------------------
fig = plt.figure(1)
plt.plot(kratos_w_x,kratos_w, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_w_x,analytical_w,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
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




### CURVATURES   ---------------------------------------------------------------

# Kappa X
fig = plt.figure(11)
plt.plot(kratos_x,kratos_kappa_x, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_kappa_x,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate')
plt.ylabel('Kappa X')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")


# Kappa y
fig = plt.figure(12)
plt.plot(kratos_x,kratos_kappa_y, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_kappa_y,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate')
plt.ylabel('Kappa Y')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")


# Kappa xy
fig = plt.figure(13)
plt.plot(kratos_x,kratos_kappa_xy, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_kappa_xy,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate')
plt.ylabel('Kappa XY')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")



# MOMENTS      ---------------------------------------------------------------

# M xx
fig = plt.figure(21)
plt.plot(kratos_x,kratos_M_xx, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_M_xx,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate')
plt.ylabel('M_xx')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")


# M YY
fig = plt.figure(22)
plt.plot(kratos_x,kratos_M_yy, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_M_yy,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate')
plt.ylabel('M_yy')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")


# M XY
fig = plt.figure(23)
plt.plot(kratos_x,kratos_M_xy, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_M_xy,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate')
plt.ylabel('M_xy')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")



# STRESSES      ---------------------------------------------------------------

# Stress XX top surface
fig = plt.figure(31)
plt.plot(kratos_x,kratos_s_xx_top, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_s_x_top,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X coordinate',fontsize=myfontsize)
plt.ylabel('Stress (XX) @ top surface',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('navier_plate_quad_s_xx_top.pdf',bbox_inches="tight")



# Stress YY top surface
fig = plt.figure(32)
plt.plot(kratos_x,kratos_s_yy_top, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_s_y_top,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate',fontsize=myfontsize)
plt.ylabel('S_yy top surface',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")


# Stress XY top surface
fig = plt.figure(33)
plt.plot(kratos_x,kratos_s_xy_top, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_s_xy_top,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X Coordinate',fontsize=myfontsize)
plt.ylabel('S_xy top surface',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")





# Stress XX bottom surface
fig = plt.figure(34)
plt.plot(kratos_x,kratos_s_xx_bottom, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_s_x_bottom,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X coordinate',fontsize=myfontsize)
plt.ylabel('Stress_XX @ bottom surface',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('navier_plate_quad_s_xx_bot.pdf',bbox_inches="tight")



# Stress YY bottom surface
fig = plt.figure(35)
plt.plot(kratos_x,kratos_s_yy_bottom, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_s_y_bottom,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X coordinate',fontsize=myfontsize)
plt.ylabel('Stress_YY @ bottom surface',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('navier_plate_quad_s_yy_bot.pdf',bbox_inches="tight")



# Stress XY bottom surface
fig = plt.figure(36)
plt.plot(kratos_x,kratos_s_xy_bottom, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_s_xy_bottom,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X coordinate',fontsize=myfontsize)
plt.ylabel('Stress_XY @ bottom surface',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('navier_plate_quad_s_xy_bot.pdf',bbox_inches="tight")




# Stress von mises bottom surface
fig = plt.figure(37)
plt.plot(kratos_x,kratos_s_vm_bottom, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_s_vm_bottom,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,100])
plt.xlabel('X coordinate',fontsize=myfontsize)
plt.ylabel('Von Mises stress @ bottom surface',fontsize=myfontsize)
plt.grid()
plt.tick_params(labelsize=labelfontsize)
plt.savefig('navier_plate_quad_s_vm_bot.pdf',bbox_inches="tight")





plt.show()