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
        'size'   : myfontsize}

plt.rc('font', **font)
plt.rc('font', serif='Helvetica Neue') 
#===============================================================================








import numpy as np
import sympy as sp
from copy import copy, deepcopy
#from sp import *

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
kappa_xy = sp.diff(phi_x,y_s) + sp.diff(phi_y,x_s)



# max central disp
print("Analytical w_max = ",evaluateExpression(w,a/2.0,b/2.0))
print("\n\n")
print("phi_x w_max = ",evaluateExpression(phi_x,a/2.0,b/2.0))
print("kappa x w_max = ",evaluateExpression(kappa_x,a/2.0,b/2.0))
print("Analytical w_max = ",evaluateExpression(w,a/2.0,b/2.0))


phi_x = Wmn*alpha*(sp.cos(alpha*x_s))*(sp.sin(beta*y_s))
kappa_x = -1.0*Wmn*alpha*alpha*(sp.sin(alpha*x_s))*(sp.sin(beta*y_s))
print("phi_x w_max = ",evaluateExpression(phi_x,a/2.0,b/2.0))
print("kappa x w_max = ",evaluateExpression(kappa_x,a/2.0,b/2.0))

# -----------------------------------------------------------------------------
#       READ KRATOS DISPLACEMENTS
# -----------------------------------------------------------------------------

kratos_w = []
kratos_w_x = []
with open("displacements.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      kratos_w_x.append(float(line[0]))
      kratos_w.append(float(line[1]))



# -----------------------------------------------------------------------------
#       COMPARE THEORETICAL STRAINS WITH KRATOS OUTPUT
# -----------------------------------------------------------------------------

kratos_x =[]
kratos_y =[]
kratos_kappa_y = []

with open("results.txt", "r") as kratos_file:
  for line in kratos_file:
      line = line.split('\t')
      kratos_x.append(float(line[0]))
      kratos_y.append(float(line[1]))
      kratos_kappa_y.append(float(line[5]))
      #print(float(line[5]))
          
print("\n\nKratos kappa y:\n",kratos_kappa_y)

analytical_w = []
analytical_kappa_y =[] 
analytical_kappa_x =[]
analytical_kappa_xy =[]
analytical_phi_x = []
analytical_phi_y = []
e_x =[]
e_y =[]
e_xy =[]
s_x = []
s_y = []
s_xy = []

for i in range(len(kratos_x)):
    analytical_w.append(evaluateExpression(w,kratos_x[i],kratos_y[i]))
    analytical_kappa_x.append(evaluateExpression(kappa_x,kratos_x[i],kratos_y[i]))
    analytical_kappa_y.append(evaluateExpression(kappa_y,kratos_x[i],kratos_y[i]))
    analytical_kappa_xy.append(evaluateExpression(kappa_xy,kratos_x[i],kratos_y[i]))
    analytical_phi_x.append(evaluateExpression(phi_x,kratos_x[i],kratos_y[i]))
    analytical_phi_y.append(evaluateExpression(phi_y,kratos_x[i],kratos_y[i]))
    e_x.append(h/2.0*analytical_kappa_x[i])
    e_y.append(h/2.0*analytical_kappa_y[i])
    e_xy.append(h/2.0*analytical_kappa_xy[i]/2.0)
    Mat_e = sp.zeros(3,1)
    Mat_e[0,0] = e_x[i]
    Mat_e[1,0] = e_y[i]
    Mat_e[2,0] = e_xy[i]
    Mat_s = C*Mat_e
    s_x.append(Mat_s[0,0])
    s_y.append(Mat_s[1,0])
    s_xy.append(Mat_s[2,0])
    
print("\n\nTheoretical kappa Y:\n",analytical_kappa_y)         



# -----------------------------------------------------------------------------
#       PLOT GRAPHS
# -----------------------------------------------------------------------------               
          
# DISPLACEMENTS
fig = plt.figure(1)
plt.plot(kratos_w_x,kratos_w, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_w,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
plt.xlim([0,200])
plt.xlabel('X Coordinate')
plt.ylabel('Transverse displacment')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")


#
## PHI X
#fig = plt.figure(2)
##plt.plot(kratos_w_x,kratos_w, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
##plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
#plt.plot(kratos_x,analytical_phi_x,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
##plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
#plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
##plt.legend()
##lg = plt.legend()
##lg.draw_frame(False)
##lg.loc(2)
##plt.xlim([0,90])
#plt.xlabel('X Coordinate')
#plt.ylabel('Transverse displacment')
#plt.grid()
#plt.tick_params(labelsize=labelfontsize)
##plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")
#
#

# PHI X
fig = plt.figure(3)
plt.plot(kratos_x,kratos_kappa_y, color = '#77B5FE', linewidth=3.0, label = 'ANDES-DKQ',antialiased=True)
#plt.plot(phi_dsg,n_theta_dsg, color = '#FF91A4', linewidth=3.0, label = 'DSG',antialiased=True)
plt.plot(kratos_x,analytical_kappa_y,color = 'grey', linestyle='None', markerfacecolor= 'None', markersize = 7.0, marker='o', label = 'Ref',antialiased=True)
#plt.plot(disp_ref,load_ref,color = 'grey', linewidth=2.0, linestyle='--', label = 'Ref',antialiased=True)
plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False,fontsize=labelfontsize+2)
#plt.legend()
#lg = plt.legend()
#lg.draw_frame(False)
#lg.loc(2)
#plt.xlim([0,90])
plt.xlabel('X Coordinate')
plt.ylabel('Transverse displacment')
plt.grid()
plt.tick_params(labelsize=labelfontsize)
#plt.savefig('Simply_support_dome_n_theta.pdf',bbox_inches="tight")



















plt.show()