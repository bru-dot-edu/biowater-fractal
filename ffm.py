# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 10:12:05 2015

@author: Christopher L Brueck
"""

import numpy as np
import os
from matplotlib import pyplot as pl
import matplotlib as mpl
import math as m
from scipy.optimize import curve_fit


mpl.rcParams['font.family']='Helvetica'
mpl.rcParams['font.size']=10

#pl.style.use('ggplot')
#pl.style.use('bmh')
#pl.style.use('dark_background')
#pl.style.use('fivethirtyeight')
#pl.style.use('grayscale')

#Define the working directory
print(os.getcwd())
os.chdir("/Users/chrisbrueck/Dropbox/SOIL_635/Paper/Post_Review_Model/For_Github_or_Supplementary_Material")

#Create the particle distribution array
part_array = np.loadtxt("binned_PSD_mass_fraction.txt",skiprows=1,delimiter=',')
moisture_data = np.loadtxt("quincy_moisture_data.txt",skiprows=1,delimiter=',')


#==============================================================================
# Constants
#==============================================================================

imax = 5                            #max number of iterations for each box + 1
sigma = 72.7                        #dyne/cm surface tension
rho = 1                             #g/cm3 density
g = 980                             #cm/s2 gravity
unit = 10**4                        #unit conversion 

theta_r = 0.015                     #van genuchten parameter
theta_s = 0.50136125949683663       #van genuchten parameter
alpha = 0.078022347398673184        #van genuchten parameter
n_van = 1.4023848222981605          #van genuchten parameter

#==============================================================================
# create iteration array
#==============================================================================

num_part = len(part_array)
n = np.arange(2,imax)
num_iter = np.tile(n,(num_part,1))
num_iter2 = num_iter.reshape(-1)

#==============================================================================
# Create the sidelength array
#==============================================================================

S_list=[]
for i in range(0,num_part):
    s_tmp = part_array[i,0]/3.0
    for step in range(2,imax):
        s_tmp = s_tmp/3.0
        S_list.append(s_tmp)
        print(s_tmp)
        
S_list=np.asarray(S_list)

#==============================================================================
# Equations for numbers of hexagons and triangles per iteration
#==============================================================================

a = []
b = []

for step in range(0,num_part):
    a_tmp = 1
    b_tmp = 2
    for i in range(0,imax-2):
        a.append(a_tmp)
        b.append(b_tmp)
        a_tmp = 3.0*a_tmp + b_tmp
        b_tmp = 2*a[i] + 2*b_tmp

a = np.asarray(a)
b = np.asarray(b)

#########################################################################
# Calculate Total Solid Area for Each Unit Pore 
#########################################################################

Solid_Area = []

for i in range(0,len(part_array)):
    Solid_Area_tmp = m.sqrt(3)/4*((part_array[i,0]/3)**2)*4
    Solid_Area2 = []
    for j in range(2,imax):
        Solid_Area2_tmp = m.sqrt(3)/4*((part_array[i,0]/(3**j))**2)*(4**j)
        Solid_Area2.append(Solid_Area2_tmp)
    Solid_Area2_sum = np.sum(Solid_Area2)
    tot = Solid_Area_tmp + Solid_Area2_sum
    Solid_Area.append(tot)
Solid_Area = np.asarray(Solid_Area)

#########################################################################
# Expand Solid_Area to play nice with other matrices
#########################################################################

Solid_Area_ext = []

for i in range(0,len(part_array)):
    SolArea_tmp = Solid_Area[i]
    for j in range(2,imax):
        Solid_Area_ext.append(SolArea_tmp)

#==============================================================================
# Fractal Model
#==============================================================================

d = m.sqrt(3)*S_list                                #equation 5 in paper

rc = d/2.0*3.0           #equivalent to equation 18 but uses dn instead of dn-1

Ahex = 3.0/2.0*m.sqrt(3)*(S_list**2.0)              #equation 9 in paper

s = (m.sqrt(3.0)+2.0)/2.0*S_list                    #equation 5 in paper

Ab = np.sqrt(s*(s-S_list)*(s-d)*(s-S_list))         #equation 6 in papaer

Aa = Ahex - Ab                                      #equation 8 in paper

At = 3.0/8.0*(S_list**2.0)                          #equation 11 in paper

Af = a*Aa+b*Ab+At                                   #equation 12 in paper

F_excess = 8*Af                                     #equation 13 in paper

F = []
init = 0
init1 = 1
init2 = 2
for i in range(0,int(len(F_excess)/3)):
    n_1 = m.sqrt(3)/4*((part_array[i,0]/(3**1))**2)*(4**1)
    n_2 = m.sqrt(3)/4*((part_array[i,0]/(3**2))**2)*(4**2)
    n_3 = m.sqrt(3)/4*((part_array[i,0]/(3**3))**2)*(4**3)
    tmp1 = F_excess[init] - (Solid_Area_ext[init] - n_1 - n_2)
    tmp2 = F_excess[init1] - (Solid_Area_ext[init1] - n_1 - n_2 - n_3)
    tmp3 = F_excess[init2]
    F.append(tmp1)
    F.append(tmp2)
    F.append(tmp3)
    init = init + 3
    init1 = init1 + 3 
    init2 = init2 + 3
F = np.asarray(F)    

H = sigma/rho/g/rc*unit                             #equation 19 in paper

#########################################################################
# Calculate area of empty pore units
#########################################################################

BoxArea=[]
for i in range(0,num_part):
    for j in range(0,imax-2):
        BoxArea.append(part_array[i,0])
BoxArea = np.asarray(BoxArea)**2

#########################################################################
# Calculate void area 
#########################################################################
        
VoidArea = BoxArea - Solid_Area_ext

#########################################################################
# Weighting the film area with volume fraction
#########################################################################

F_adjusted=[]
init = 0
for i in range(0,num_part):
    F_adjusted[init:imax-2+init] = F[init:imax-2+init]*part_array[i,1]
    init = init + imax - 2
F_adjusted=np.asarray(F_adjusted)

Sat = F_adjusted/VoidArea

#==============================================================================
# Define the fitting function
#==============================================================================
    
def fit_func_x(x, a, b, c):
    return a*(x**b)+c

#==============================================================================
# Plotting the H-S data from each fractal and fitting to the function
# Saves the parameters for each fitted function into a matrix
#==============================================================================

init = 0 
params_save = []    
fig = pl.figure()
ax = fig.add_subplot(1,1,1)
H_min_max = []

for i in range(0,num_part):
    x_data = H[init:imax-2+init]
    y_data = Sat[init:imax-2+init]
    fit_range = np.linspace(H[init],H[imax-3+init], 50)
    params, pcov = curve_fit(fit_func_x, x_data, y_data, p0=(2, -.67, 0.0))
    pl.plot(x_data, y_data, 'o', fit_range, fit_func_x(fit_range, params[0], params[1], params[2]))
    params_save.append(params)
    h_minmax = [min(x_data),max(x_data)]
    H_min_max.append(h_minmax)
    init = init + imax - 2
    pl.xlabel("Capillary Pressure, H (cm)")
    pl.ylabel("Saturation, S")
    
#pl.savefig("X:/Chris/Research/Fractal_Film_Model/individual_S-P.pdf", bbox_inches='tight')
params_save = np.asarray(params_save)

#==============================================================================
# Plotting the fitted functions from each fractal box over a wider range of H values
# Saves an array of the sum of all saturations for each fitted function (Sat_sum)
# Saves an array of each of the fitted functions (Sat_f)
#==============================================================================

fig2 = pl.figure()
ax = fig2.add_subplot(1,1,1)
Sat_sum = np.empty(100)
Sat_f=[]
init = 0

for i in range(0,num_part):   
    y = np.linspace(min(H),max(H),100)
    y_fit = np.linspace(H[init],H[imax-3+init],100)
    Sat_final_tmp = fit_func_x(y, params_save[i,0], params_save[i,1], params_save[i,2])
    print(Sat_final_tmp)
    Sat_sum = Sat_sum + Sat_final_tmp
    Sat_f.append(Sat_final_tmp)
    #pl.plot(Sat_sum,y_fit)
    pl.plot(Sat_final_tmp,y)
    #pl.plot(Sat_sum,y_fit)
    pl.ylabel("Capillary Pressure, H (cm)")
    pl.xlabel("Saturation, S") 
    init = init + imax - 2
    
#pl.savefig("X:/Chris/Research/Fractal_Film_Model/sum_P-S.pdf", bbox_inches='tight') 

Sat_f = np.asarray(Sat_f)
Sat_f = np.transpose(Sat_f)

#Add each fitted equation to one another without extending limits

S_fit = np.zeros((round(max(H))-round(min(H))+1))
o = round(min(H))
init = 0
for i in range (0,num_part):
    k = round(H[init])
    l = round(H[init+imax-3])
    y_round=np.linspace(k,l,l-k+1)
    S_f_tmp = fit_func_x(y_round, params_save[i,0], params_save[i,1], params_save[i,2])
    S_fit[y_round[0]-o:len(y_round)+y_round[0]-o] = S_fit[y_round[0]-o:len(y_round)+y_round[0]-o] + S_f_tmp
    init = init + imax - 2

y = np.linspace(round(min(H)),round(max(H)),len(S_fit))

#==============================================================================
# Van Genuchten model
#==============================================================================

y_fit_van = np.linspace(105000,0,100000)
#commented because it will crash a computer unless it has ~40 GB of ram
S_van = 1/((1+(alpha*y_fit_van)**n_van)**(1-1/n_van))
S_van_limit = 1/((1+(alpha*y)**n_van)**(1-1/n_van))
quincy_sat = (moisture_data[:,0]-theta_r)/(theta_s-theta_r)
#quincy_sat = moisture_data[:,0]/theta_s
quincy_pressure = moisture_data[:,1]

S_van_subtract = S_van_limit - S_fit
S_ratio = S_fit/S_van_limit*100
S_sub_ratio = S_fit/S_van_subtract*100
S_cap_ratio = S_van_subtract/S_van_limit*100

#==============================================================================
# Create plots/figures
#==============================================================================

fig3 = pl.figure(figsize=(6,4))
ax = fig3.add_subplot(1,1,1)
pl.plot(S_van,y_fit_van,'-.',label='Van Genuchten model')
pl.plot(S_fit,y,label='Thick film model', color="green")
pl.plot(quincy_sat, quincy_pressure, 'ro', alpha=0.7, label='Moisture retention data')
pl.plot(S_van_subtract,y, 'k--', label='Capillary contribution')
#pl.hlines(y=max(y),xmin=0,xmax=min(S_van_limit),color='k',linestyle='dashed')
#pl.hlines(y=min(y),xmin=0,xmax=max(S_van_limit),color='k',linestyle='dashed')
pl.xlabel('Saturation [-]')
pl.ylabel('Capillary Pressure [cm]')
pl.legend(loc="upper right")
ax.set_yscale('log')
ax.set_ylim(0.08,1000000)
ax.set_xlim(0,1.05)
fig3.savefig('fig3.png', bbox_inches='tight', dpi=1200)
fig3.savefig('fig3.pdf', bbox_inches='tight')

fig4 = pl.figure(figsize=(6,4))
ax1 = fig4.add_subplot(1,1,1)
pl.plot(y_fit_van,S_van,'-.',label='Van Genuchten model')
pl.plot(y,S_fit,label='Thick film model', color="green")
pl.plot(quincy_pressure,quincy_sat,'ro', alpha=0.7, label='Moisture retention data')
pl.plot(y,S_van_subtract,'k--', label='Capillary contribution')
#pl.hlines(y=max(y),xmin=0,xmax=min(S_van_limit),color='k',linestyle='dashed')
#pl.hlines(y=min(y),xmin=0,xmax=max(S_van_limit),color='k',linestyle='dashed')
pl.ylabel('Saturation [-]')
pl.xlabel('Capillary Pressure [cm]')
pl.legend(loc="upper right")
ax1.set_xscale('log')
ax1.set_xlim(0.08,1000000)
ax1.set_ylim(0,1.05)
fig4.savefig('fig4.png', bbox_inches='tight', dpi=1200)
fig4.savefig('fig4.pdf', bbox_inches='tight')

fig5 = pl.figure(figsize=(6,4))
ax5 = fig5.add_subplot(1,1,1)
pl.plot(y,S_cap_ratio,':',  label='Pore-filled capillary water pool')
pl.plot(y,S_ratio, 'm',label='Thick film water pool')
pl.ylabel('% Contribution')
pl.xlabel('Capillary Pressure [cm]')
pl.legend(loc='center right')
ax5.set_ylim(0,100)
pl.xlim((10,1000))
#ax4.set_xticks([10,20,30, 100])
#pl.xscale('log')
fig5.savefig('fig5.png', bbox_inches='tight', dpi=1200)
fig5.savefig('fig5.pdf', bbox_inches='tight')

fig6 = pl.figure(figsize=(6,4))
ax6 = fig6.add_subplot(1,1,1)
pl.plot(y,S_van_subtract, 'b')
pl.plot(y_fit_van[99989:99999], S_van[99989:99999], 'b')
pl.plot(y_fit_van[0:99520], S_van[0:99520], 'b')
pl.plot(y,S_fit, 'g')
pl.plot(y,S_van_limit, 'k--')
ax6.set_xscale('log')
pl.ylim(0,1)
pl.ylabel('Saturation [-]')
pl.xlabel('Capillary Pressure [cm]')
