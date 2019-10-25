# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 15:53:24 2016
@author: Christopher L Brueck
"""

import numpy as np
from scipy.optimize import curve_fit
import os
from matplotlib import pyplot as pl

os.chdir("/Users/chrisbrueck/Dropbox/SOIL_635/Paper/Post_Review_Model")

moisture_data = np.loadtxt("quincy_moisture_data.txt",skiprows=1,delimiter=',')

quincy_vol_moist_cont = moisture_data[:,0]
quincy_pressure_cm = moisture_data[:,1]
theta_r = 0.015

def func(x, n, alpha, theta_s):
    return 1/((1+(alpha*x)**n)**(1-1/n))*(theta_s-theta_r)+theta_r

popt, pcov = curve_fit(func, quincy_pressure_cm, quincy_vol_moist_cont, p0=(1.5,0.6,0.5))

n = popt[0]
alpha = popt[1]
theta_s = popt[2]

#n = 1.5
#alpha = 0.6
#theta_r = 0.02
#theta_s = 0.5
xs = np.linspace(0,4000,10000)
y = 1/((1+(alpha*xs)**n)**(1-1/n))*(theta_s-theta_r)+theta_r

fig = pl.figure()
axs = fig.add_subplot(1,1,1)
axs.set_yscale('log')
axs.set_ylim(0.1,10000)
axs.set_xlim(0.03,0.5)
pl.plot(y,xs)
pl.plot(moisture_data[:,0],moisture_data[:,1], 'ko')
