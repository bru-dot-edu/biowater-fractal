# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 21:24:59 2016

@author: Christopher L Brueck
"""

import numpy as np
import os
from matplotlib import pyplot as pl
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic

#Define the working directory
print(os.getcwd())
#os.chdir("Y:/Chris/Research/Fractal_Film_Model")
os.chdir("/Users/chrisbrueck/Dropbox/SOIL_635/Paper/Model")

#Create the particle distribution array
w = np.loadtxt("particle_size_data.txt",skiprows=1,delimiter=',')
w_x = w[:,0]
w_y = w[:,1]


def fit_func_x(x, a, b, c):
    return a*np.arctan(b*x)+c
    
xdata = np.linspace(min(w_x),max(w_x),100)
ydata = np.linspace(min(w_y),max(w_y),100)

popt, pcov = curve_fit(fit_func_x, w_x, w_y, p0=(0.6, 0.0169, -0.3518))

bin_means, bin_edges, binnumber = binned_statistic(xdata, fit_func_x(xdata, *popt),
         statistic='mean', bins=25)
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2
PSD1 = np.vstack((bin_centers, bin_means))
new_PSD1 = np.transpose(PSD1)
newest_PSD1 = np.flipud(new_PSD1)

pl.plot(xdata, fit_func_x(xdata, *popt), label="Fitted Curve")
pl.plot(w_x,w_y, 'ro', label="original data")
pl.plot(newest_PSD1[:,0],newest_PSD1[:,1], 'g^')
pl.legend(loc="lower right")
pl.show()




bin_means_normal = bin_means/sum(bin_means)
           
PSD = np.vstack((bin_centers, bin_means_normal))
new_PSD = np.transpose(PSD)
newest_PSD = np.flipud(new_PSD)


np.savetxt("binned_PSD_mass_fraction.txt",newest_PSD,delimiter = ',', header="diameter (micron), mass fraction", comments='')
