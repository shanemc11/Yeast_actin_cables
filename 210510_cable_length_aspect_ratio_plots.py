# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 16:39:25 2019

@author: Shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy
import scipy.stats
import glob
import statsmodels.stats.api as sms
#import matplotlib for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import ScalarFormatter, LogFormatter
import seaborn as sns
import math
from scipy.spatial import distance

#import os to handle operating system
import os
#=============================================================================
#Goal: Compare measured cable architecture parameters from different strains/
#mutants. 
#=============================================================================
#indicate the data directory
datadir = "D:\\Goode_Lab\\Projects\\actin_cables\\data\\summary_data\\"

#initalize data frame to append all data 
df = pd.DataFrame()

#import data to dataframe
df = pd.read_csv(datadir + '210501_cdc28-13ts_t-8_t1_yBG12_yBG9_all_cable_analysis.csv')

#calculate the ratio of length versus diameter and volume for each cable
df['len_ratio'] = df['L'] / df['cell_diameter']

#calculate the mean for each replicate and strain for defining bins
df_all_mean = df.groupby(['strain', 'cell_number'],\
                                        sort=False).mean().reset_index()
    
df_all_mean['AR_round'] = df_all_mean['aspect_ratio'].round(2)
    
#define bins by rounding the aspect ratio to the nearest quarter value
df_1 = df_all_mean[(df_all_mean['AR_round'] < 1.12)].reset_index()  
df_1['Bin'] = 1

df_2 = df_all_mean[(df_all_mean['AR_round'] >=1.12) & (df_all_mean['AR_round']<=1.37)].reset_index()
df_2['Bin'] = 2  

df_3 = df_all_mean[(df_all_mean['AR_round'] >=1.38) & (df_all_mean['AR_round']<=1.62)].reset_index()
df_3['Bin'] = 3 

df_4 = df_all_mean[(df_all_mean['AR_round'] > 1.63 )].reset_index() 
df_4['Bin'] = 4 
  
frames = [df_1, df_2, df_3, df_4]

df_bins = pd.concat(frames)
    
#=============================================================================
     
#import curve_fit from scipy
from scipy.optimize import curve_fit

#write a function to fit using the power law with curve fit
def powerlaw(x,a,b):
    '''
    Parameters
    ----------
    x : The size dimension to scale intensity/length data with.
    a : Constant.
    b : Scaling exponent.

    Returns
    -------
    Use to fit data and extract scaling exponent for various cellular
    dimensions.

    '''
    y = a*(x**b)
    
    return y

#==============================================================================
#write a function to calculate the coefficient of determination for powerlaw
#fits

def cof_det(y, x, z):
    '''
    Parameters
    ----------
    y : dependent variable from powerlaw fits (ex: intensity, length).
    x : independent varibale from powerlaw fits (cells size dimensions).
    z : fitting parameter from powerlaw fits (scaling exponent and constant).

    Returns
    -------
    r_squared : coefficient of determination of powerlaw fit.

    '''
    
    res = y - powerlaw(x, *z)
    
    ss_res = np.sum(res**2)
    
    ss_tot = np.sum(( y - np.mean(y))**2)

                    
    r_squared = 1 - (ss_res / ss_tot)

    return (r_squared)

#==============================================================================
#write a function to calculate the powerlaw scaling parameters, the 
#coefficient of determination, and plot the results on a loglog plot

def scaling_plot(powerlaw, cof_det, x, y, c, txt, sv):
    '''
    

    Parameters
    ----------
    powerlaw : function
        Fits the data with the powerlaw to measure the scaling coeffcient and
        constant.
    cof_det : function
        Calculates the coefficient of determination.
    x : variable
        Dependent variable to use in the above functions.
    y : variable
        Independent variable to use in the above functions.
    c : string
        Color for markers.        
    txt : string
        Label for plot.
    sv : string
        Label for plot file during save.        
                

    Returns
    -------
    Results from fitting the data with the powerlaw and a loglog plot.

    '''
    pars, covar = curve_fit(powerlaw,x,y)
    
    r2 = cof_det(y, x, pars)
    
    #calculate the standard deviation for the fit parameters
    sigma = np.sqrt(np.diagonal(covar))
    
    with sns.axes_style('ticks'):
        f, ax = plt.subplots(figsize=(8, 8))
        plt.loglog(x,powerlaw(x,*pars),'k', lw=3, \
                   label= r"Fit, Slope = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
                   .format(pars[1], r2, sigma[1]))
        plt.loglog(x, y, mew=1, marker='o', markersize=13, \
                   linestyle='None', mfc=c, mec='k', label=txt)
        
        for axis in [ax.xaxis, ax.yaxis]:
            formatter_min = LogFormatter(labelOnlyBase=True)
            axis.set_minor_formatter(formatter_min)            
        ax.tick_params('both', length=10, which='both')
        plt.ylabel(u'Average cable length (${\mu}m$)', fontsize=24)
        plt.xlabel(u'Cell length (${\mu}m$)', fontsize=24)
        plt.legend(loc='upper left', prop={'size': 20})
        plt.rc('xtick', labelsize=24)
        plt.rc('ytick', labelsize=24)
        plt.ylim([1e0, 3e1]) 
        plt.xlim([3, 12])        
        plt.tight_layout()
        plt.savefig(sv) 

#==============================================================================
#make plots for each bin
#colors to use: '#441853', '#3C538C', '#1E928D', '#68BE63'
             
scaling_plot(powerlaw, cof_det,\
              df_1['cell_diameter'],df_1['L'],\
              '#441853', 'Bin 1',\
                  '210504_AR_bin1_scaling_plot.svg')
    
    
#==============================================================================
#import joypy for ridgeline plots
import joypy

# plot distributions of L vs diameter 
fig, axes = joypy.joyplot(df_bins, column=['aspect_ratio'], overlap=1.0, \
                          ylim='own', x_range=(0,3), fill=True, by="Bin",\
                          figsize=(10,12), legend=False, xlabels=True,\
                          ylabels=False, bw_method=1, colormap=cm.viridis,\
                          alpha=1.0, linewidth=2, linecolor='k')#, fade=True)
# plt.axvline(x=1.0, lw=2, color='k')    
plt.tick_params('both', length=10, width=3, which='both')
plt.rc("font", size=30)
plt.xlabel('Aspect ratio', fontsize=32, color='k', alpha=1)
# plt.savefig('210505_aspect_ratio_bin_distributions.svg')  

# plot distributions of L vs diameter 
fig, axes = joypy.joyplot(df_bins, column=['len_ratio'], overlap=1.0, by="Bin",\
                          ylim='own', x_range=(0,2), fill=True,\
                          figsize=(10,12), legend=False, xlabels=True,\
                          ylabels=False, bw_method=1,colormap=cm.viridis, \
                          alpha=1.0, linewidth=2, linecolor='k')#, fade=True)
plt.axvline(x=1.0, lw=2, color='k')    
plt.tick_params('both', length=10, width=3, which='both')
plt.rc("font", size=30)
plt.xlabel('Avg. Cable length / Cell length ratio', fontsize=32, color='k',\
           alpha=1)    
# plt.savefig('210505_cable_len_cell_len_bin_distributions.svg')  






