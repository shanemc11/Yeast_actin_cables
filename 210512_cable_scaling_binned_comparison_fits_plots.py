# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 16:39:25 2019

@author: Shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy
import scipy.stats as stats
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
#import files to analyze
datadir = "D:\\Goode_Lab\\Projects\\actin_cables\\data\\summary_data\\"

#initalize data frame to append all data 
df = pd.DataFrame()
#import data to dataframe
df = pd.read_csv(datadir + '201210_cdc28-13ts_t-8_t1_yBG12_yBG9_all_cable_analysis.csv')

#=============================================================================
#parse the data into the necessary strain types for plotting
#setup df with only yBG12 cells
df_hap = pd.DataFrame()

df_hap = df.loc[(df['strain']=='yBG12')].reset_index()

#setup df with only yBG9 cells
df_dip = pd.DataFrame()

df_dip = df.loc[(df['strain']=='yBG9')].reset_index()

#setup df with only uninduced cdc28 cells
df_cdcu = pd.DataFrame()

df_cdcu = df.loc[(df['strain']=='cdc28-13ts, t0')].reset_index()

#setup df with only induced cdc28 cells
df_cdci = pd.DataFrame()

df_cdci = df.loc[(df['strain']=='cdc28-13ts, t8')].reset_index()

#==============================================================================
# Bin the data frame by "cell diameter" with 10 bins...
d_bins = np.linspace(df.cell_diameter.min(), df.cell_diameter.max(), 11)
# Get the mean of parameters, binned by the values in cell diameter
d_binned_data = pd.DataFrame()
d_binned_data = df.groupby(pd.cut(df.cell_diameter, d_bins)).mean()

d_binned_err = pd.DataFrame()
d_binned_err = df.groupby(pd.cut(df.cell_diameter, d_bins)).std()

# Bin the data frame by "cell volume" with 10 bins...
v_bins = np.linspace(df.cell_volume.min(), df.cell_volume.max(), 11)
# Get the mean of parameters, binned by the values in cell diameter
v_binned_data = pd.DataFrame()
v_binned_data = df.groupby(pd.cut(df.cell_volume, v_bins)).mean()

v_binned_err = pd.DataFrame()
v_binned_err = df.groupby(pd.cut(df.cell_volume, v_bins)).std()

# Bin the data frame by "cell width" with 10 bins...
w_bins = np.linspace(df.cell_diameter_2.min(), df.cell_diameter_2.max(), 11)
# Get the mean of parameters, binned by the values in cell diameter
w_binned_data = pd.DataFrame()
w_binned_data = df.groupby(pd.cut(df.cell_diameter_2, w_bins)).mean()

w_binned_err = pd.DataFrame()
w_binned_err = df.groupby(pd.cut(df.cell_diameter_2, w_bins)).std()

#==============================================================================

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

def scaling_plot(powerlaw, cof_det, x, x_bin, x_bin_err, y, y_bin, y_bin_err,\
                 c, txt, sv):
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
    x_bin : variable
        Binned values for x.
    x_bin_err: variable
        Errors for x_bin values.        
    y : variable
        Independent variable to use in the above functions.
    y_bin : variable
        Binned values for y.
    y_bin_err: variable
        Errors for y_bin values.        
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

        plt.loglog(x, y, mew=1, marker='o', markersize=13, \
                    linestyle='None', mfc=c, mec='k', label=txt)

        plt.loglog(x, powerlaw(x,*pars), 'k',lw=3, \
               label= r"Fit, Slope = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
                   .format(pars[1],r2,sigma[1]))
       
        ax.loglog(x_bin, y_bin, marker='s', \
                mfc='k', mew=3, mec='k', markersize=10,\
                    linestyle='None', label= 'Binned data')

        plt.errorbar(x_bin, y_bin,\
                 xerr = x_bin_err,\
                     yerr = y_bin_err,linestyle='none', \
                         capsize=5, markersize="12", color='k')          

        ax.tick_params('both', length=10, which='both')
        for axis in [ax.xaxis, ax.yaxis]:
            formatter_min = LogFormatter(labelOnlyBase=True)
            axis.set_minor_formatter(formatter_min)                       
        plt.legend(loc='upper left', prop={'size': 25})
        plt.ylim([1e0, 3e1]) 
        
        plt.ylabel(u'Cable length (${\mu}m$)', fontsize=24)
        plt.xlabel(u'Mother cell width (${\mu}m$)', fontsize=24)
        plt.legend(loc='upper left')
        plt.rc('xtick', labelsize=24)
        plt.rc('ytick', labelsize=24)        
        plt.tight_layout()
        plt.savefig(sv) 

#==============================================================================
#make plots for bud/mom, total_cell_corr_int, cable_corr_int, patch_corr_int
#colors to use: '#f6511d', '#00a6ed', '#7fb800', '#ffb400', '#1CCAD8'
             
scaling_plot(powerlaw, cof_det, df['cell_diameter'], 
             d_binned_data['cell_diameter'], 
             d_binned_err['cell_diameter'], 
             df['L'],
             d_binned_data['L'], 
             d_binned_err['L'],\
             '#00a6ed', 'Cell length', '210512_test.svg')
    
#==============================================================================








