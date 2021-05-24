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

# equal_bins = pd.DataFrame()
# equal_bins['cell_diamter'] = stats.mstats.mquantiles(df.cell_diameter,\
#                                     [0, 1./10, 2./10, 3./10, 4./10, 5./10,\
#                                      6./10, 7./10, 8./10, 9./10, 1])
# equal_bins['L'] = scipy.mean(df.L,\
#                                     [0, 1./10, 2./10, 3./10, 4./10, 5./10,\
#                                      6./10, 7./10, 8./10, 9./10, 1])    

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
#=============================================================================
#use curve fit to fit the data to the function defined above

#cell length
pars_d,covar_d =  curve_fit(powerlaw,df['cell_diameter'],df['L'],\
                                  absolute_sigma=False)
d_r = cof_det(df['L'], df['cell_diameter'], pars_d)

#calculate the standard deviation for the fit parameters
sigma_d = np.sqrt(np.diagonal(covar_d))

#=============================================================================
#set formatting parameters

fg = 8 #figure size
ft = 24 #font size for plots
ms = 13 #marker size
st = 'ticks' #set the style of ticks
a = 0.1 #alpha value

l = np.linspace(3, 13, 21) #lengths to plot 
v = np.linspace(20, 350, 12) #volumes to plot 
w = np.linspace(2, 11, 21) #width to plot 

#=============================================================================
#plot cable length vs cell length with double log plots

with sns.axes_style(st):
    f, ax = plt.subplots(figsize=(fg, fg)) 
    
    ax.loglog(df_hap['cell_diameter'],df_hap['L'], marker='o', \
                mfc='#f6511d', mew=1, mec=None, markersize=ms,\
                    linestyle='None', label=None, alpha=a)
    
    ax.loglog(df_dip['cell_diameter'],df_dip['L'], marker='o', \
                mfc='#00a6ed', mew=1, mec=None, markersize=ms,\
                    linestyle='None', label=None, alpha=a)
        
    ax.loglog(df_cdcu['cell_diameter'],df_cdcu['L'], marker='o', \
                mfc='#7fb800', mew=1, mec=None, markersize=ms, \
                    linestyle='None', label= None, alpha=a)
    
    ax.loglog(df_cdci['cell_diameter'],df_cdci['L'], marker='o', \
                mfc='#ffb400', mew=1, mec=None, markersize=ms,\
                    linestyle='None', label= None, alpha=a)        

    ax.loglog(l, powerlaw(l,*pars_d), 'k', lw=3, \
               label= r"Fit, Slope = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
                   .format(pars_d[1],d_r,sigma_d[1]))
        
    ax.loglog(l,powerlaw(l,1.01,1), 'k', ls="--", lw=2, \
               label='Isometric scaling, Slope = 1')
        
    ax.loglog(d_binned_data['cell_diameter'], d_binned_data['L'], marker='s', \
                mfc='k', mew=3, mec='k', markersize=10,\
                    linestyle='None', label= 'Binned data')

    plt.errorbar(d_binned_data['cell_diameter'], d_binned_data['L'],\
                 xerr = d_binned_err['cell_diameter'],\
                     yerr = d_binned_err['L'],linestyle='none', \
                         capsize=5, markersize="12", color='k')        
               
    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    plt.xlabel(u'Mother cell length (${\mu}m$)', fontsize=ft)
    ax.tick_params('both', length=10, which='both')
    for axis in [ax.xaxis, ax.yaxis]:
        formatter_min = LogFormatter(labelOnlyBase=True)
        axis.set_minor_formatter(formatter_min)                       
    plt.legend(loc='upper left', prop={'size': 16})
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([1e0, 3e1]) 
    # plt.xlim([2, 2e1])        
    plt.tight_layout()
    # plt.savefig('210221_binned_cable_length_powerlaw_scaling.svg')

#=============================================================================
#use curve fit to fit the data to the function defined above

#cell volume
pars_v,covar_v =  curve_fit(powerlaw,df['cell_volume'],df['L'],\
                                  absolute_sigma=False)
v_r = cof_det(df['L'], df['cell_volume'], pars_v)

#calculate the standard deviation for the fit parameters
sigma_v = np.sqrt(np.diagonal(covar_v))

#=============================================================================
#plot cable length vs cell volume with double log plots

with sns.axes_style(st):
    f, ax = plt.subplots(figsize=(fg, fg)) 
    
    plt.loglog(df_hap['cell_volume'],df_hap['L'], marker='o', \
               mfc='#f6511d', mew=1, mec=None, markersize=ms,\
                   linestyle='None', label=None, alpha=a)
    
    plt.loglog(df_dip['cell_volume'],df_dip['L'], marker='o', \
               mfc='#00a6ed', mew=1, mec=None, markersize=ms,\
                   linestyle='None', label=None, alpha=a)
        
    plt.loglog(df_cdcu['cell_volume'],df_cdcu['L'], marker='o', \
               mfc='#7fb800', mew=1, mec=None, linestyle='None',\
               markersize=ms, label= None, alpha=a)
    
    plt.loglog(df_cdci['cell_volume'],df_cdci['L'], marker='o', \
               mfc='#ffb400', mew=1, mec=None, linestyle='None',\
               markersize=ms, label= None, alpha=a)        

    plt.loglog(v,powerlaw(v,*pars_v), 'k', lw=3, \
               label= r"Fit, Slope = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
                   .format(pars_v[1],v_r,sigma_v[1]))
        
    plt.loglog(v,powerlaw(v,0.18,1), \
                'k', ls="--", lw=2, label='Isometric scaling, Slope = 1')
       
    ax.loglog(v_binned_data['cell_volume'], v_binned_data['L'], marker='s', \
                mfc='k', mew=3, mec='k', markersize=10,\
                    linestyle='None', label= 'Binned data')

    plt.errorbar(v_binned_data['cell_volume'], v_binned_data['L'],\
                 xerr = v_binned_err['cell_volume'],\
                     yerr = v_binned_err['L'],linestyle='none', \
                         capsize=5, markersize="12", color='k')         

    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    plt.xlabel(u'Mother cell volume (${\mu}m^3$)', fontsize=ft)
    ax.tick_params('both', length=10, which='both')
    for axis in [ax.xaxis, ax.yaxis]:
        formatter_min = LogFormatter(labelOnlyBase=True)
        axis.set_minor_formatter(formatter_min)                       
    plt.legend(loc='upper left', prop={'size': 16})
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([1e0, 3e1]) 
#    plt.xlim([1.5e0, 4e0])        
    plt.tight_layout()
    # plt.savefig('210221_binned_cable_volume_powerlaw_scaling.svg')
    
#=============================================================================
#use curve fit to fit the data to the function defined above

#cell width
pars_w,covar_w =  curve_fit(powerlaw,df['cell_diameter_2'],df['L'],\
                                  absolute_sigma=False)
w_r = cof_det(df['L'], df['cell_diameter_2'], pars_w)

#calculate the standard deviation for the fit parameters
sigma_w = np.sqrt(np.diagonal(covar_w))

#=============================================================================
#plot cable length vs cell width with double log plots

with sns.axes_style(st):
    f, ax = plt.subplots(figsize=(fg, fg)) 
    
    plt.loglog(df_hap['cell_diameter_2'],df_hap['L'], marker='o', \
               mfc='#f6511d', mew=1, mec=None, markersize=ms,\
                   linestyle='None', label=None, alpha=a)
    
    plt.loglog(df_dip['cell_diameter_2'],df_dip['L'], marker='o', \
               mfc='#00a6ed', mew=1, mec=None, markersize=ms,\
                   linestyle='None', label=None, alpha=a)
        
    plt.loglog(df_cdcu['cell_diameter_2'],df_cdcu['L'], marker='o', \
               mfc='#7fb800', mew=1, mec=None, linestyle='None',\
               markersize=ms, label= None, alpha=a)
    
    plt.loglog(df_cdci['cell_diameter_2'],df_cdci['L'], marker='o', \
               mfc='#ffb400', mew=1, mec=None, linestyle='None',\
               markersize=ms, label= None, alpha=a)        

    plt.loglog(w, powerlaw(w,*pars_w), 'k',lw=3, \
               label= r"Fit, Slope = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
                   .format(pars_w[1],w_r,sigma_w[1]))
        
    plt.loglog(w,powerlaw(w,1.45,1), 'k', ls="--", lw=2, \
               label='Isometric scaling, Slope = 1')
        
       
    ax.loglog(w_binned_data['cell_diameter_2'], w_binned_data['L'], marker='s', \
                mfc='k', mew=3, mec='k', markersize=10,\
                    linestyle='None', label= 'Binned data')

    plt.errorbar(w_binned_data['cell_diameter_2'], w_binned_data['L'],\
                 xerr = w_binned_err['cell_diameter_2'],\
                     yerr = w_binned_err['L'],linestyle='none', \
                         capsize=5, markersize="12", color='k')          

    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    plt.xlabel(u'Mother cell width (${\mu}m$)', fontsize=ft)
    ax.tick_params('both', length=10, which='both')
    for axis in [ax.xaxis, ax.yaxis]:
        formatter_min = LogFormatter(labelOnlyBase=True)
        axis.set_minor_formatter(formatter_min)                       
    plt.legend(loc='upper left', prop={'size': 16})
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([1e0, 3e1]) 
#    plt.xlim([1.5e0, 4e0])        
    plt.tight_layout()
    # plt.savefig('210221_binned_cable_width_powerlaw_scaling.svg')    










