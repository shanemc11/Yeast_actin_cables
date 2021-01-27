# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 09:08:22 2020

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
import matplotlib.ticker as ticker
import seaborn as sns
import math
from scipy.spatial import distance

#import os to handle operating system
import os
#============================================================================= 
#Goal: 
#==============================================================================
#setup the data directory
datadir = "D:\\Goode_Lab\\Projects\\actin_cables\\data\\cable_segmented_integration\\"

#initalize data frame to append all data 
df = pd.DataFrame()

#import data to dataframe
df = pd.read_excel(datadir + '201218_cdc28-13ts_only_timecourse_d_v.xlsx')

df_mean = df.groupby(['strain']).mean().reset_index()

#=============================================================================
#set color palette to use for plots
ft = 18 #font size for x axis
ms = 60 #marker size

t = np.linspace(0, 8, 9) #lengths to plot model

st = 'ticks'
  
from scipy.stats import linregress

slope_d, intercept_d, r_value_d, p_value_d,\
    std_err_d = \
scipy.stats.linregress(df['time'],\
                       df['d1'])

slope_v, intercept_v, r_value_v, p_value_v,\
    std_err_v = \
scipy.stats.linregress(df['time'],\
                       df['volume'])
    
#=============================================================================

with sns.axes_style(st):
    plt.figure(figsize=(5,5))
            
    ax = sns.scatterplot(x=df['time'], y=df['d1'], color='#8d8d8d',\
              s=ms, linewidth=0.5, edgecolor='k')
        
    # ax = sns.lineplot(x=df['time'], y=df['d1'], \
    #           color='k', ci=95, lw=3)
    
    plt.plot(t,(intercept_d + slope_d*t),\
              'k--', lw=3,\
              label=r"Slope = {0:.2f}+/-{2:.2f}, R$^2$ = {1:.2f}".\
                  format(slope_d, r_value_d**2, 1.96*std_err_d)) 
 
    plt.xlabel('Time (hours)', fontsize=ft)
    plt.ylabel(u'Mother cell length (${\mu}m$)', fontsize=ft)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))        
    plt.legend(loc='upper left')    
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    # plt.ylim([0, 15])
    # plt.xlim([0, 10])        
    plt.tight_layout()
    plt.savefig('201217_cdc28-13ts_timecourse_cell_length.svg')

    

with sns.axes_style(st):
    plt.figure(figsize=(5,5))
            
    ax = sns.scatterplot(x=df['time'], y=df['volume'], color='#8d8d8d',\
              s= ms, linewidth=0.5, edgecolor='k')
        
    # ax = sns.lineplot(x=df['time'], y=df['volume'], \
    #           color='k', ci=95, lw=3)
        
    plt.plot(t,(intercept_v + slope_v*t),\
              'k--', lw=3,\
              label=r"Slope = {0:.1f}+/-{2:.1f}, R$^2$ = {1:.2f}".\
                  format(slope_v, r_value_v**2, 1.96*std_err_v)) 
 
    plt.xlabel('Time (hours)', fontsize=ft)
    plt.ylabel(u'Mother cell volume (${\mu}m^3$)', fontsize=ft)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))        
    plt.legend(loc='upper left')    
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    # plt.ylim([0, 15])
    # plt.xlim([0, 10])        
    plt.tight_layout()
    plt.savefig('201217_cdc28-13ts_timecourse_cell_volume.svg')
    
    


