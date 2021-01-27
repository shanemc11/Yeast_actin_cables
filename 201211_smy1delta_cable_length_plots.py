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
import matplotlib.ticker as ticker
import seaborn as sns
import math
from scipy.spatial import distance

#import os to handle operating system
import os
#=============================================================================
#Goal: Compare measured cable lengths between wt and mutant cells.
 
#=============================================================================
#import files to analyze
datadir = "D:\\Goode_Lab\\Projects\\actin_cables\\data\\summary_data\\"

#initalize data frame to append all data 
df = pd.DataFrame()
#import data to dataframe
df = pd.read_csv(datadir + '200730_Abp140Smy1delta_alldata_summary_edit.csv')  

df_means = df.groupby(['strain', 'n'], sort=False).mean().reset_index()

df_means_sort = df.groupby(['n', 'strain']).mean().reset_index()

#============================================================================= 
#initialize plotting parameters

o = ['Abp140Envy', 'Abp140EnvySmy1d']
cmap = ["#C200FB", "#F28D35"] 
ft = 22 #font size for x axis
st = 'ticks' #set the style of ticks

#============================================================================= 
#plot the cable length in wt and smy1delta cells
    
with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    sns.set_palette(cmap)
    
    sns.swarmplot(x='strain', y='cell_volume', data = df, linewidth=0.5,size=10,\
                  alpha=1, edgecolor='k', zorder=0, dodge=True)      
        
    ax = sns.stripplot(x='strain', y='L', data = df_means_sort[:2], size=15,\
                       color='grey', edgecolor='k', marker="s",\
                       linewidth=1, dodge=True, \
                       order = o)
        
    ax = sns.stripplot(x='strain', y='L', data = df_means_sort[2:4], size=15,\
                       color='grey', edgecolor='k', marker="o",\
                       linewidth=1, dodge=True,\
                       order = o)
        
    ax = sns.stripplot(x='strain', y='L',\
                       data = df_means_sort[4:], size=15,\
                       color='grey', edgecolor='k', marker="^",\
                       linewidth=1, dodge=True,\
                       order = o)        
        
    ax = sns.pointplot(x='strain', y='L', data = df_means,\
                       capsize = 0.8, join=False, color='k')
        
    plt.ylabel(u'Cable length(${\mu}m$)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2))            
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=ft)
    plt.xlabel(None)
    ax.set(xticks=[])
    ax.tick_params('both', length=5, which='both')    
    plt.tight_layout()
    plt.ylim([0, 80])
    plt.legend([],[], frameon=False)
    # plt.savefig('201217_smy1d_wt_cells_cables_L.svg') 
    
    
#============================================================================
#do statistical tests on the data

#organize cable data by stain and expt number in a new dataframe
df_cable_mean = df.groupby(['strain', 'n'], sort=True).mean().reset_index()

df_cable_stats = df_cable_mean[['D', 'L', 'L/D']]

#organize cell data by stain and expt number in a new dataframe
df_cell_mean = df.groupby(['cell_number', 'strain', 'n'],\
                                       sort=True).mean().reset_index()
df_cell_stats = df_cell_mean[['cell_diameter', 'cell_volume']]


#use ttest to determine statistical signficance
stat_cable, pval_cable = scipy.stats.ttest_ind(df_cable_stats[:3],\
                                                         df_cable_stats[3:])
    
stat_cell, pval_cell = scipy.stats.ttest_ind(df_cell_stats[:3],\
                                                         df_cell_stats[3:])
#============================================================================
    
    












