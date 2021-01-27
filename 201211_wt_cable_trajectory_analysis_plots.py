# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:29:57 2020

@author: Shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy
import scipy.stats
import operator
from operator import truediv
import glob
import statsmodels.stats.api as sms
#import matplotlib for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import seaborn as sns
import math
from math import sqrt
from scipy.spatial import distance

#import os to handle operating system
import os
#============================================================================= 
#Goal: Import appended datasets to generate summary plots.

#==============================================================================
#setup the data directory
datadir = "D:\\Goode_Lab\\Projects\\actin_cables\\data\\cable_trajectory_data\\"


#initalize data frame to append all data 
df_wt = pd.DataFrame()

#read in the summary data file
df_wt = pd.read_csv(datadir + \
                    "200901_all_yBG12-Abp140Envy_trajectories_cutoff.csv")    

#group data by experimental replicate and time to calculate mean values
df_expt_mean = df_wt.groupby(['lifetime', 'n'],\
                                       sort=False).mean().reset_index()
#=============================================================================    
# initiate plotting parameters
ft=22    
t_tot = np.linspace(0,40,41)
st = 'ticks' #set the style of ticks
    
#=============================================================================

#plot the change in extension rate over time

with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    
    #plot the theoretical constant extension rate
    plt.hlines(y=0.36, xmin=0, xmax=45, color='#fefe22', linestyle='--',\
                label='Constant extension', lw=3)  

    #plot the mean from each expt replicate
    ax = sns.scatterplot(x=df_expt_mean['lifetime'], y=df_expt_mean['vel'],
                      color = '#C200FB', label='Experiment', edgecolor='k',\
                       linewidth=1, alpha=1, s=80)
        
    #plot the mean of the replicates with the 95%CI    
    ax = sns.lineplot(x=df_expt_mean['lifetime'], y=df_expt_mean['vel'],
                      color='k', ci=95, label='Mean', lw=3)        
      
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft)
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))    
    ax.tick_params('both', length=5, which='both')    
    plt.legend([],[], frameon=False)    
    plt.ylim([0, 0.4])
    plt.xlim([0, 45])        
    plt.tight_layout()
    # plt.savefig('201217_wt_cable_ext_vs_lifetime.svg') 

#=============================================================================
    
#plot the change in extension rate over time for each replicate
with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    
    #plot the mean from each expt replicate        
    ax = sns.lineplot(x=df_expt_mean['lifetime'], y=df_expt_mean['vel'],
                       hue=df_expt_mean['n'], palette='Set1', lw=3)        
      
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft)
    plt.rc('xtick', labelsize=20)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))    
    plt.rc('ytick', labelsize=20)
    ax.tick_params('both', length=5, which='both')    
    plt.ylim([0, 0.4])
    plt.xlim([0, 45])        
    plt.tight_layout()
    # plt.savefig('201217_wt_cable_ext_vs_lifetime_all_replicate_means.svg')

#=============================================================================

#calculate the mean values as a function of time
df_wt_t_mean = pd.DataFrame()
df_wt_t_mean = df_wt.groupby(['lifetime']).mean().reset_index() 

#use linear regression to fit the first ~10seconds of data to estimate
#the theoretical constant extension rate
from scipy.stats import linregress
slope_r0, intercept_r0, r_value_r0, p_value_r0, std_err_r0 = \
scipy.stats.linregress(df_wt_t_mean['lifetime'][:5], \
                       df_wt_t_mean['neck_dist'][:5])
    
    
#plot the change in cable length over time
    
with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    
    #plot the change in length as a function of time for all cables tracked
    sns.scatterplot(x=df_expt_mean['lifetime'], y=df_expt_mean['neck_dist'], \
              color='#C200FB',alpha=1, linewidth=1, edgecolor='k',\
                  label = 'Experiment', s=80)
        
    #plot the mean change in length as a function of time 
    #for all expt replicates    
    sns.lineplot(x=df_expt_mean['lifetime'], y=df_expt_mean['neck_dist'],\
                 ci=95, color = 'k', label = 'Mean', lw=3) 
    
    #plot the theoretical linear increase in cable length using the parameters
    #from linear regression    
    plt.plot(t_tot, (slope_r0*t_tot), linestyle='--', lw=3,\
             label = 'Boundary sensing model', color='#fefe22')        
              
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    ax.tick_params('both', length=5, which='both')
    plt.legend([],[], frameon=False)
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 10])
    plt.xlim([0, 45])        
    plt.tight_layout()
    # plt.savefig('201217_wt_len_vs_lifetime.svg') 

#=============================================================================
    
#plot the change in cable length over time for each replicate
with sns.axes_style(st):
    plt.figure(figsize=(5,5))
            
    ax = sns.lineplot(x=df_expt_mean['lifetime'], y=df_expt_mean['neck_dist'],
                       hue=df_expt_mean['n'], palette='Set1', lw=3)        
      
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    plt.legend(loc='upper left')
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    ax.tick_params('both', length=5, which='both')
    plt.ylim([0, 10])
    plt.xlim([0, 45])        
    plt.tight_layout()
    # plt.savefig('201217_wt_len_vs_lifetime_all_replicate_means.svg')    
    

#=============================================================================
    

