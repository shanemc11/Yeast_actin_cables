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
#Goal: Import appended datasets to generate summary plots and conduct fits, 
# etc.
#==============================================================================
#setup the data directory
datadir = "D:\\Goode_Lab\\Projects\\actin_cables\\data\\cable_trajectory_data\\"

#initialize empty dataframes for the data
df_wt = pd.DataFrame()

df_smy1 = pd.DataFrame()

all_df = pd.DataFrame()

#import data to dataframes
df_wt = pd.read_csv(datadir + \
                    "200805_yBG12-Abp140Envy_cable_trajectory_analysis_cutoff.csv")

df_smy1 = pd.read_csv(datadir + \
                    "200805_yBG12-Abp140EnvySmy1delta_cable_trajectory_analysis_cutoff.csv")

all_df = pd.read_csv(datadir + \
                    "200805_smy1delta_wt_all_cable_extension_rate_analysis.csv")  
    
#=============================================================================    
#calculate means for each cell
df_wt_cell_mean = pd.DataFrame()
df_wt_cell_mean = df_wt.groupby(['cell']).mean() 

df_smy1_cell_mean = pd.DataFrame()
df_smy1_cell_mean = df_smy1.groupby(['cell']).mean()

#calculate means and std deviation for each time point
df_wt_t_mean = pd.DataFrame()
df_wt_t_mean = df_wt.groupby(['lifetime']).mean().reset_index() 

df_wt_t_std = pd.DataFrame()
df_wt_t_std = df_wt.groupby(['lifetime']).std().reset_index() 

 
df_smy1_t_mean = pd.DataFrame()
df_smy1_t_mean = df_smy1.groupby(['lifetime']).mean().reset_index()  
 
df_smy1_t_std = pd.DataFrame()
df_smy1_t_std = df_smy1.groupby(['lifetime']).std().reset_index()   

#calculate means for each timepoint for each replicate
df_wt_expt_mean = df_wt.groupby(['lifetime', 'n'],\
                                       sort=False).mean().reset_index()
    
df_smy1d_expt_mean = df_smy1.groupby(['lifetime', 'n'],\
                                       sort=False).mean().reset_index()
    
#=============================================================================
#set plotting parameters
cmap = ["#C200FB", "#F28D35"] 
ft = 22 #font size for x axis
st = 'ticks' #set the style of ticks

#=============================================================================
#plot the change in rate as a function of length 
#use the time averaged lengths and rates for this or the plots don't plot well
#due to the different values for each length

#calculate the errors (95%CI) for each strain   
wt_lower_ci = df_wt_t_mean - (1.96 * (df_wt_t_std / np.sqrt(65)))
wt_upper_ci = df_wt_t_mean + (1.96 * (df_wt_t_std / np.sqrt(65)))

smy1_lower_ci = df_smy1_t_mean - (1.96 * (df_smy1_t_std / np.sqrt(39)))
smy1_upper_ci = df_smy1_t_mean + (1.96 * (df_smy1_t_std / np.sqrt(39)))
   
#plot the change in rate as a function of length 
with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    sns.set_palette(cmap)   
        
    ax = sns.lineplot(x=df_wt_t_mean['neck_dist'], y=df_wt_t_mean['vel'], \
              ci=95, label='Abp140Envy', color='#C200FB', lw=3)
    
    ax = sns.lineplot(x=df_smy1_t_mean['neck_dist'], y=df_smy1_t_mean['vel'],\
                      color='#F28D35', ci=95, \
                          label = u'Abp140Envy;Smy1$\Delta$',\
                          lw=3)

    plt.fill_between(df_wt_t_mean['neck_dist'], wt_lower_ci['vel'],\
                     wt_upper_ci['vel'],\
                     color='#C200FB', alpha=0.3) 
        
    plt.fill_between(df_smy1_t_mean['neck_dist'], smy1_lower_ci['vel'],\
                     smy1_upper_ci['vel'],\
                     color='#F28D35', alpha=0.3)        
        
 
    plt.xlabel(u'Cable length (${\mu}m$)', fontsize=ft)
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))        
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))        
    ax.tick_params('both', length=5, which='both')    
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 0.5])
    plt.xlim([0, 6])        
    plt.tight_layout()
    # plt.savefig('201217_wt_v_smy1d_extrate_v_len.svg') 
    
#=============================================================================
#plot the extension rates using the mean of each wildtype replicate 

with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    #plot the mean and 95%CI for the replicates
    ax = sns.lineplot(x=df_wt_expt_mean['lifetime'],\
                      y=df_wt_expt_mean['vel'],
                      color='#C200FB', ci=95, label='Mean', lw=3) 
    #plot the mean for each replicates
    ax = sns.scatterplot(x=df_wt_expt_mean['lifetime'], \
                         y=df_wt_expt_mean['vel'],
                         color = 'grey', label='Experiment', edgecolor='k',\
                         linewidth=1, alpha=1, s=80)
      
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.tick_params('both', length=5, which='both')         
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 0.5])
    plt.xlim([0, 40])        
    plt.tight_layout()
    # plt.savefig('201217_wt_cont_cable_ext_vs_lifetime.svg')   
    
#plot the change in length as a function of time  

with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    #plot the mean for each replicates    
    sns.scatterplot(x=df_wt_expt_mean['lifetime'], \
                    y=df_wt_expt_mean['neck_dist'], \
                        color='#C200FB',alpha=1, linewidth=1, edgecolor='k',\
                            label = 'Experiment', s=80)
        
    #plot the mean and 95%CI for the replicates   
    sns.lineplot(x=df_wt_expt_mean['lifetime'],\
                 y=df_wt_expt_mean['neck_dist'],\
                 ci=95, color = 'k', label = 'Mean', lw=3) 
                      
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    ax.tick_params('both', length=5, which='both')        
    plt.legend(loc='upper left')
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 10])
    plt.xlim([0, 40])        
    plt.tight_layout()
    # plt.savefig('201217_wt_cont_len_vs_lifetime.svg')


#=============================================================================    
#plot the extension rates using the mean of each smy1delta replicate 
 
with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    #plot the mean and 95%CI for the replicates
    ax = sns.lineplot(x=df_smy1d_expt_mean['lifetime'],\
                      y=df_smy1d_expt_mean['vel'],
                      color='#F28D35', ci=95, label='Mean', lw=3) 
    #plot the mean for each replicates    
    ax = sns.scatterplot(x=df_smy1d_expt_mean['lifetime'], \
                         y=df_smy1d_expt_mean['vel'],\
                         color = 'grey', label='Experiment', edgecolor='k',\
                         linewidth=1, alpha=1, s=80) 
        
       
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.tick_params('both', length=5, which='both')         
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 0.5])
    plt.xlim([0, 40])        
    plt.tight_layout()
    # plt.savefig('201217_smy1delta_cable_ext_vs_lifetime.svg')   
    
     
with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    #plot the mean for each replicates    
    sns.scatterplot(x=df_smy1d_expt_mean['lifetime'],\
                    y=df_smy1d_expt_mean['neck_dist'], \
                        color='#F28D35',alpha=1, linewidth=1, edgecolor='k',\
                            label = 'Experiment', s=80)
        
     #plot the mean and 95%CI for the replicates       
    sns.lineplot(x=df_smy1d_expt_mean['lifetime'], \
                 y=df_smy1d_expt_mean['neck_dist'],\
                 ci=95, color = 'k', label = 'Mean', lw=3) 
        
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    ax.tick_params('both', length=5, which='both')    
    plt.legend(loc='upper left')
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 10])
    plt.xlim([0, 40])        
    plt.tight_layout()
    # plt.savefig('201217_smy1d_len_vs_lifetime.svg')     
    
#=============================================================================
#setup df with the initial extension rate
df2 = pd.DataFrame()
df2 = all_df.loc[(all_df['frame']== 0)]

#calculate the mean for each replicate and each strain
df3 = df2.groupby(['strain', 'n']).mean().reset_index()  
df3 = df3.sort_values('strain', ascending=False).reset_index()

df_means_sort = df2.groupby(['n', 'strain']).mean().reset_index()
df_means_sort = df_means_sort.sort_values(['n', 'strain'], \
                                          ascending=[True, False]).reset_index()

#plot the initial extension rate
with sns.axes_style(st):
    plt.figure(figsize=(5,6))
    sns.set_palette(cmap)
    
    sns.swarmplot(x='strain', y='vel', data = df2, \
                  linewidth=0.5,\
                  alpha=1, edgecolor='k', size=10, dodge=True)  
        
    ax = sns.stripplot(x='strain', y='vel', data = df_means_sort[:2], size=15,\
                       color='grey', edgecolor='k', marker="s",\
                       linewidth=1, dodge=True, \
                       order = ['wt', 'smy1delta'])
        
    ax = sns.stripplot(x='strain', y='vel', data = df_means_sort[2:4], size=15,\
                       color='grey', edgecolor='k', marker="o",\
                       linewidth=1, dodge=True,\
                       order = ['wt', 'smy1delta'])
        
    ax = sns.stripplot(x='strain', y='vel',\
                       data = df_means_sort[4:], size=15,\
                       color='grey', edgecolor='k', marker="^",\
                       linewidth=1, dodge=True,\
                       order = ['wt', 'smy1delta'])        
        
    ax = sns.pointplot(x='strain', y='vel', data = df3,\
                       capsize = 0.8, join=False, color='k')
        
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2))            
    plt.ylabel(u'Initial rate (${\mu}m/sec$)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.tick_params('both', length=5, which='both')    
    plt.xlabel(None)
    ax.set(xticks=[])     
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.ylim([0, 1])
    plt.legend([],[], frameon=False)
    # plt.savefig('201217_smy1d_wt_cells_cables_Ro.svg') 

#=============================================================================
#use ttest to determine statistical significance    
df4 = df2.groupby(['strain', 'n']).mean()    

print(scipy.stats.ttest_ind(df4['vel'][:3], df4['vel'][3:]))    
    
    





