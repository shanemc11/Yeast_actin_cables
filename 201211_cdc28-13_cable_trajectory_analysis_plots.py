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
df_t0 = pd.DataFrame()

df_t8 = pd.DataFrame()

all_df = pd.DataFrame()


#read in the summary data files to compare t0 to t8 cells
df_t0 = pd.read_csv(datadir + \
                    "200826_t0_all_cable_extension_rate_analysis_cutoff.csv")

df_t8 = pd.read_csv(datadir + \
                    "200826_t8_all_cable_extension_rate_analysis_cutoff.csv")
    
#combine data into a single dataframe for some of the plotting/stats
frames = [df_t0, df_t8]

all_df = pd.concat(frames)    
    
#=============================================================================    
#calculate means and std deviation for each time point

df_t0_t_mean = pd.DataFrame()
df_t0_t_mean = df_t0.groupby(['lifetime']).mean().reset_index()  

df_t0_t_std = pd.DataFrame()
df_t0_t_std =  df_t0.groupby(['lifetime']).std().reset_index()

df_t8_t_mean = pd.DataFrame()
df_t8_t_mean = df_t8.groupby(['lifetime']).mean().reset_index()

df_t8_t_std = pd.DataFrame()
df_t8_t_std =  df_t8.groupby(['lifetime']).std().reset_index()

#calculate means for each timepoint for each replicate
df_t0_expt_mean = df_t0.groupby(['lifetime', 'n'],\
                                       sort=False).mean().reset_index()
    
df_t8_expt_mean = df_t8.groupby(['lifetime', 'n'],\
                                       sort=False).mean().reset_index()
 
#=============================================================================
#initialize plotting parameters
cmap = ["#7fb800", "#ffb400"] 
ft = 22 #font size for axis
ft2 = 30 #font size for axis
t_tot = np.linspace(0,60,61) #range of times to plot
o = ['t0', 't8'] #order to plot initial rate
st = 'ticks' #set the style of ticks

#=============================================================================
#plot the extension rates using the mean of each replicate 

with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    #plot the mean and 95%CI for the replicates   
    ax = sns.lineplot(x=df_t0_expt_mean['lifetime'],\
                      y=df_t0_expt_mean['vel'],
                      color='#7fb800', ci=95, label='Mean', lw=3) 
    #plot the mean for each replicates          
    ax = sns.scatterplot(x=df_t0_expt_mean['lifetime'], \
                         y=df_t0_expt_mean['vel'],
                         color = 'grey', label='Experiment', edgecolor='k',\
                         linewidth=1, alpha=1, s=80) 
        
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))     
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 0.5])
    plt.xlim([0, 60])        
    plt.tight_layout()
    # plt.savefig('201217_uninduced_cable_ext_vs_lifetime_exptN.svg')  
 

with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    #plot the mean and 95%CI for the replicates   
    ax = sns.lineplot(x=df_t8_expt_mean['lifetime'],\
                      y=df_t8_expt_mean['vel'],
                      color='#ffb400', ci=95, label='Mean', lw=3) 
    #plot the mean for each replicates       
    ax = sns.scatterplot(x=df_t8_expt_mean['lifetime'], \
                         y=df_t8_expt_mean['vel'],
                         color = 'grey', label='Experiment', edgecolor='k',\
                         linewidth=1, alpha=1, s=80) 
        
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.tick_params('both', length=5, which='both')         
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 0.5])
    plt.xlim([0, 60])        
    plt.tight_layout()
    # plt.savefig('201217_induced_cable_ext_vs_lifetime_exptN.svg')     

#=============================================================================
#plot the change in length as a function of time      
    
with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    sns.set_palette(cmap)   
    #plot each cable   
    ax = sns.scatterplot(x=df_t8_expt_mean['lifetime'], \
                         y=df_t8_expt_mean['neck_dist'], \
                             color='#ffb400',alpha=1, \
                                 linewidth=0.7, edgecolor='k')
  
    ax = sns.scatterplot(x=df_t0_expt_mean['lifetime'],\
                         y=df_t0_expt_mean['neck_dist'], \
                             color='#7fb800',alpha=1,\
                                 linewidth=0.7, edgecolor='k')
        
    #plot the mean and 95%CI for the replicates   
    ax = sns.lineplot(x=df_t0_expt_mean['lifetime'], \
                      y=df_t0_expt_mean['neck_dist'], \
                          ci=95, label='cdc28-13, uninduced',\
                              color='#7fb800', lw=3)
    
    ax = sns.lineplot(x=df_t8_expt_mean['lifetime'], \
                      y=df_t8_expt_mean['neck_dist'],\
                      color='#ffb400', ci=95, label = 'cdc28-13, induced',\
                          lw=3)
 
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))        
    ax.tick_params('both', length=5, which='both')
    ax.get_legend().remove()    
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 15])
    plt.xlim([0, 60])        
    plt.tight_layout()
    # plt.savefig('201217_cdct0_cdct8_len_vs_lifetime.svg') 


#=============================================================================
#plot the change in rate as a function of length and the ratio of lengths
#use the time averaged lengths and rates for this or the plots don't plot well
#due to the different values for each length

#calculate the errors (95%CI) for each strain   
t0_lower_ci = df_t0_t_mean - (1.96 * (df_t0_t_std / np.sqrt(64)))
t0_upper_ci = df_t0_t_mean + (1.96 * (df_t0_t_std / np.sqrt(64)))

t8_lower_ci = df_t8_t_mean - (1.96 * (df_t8_t_std / np.sqrt(57)))
t8_upper_ci = df_t8_t_mean + (1.96 * (df_t8_t_std / np.sqrt(57)))


#plot the change in extension rate as a function of the relative length
with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    sns.set_palette(cmap)   
            
    ax = sns.lineplot(x=df_t0_t_mean['neck_dist']/4.9, y=df_t0_t_mean['vel'],\
              ci=95, label='cdc28-13, uninduced', color='#7fb800', lw=3)
    
    ax = sns.lineplot(x=df_t8_t_mean['neck_dist']/8.7, y=df_t8_t_mean['vel'],\
                      color='#ffb400', ci=95, label = 'cdc28-13, induced',\
                          lw=3)

    plt.fill_between(df_t0_t_mean['neck_dist']/4.9, t0_lower_ci['vel'],\
                     t0_upper_ci['vel'],\
                     color='#7fb800', alpha=0.3) 
        
    plt.fill_between(df_t8_t_mean['neck_dist']/8.7, t8_lower_ci['vel'],\
                     t8_upper_ci['vel'],\
                     color='#ffb400', alpha=0.3)        
        
 
    plt.xlabel('Cable length / Mother cell length', fontsize=ft-6)
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))        
    ax.tick_params('both', length=5, which='both')    
    ax.get_legend().remove()
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 0.4])
    plt.xlim([0, 1.1])        
    plt.tight_layout()
    # plt.savefig('201217_cdct0_cdct8_extrate_v_len_norm.svg')
    
#plot the change in extension rate as a function of the length
with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    sns.set_palette(cmap)   
        
    ax = sns.lineplot(x=df_t0_t_mean['neck_dist'], y=df_t0_t_mean['vel'], \
              ci=95, label='cdc28-13, uninduced', color='#7fb800', lw=3)
    
    ax = sns.lineplot(x=df_t8_t_mean['neck_dist'], y=df_t8_t_mean['vel'],\
                      color='#ffb400', ci=95, label = 'cdc28-13, induced',\
                          lw=3)

    plt.fill_between(df_t0_t_mean['neck_dist'], t0_lower_ci['vel'],\
                     t0_upper_ci['vel'],\
                     color='#7fb800', alpha=0.3) 
        
    plt.fill_between(df_t8_t_mean['neck_dist'], t8_lower_ci['vel'],\
                     t8_upper_ci['vel'],\
                     color='#ffb400', alpha=0.3)        
        
        
 
    plt.xlabel(u'Cable length (${\mu}m$)', fontsize=ft)
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))        
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))        
    ax.tick_params('both', length=5, which='both')        
    ax.get_legend().remove()
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 0.4])
    plt.xlim([0, 10])        
    plt.tight_layout()
    # plt.savefig('201217_cdct0_cdct8_extrate_v_len.svg')
    

#plot the change in extension rate as a function of the theoretical time to 
#reach the end of the cell compartment

with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    sns.set_palette(cmap)   
        
    ax = sns.scatterplot(x=df_t8_expt_mean['lifetime']/(8.7/0.32),\
                         y=df_t8_expt_mean['neck_dist']/8.7, \
              color='#ffb400',alpha=1,linewidth=0.7, edgecolor='k')
  
    ax = sns.scatterplot(x=df_t0_expt_mean['lifetime']/(4.9/0.35),\
                         y=df_t0_expt_mean['neck_dist']/4.9, \
              color='#7fb800',alpha=1,linewidth=0.7, edgecolor='k')
    
    ax = sns.lineplot(x=df_t0_expt_mean['lifetime']/(4.9/0.35),\
                      y=df_t0_expt_mean['neck_dist']/4.9, \
              ci=95, label='cdc28-13, uninduced', color='#7fb800', lw=3)
    
    ax = sns.lineplot(x=df_t8_expt_mean['lifetime']/(8.7/0.32),\
                      y=df_t8_expt_mean['neck_dist']/8.7,\
                      color='#ffb400', ci=95, label = 'cdc28-13, induced',\
                          lw=3)
 
    plt.xlabel('Time / Time$_{max}$', fontsize=ft)    
    plt.ylabel('Cable length / Mother cell length', fontsize=ft-6)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))        
    ax.tick_params('both', length=5, which='both')        
    ax.get_legend().remove()
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 1.6])
    plt.xlim([0, 2.4])        
    plt.tight_layout()
    # plt.savefig('201217_cdct0_cdct8_rellen_vs_tmax_uniqueR0.svg')


#plot the change in extension rate as a function of time
with sns.axes_style(st):
    plt.figure(figsize=(5,5))
    sns.set_palette(cmap)   
        
    ax = sns.scatterplot(x=df_t8_expt_mean['lifetime'],\
                         y=df_t8_expt_mean['neck_dist']/8.7, \
              color='#ffb400',alpha=1,linewidth=0.7, edgecolor='k')
  
    ax = sns.scatterplot(x=df_t0_expt_mean['lifetime'],\
                         y=df_t0_expt_mean['neck_dist']/4.9, \
              color='#7fb800',alpha=1,linewidth=0.7, edgecolor='k')
    
    ax = sns.lineplot(x=df_t0_expt_mean['lifetime'],\
                      y=df_t0_expt_mean['neck_dist']/4.9, \
              ci=95, label='cdc28-13, uninduced', color='#7fb800', lw=3)
    
    ax = sns.lineplot(x=df_t8_expt_mean['lifetime'],\
                      y=df_t8_expt_mean['neck_dist']/8.7,\
                      color='#ffb400', ci=95, label = 'cdc28-13, induced',\
                          lw=3)
 
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel('Cable length / Mother cell length', fontsize=ft-6)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))        
    ax.tick_params('both', length=5, which='both')        
    ax.get_legend().remove()
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.ylim([0, 1.6])
    # plt.xlim([0, ])        
    plt.tight_layout()
    # plt.savefig('201217_cdct0_cdct8_rellen_vs_time.svg')    


#=============================================================================
#now fit the initial extension rates with linear regression    

#first fit the FP data
t_fit = 5 #final frame to fit, ~10s
  
from scipy.stats import linregress

#fit the uninduced cells
slope_t0, intercept_t0, r_value_t0, p_value_t0, std_err_t0 = \
scipy.stats.linregress(df_t0_t_mean['lifetime'][:t_fit], \
                       df_t0_t_mean['vel'][:t_fit])
      
print("r-squared_t0:", r_value_t0**2, "slope_t0:", slope_t0)     

#fit the induced cells
slope_t8, intercept_t8, r_value_t8, p_value_t8, std_err_t8 = \
scipy.stats.linregress(df_t8_t_mean['lifetime'][:t_fit-1], \
                       df_t8_t_mean['vel'][:t_fit-1])
      
print("r-squared_t8:", r_value_t8**2, "slope_t8:", slope_t8) 

#=============================================================================
#plot the fit from linear regression over the initial extension rates
with sns.axes_style(st):
    plt.figure(figsize=(6,7))
    
    
    plt.plot(t_tot,(intercept_t0 + slope_t0*t_tot),\
             'k--', lw=3,\
             label=r"Slope = {0:.3f}+/-{2:.3f}, R$^2$ = {1:.2f}".\
                 format(slope_t0, r_value_t0**2,1.96*std_err_t0 )) 

        
    sns.scatterplot(x=df_t0_t_mean['lifetime'], y=df_t0_t_mean['vel'], \
                    color='#7fb800', s = 300, ci=95, linewidth=0.5,\
                    label=None, edgecolor='k')
        
    ax =  sns.lineplot(x=df_t0['lifetime'], y=df_t0['vel'], \
              color='#7fb800', ci=95, label=None,\
                  lw=0)
    
    plt.xlabel('Extension time (sec)', fontsize=ft2)    
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft2)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.tick_params('both', length=5, which='both')                                        
    plt.rc('xtick', labelsize=ft2)
    plt.rc('ytick', labelsize=ft2)
    plt.ylim([0.1, 0.4])
    plt.xlim([0, 15])        
    plt.tight_layout()
    plt.legend(loc='upper right', prop={'size': 13})    
    # plt.savefig('201217_t0_cable_ext_linearfit.svg')
    
with sns.axes_style(st):
    plt.figure(figsize=(6,7))
       
    plt.plot(t_tot,(intercept_t8 + slope_t8*t_tot),\
             'k--', lw=3,\
             label=r"Slope = {0:.3f}+/-{2:.3f}, R$^2$ = {1:.2f}".\
                 format(slope_t8,r_value_t8**2,1.96*std_err_t8 )) 

        
    sns.scatterplot(x=df_t8_t_mean['lifetime'], y=df_t8_t_mean['vel'], \
                    color='#ffb400', s = 300, ci=95, linewidth=0.5,\
                    label=None, edgecolor='k')
        
    ax =  sns.lineplot(x=df_t8['lifetime'], y=df_t8['vel'], \
              color='#ffb400', ci=95, label=None,\
                  lw=0)
        
    plt.xlabel('Extension time (sec)', fontsize=ft2)    
    plt.ylabel(u'Extension rate (${\mu}m$/sec)', fontsize=ft2)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))                    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.tick_params('both', length=5, which='both')                                        
    plt.rc('xtick', labelsize=ft2)
    plt.rc('ytick', labelsize=ft2)
    plt.ylim([0.1, 0.4])
    plt.xlim([0, 15])        
    plt.tight_layout()
    plt.legend(loc='upper right', prop={'size': 13})    
    # plt.savefig('201217_t8_cable_ext_linearfit.svg')    
    
#=============================================================================
#setup of the inital rates for each strain
df2 = pd.DataFrame()

df2 = all_df.loc[(all_df['frame']== 0)]

df3 = df2.groupby(['time', 'n']).mean().reset_index()  
df3 = df3.sort_values('time', ascending=True)

df_means_sort = df2.groupby(['n', 'time']).mean().reset_index()

#plot the initial extension rates
with sns.axes_style(st):
    plt.figure(figsize=(5,6))
    sns.set_palette(cmap)
    
    sns.swarmplot(x='time', y='vel', data = df2, \
                  linewidth=0.5,\
                  alpha=1, edgecolor='k', size=10, dodge=True)  
        
    ax = sns.stripplot(x='time', y='vel', data = df_means_sort[:2], size=15,\
                       color='grey', edgecolor='k', marker="s",\
                       linewidth=1, dodge=True, \
                       order = o)
        
    ax = sns.stripplot(x='time', y='vel', data = df_means_sort[2:4], size=15,\
                       color='grey', edgecolor='k', marker="o",\
                       linewidth=1, dodge=True,\
                       order = o)
        
    ax = sns.stripplot(x='time', y='vel',\
                       data = df_means_sort[4:6], size=15,\
                       color='grey', edgecolor='k', marker="^",\
                       linewidth=1, dodge=True,\
                       order = o)  
        
    ax = sns.stripplot(x='time', y='vel',\
                        data = df_means_sort[6:], size=15,\
                        color='grey', edgecolor='k', marker="d",\
                        linewidth=1, dodge=True,\
                        order = o)        
        
    ax = sns.pointplot(x='time', y='vel', data = df3,\
                       capsize = 0.8, join=False, color='k')
        
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2))            
    plt.ylabel(u'Initial rate (${\mu}m/sec$)', fontsize=ft)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))     
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=ft)
    plt.xlabel(None)
    ax.set(xticks=[])
    ax.tick_params('both', length=5, which='both')        
    plt.tight_layout()
    plt.ylim([0, 1])
    plt.legend([],[], frameon=False)
    # plt.savefig('201217_cdc28-13_cables_Ro.svg')   
     
#use ttest to determine statistical significance
print(scipy.stats.ttest_ind(df3['vel'][:4], df3['vel'][4:]))    

