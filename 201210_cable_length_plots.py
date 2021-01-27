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
#Goal: Compare measured cable architecture parameters from different strains/
#mutants. 
#=============================================================================
#indicate the data directory
datadir = "D:\\Goode_Lab\\Projects\\actin_cables\\data\\summary_data\\"

#initalize data frame to append all data 
df = pd.DataFrame()

#import data to dataframe
df = pd.read_csv(datadir + '210105_cdc28-13ts_t-8_t1_yBG12_yBG9_all_cable_analysis.csv')

#calculate the ratio of length versus diameter and volume for each cable
df['len_ratio'] = df['L'] / df['cell_diameter']
df['vol_ratio'] = df['L'] / df['cell_volume']

#calculate the mean for each replicate and strain for swarm and pointplots
df_means = df.groupby(['strain', 'expt_num'], sort=False).mean().reset_index()

df_means_sort = df.groupby(['expt_num', 'strain']).mean().reset_index()

df_all_mean = df.groupby(['strain', 'cell_number'],\
                                        sort=False).mean().reset_index()

df_all_expt_mean = df.groupby(['cell_number', 'strain', 'expt_num'],\
                                        sort=False).mean().reset_index()
    
df_all_expt_mean_sort = df.groupby(['expt_num', 'strain'],\
                                       sort=True).mean().reset_index()
    
#=============================================================================
#set color palette to use for plots
cmap = ["#f6511d", "#00a6ed", "#7fb800", "#ffb400", "#0d2c54"] 

#order to plot experimental replicates
o = ['yBG12', 'yBG9', 'cdc28-13ts, t0', \
                                'cdc28-13ts, t8']
    
ft = 22 #font size for axis    

st = 'ticks' #set the style of ticks
#=============================================================================

#plot each of the measured cable parameters

#plot D
with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    sns.set_palette(cmap)
    
    sns.swarmplot(x='strain', y='D', data = df, linewidth=0.5,\
                  edgecolor='k', zorder=0, size=7, dodge=True)   
        
    ax = sns.stripplot(x='strain', y='D', data = df_means_sort[:4], size=15,\
                       color='grey', edgecolor='k', marker="s",\
                       linewidth=1, order = o)
        
    ax = sns.stripplot(x='strain', y='D', data = df_means_sort[4:8], size=15,\
                       color='grey', edgecolor='k', marker="o",\
                       linewidth=1, order = o)
        
    ax = sns.stripplot(x='strain', y='D',\
                       data = df_means_sort[8:12], size=15,\
                       color='grey', edgecolor='k', marker="^",\
                       linewidth=1, order = o)        
        
    ax = sns.pointplot(x='strain', y='D', data = df_means,\
                       capsize = 0.8, join=False, color='k')
        
    plt.ylabel(u'End-to-end distance(${\mu}m$)', fontsize=ft)
    plt.xlabel(None)
    ax.set(xticks=[])
    ax.tick_params('both', length=5, which='both')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2))            
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.ylim([0, 13])
    plt.legend([],[], frameon=False)
    # plt.savefig('201217_wt_cells_cables_D.svg')
    

#plot L   
with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    sns.set_palette(cmap)
    
    sns.swarmplot(x='strain', y='L', data = df, linewidth=0.5,\
                  edgecolor='k', zorder=0, size=7, dodge=True)   
        
    ax = sns.stripplot(x='strain', y='L', data = df_means_sort[:4], size=15,\
                       color='grey', edgecolor='k', marker="s",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='L', data = df_means_sort[4:8], size=15,\
                       color='grey', edgecolor='k', marker="o",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='L',\
                       data = df_means_sort[8:12], size=15,\
                       color='grey', edgecolor='k', marker="^",\
                       linewidth=1, dodge=True, order = o)        
        
    ax = sns.pointplot(x='strain', y='L', data = df_means,\
                       capsize = 0.8, join=False, color='k')
        
    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    plt.xlabel(None)
    ax.set(xticks=[])
    ax.tick_params('both', length=5, which='both')    
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2))            
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.ylim([0, 17])
    plt.legend([],[], frameon=False)
    # plt.savefig('201217_wt_cells_cables_L.svg') 
    

#plot L/D   
with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    sns.set_palette(cmap)
    
    sns.swarmplot(x='strain', y='L/D', data = df, linewidth=0.5,\
                  edgecolor='k', zorder=0, size=7, dodge=True)   
        
    ax = sns.stripplot(x='strain', y='L/D', \
                       data = df_means_sort[:4], size=15,\
                       color='grey', edgecolor='k', marker="s",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='L/D', \
                       data = df_means_sort[4:8], size=15,\
                       color='grey', edgecolor='k', marker="o",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='L/D', \
                       data = df_means_sort[8:12], size=15,\
                       color='grey', edgecolor='k', marker="^",\
                       linewidth=1, dodge=True, order = o)        
        
    ax = sns.pointplot(x='strain', y='L/D', data = df_means,\
                       capsize = 0.8, join=False, color='k')
        
    plt.ylabel('Tortuosity (L/D)', fontsize=ft)
    plt.xlabel(None)
    ax.set(xticks=[])
    ax.tick_params('both', length=5, which='both')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2))            
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.ylim([0, 8])
    plt.legend([],[], frameon=False)
    # plt.savefig('201217_wt_cells_tortuosity.svg')
    

#plot cell_diameter    
with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    sns.set_palette(cmap)
    
    sns.swarmplot(x='strain', y='cell_diameter', data = df_all_mean, \
                  linewidth=0.5, edgecolor='k', zorder=0, size=10, dodge=True)   
        
    ax = sns.stripplot(x='strain', y='cell_diameter', \
                       data = df_all_expt_mean_sort[:4], size=15,\
                       color='grey', edgecolor='k', marker="s",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='cell_diameter',\
                       data = df_all_expt_mean_sort[4:8], size=15,\
                       color='grey', edgecolor='k', marker="o",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='cell_diameter', \
                       data = df_all_expt_mean_sort[8:12], size=15,\
                       color='grey', edgecolor='k', marker="^",\
                       linewidth=1, dodge=True, order = o)        
        
    ax = sns.pointplot(x='strain', y='cell_diameter',\
                       data = df_all_expt_mean_sort,\
                       capsize = 0.8, join=False, color='k', order = o)
        
    plt.ylabel(u'Mother cell length (${\mu}m$)', fontsize=ft)
    plt.xlabel(None)
    ax.set(xticks=[])
    ax.tick_params('both', length=5, which='both')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2))            
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.ylim([0, 15])
    plt.legend([],[], frameon=False)
    # plt.savefig('201217_wt_cells_diameter.svg')
    

#plot cell_volume    
with sns.axes_style(st):
    plt.figure(figsize=(5.3,6))#use 3,4 for figures; 8,9 for terminal
    sns.set_palette(cmap)
    
    sns.swarmplot(x='strain', y='cell_volume', data = df_all_mean, \
                  linewidth=0.5, edgecolor='k', zorder=0, size=10, dodge=True)   
        
    ax = sns.stripplot(x='strain', y='cell_volume', \
                       data = df_all_expt_mean_sort[:4], size=15,\
                       color='grey', edgecolor='k', marker="s",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='cell_volume',\
                       data = df_all_expt_mean_sort[4:8], size=15,\
                       color='grey', edgecolor='k', marker="o",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='cell_volume', \
                       data = df_all_expt_mean_sort[8:12], size=15,\
                       color='grey', edgecolor='k', marker="^",\
                       linewidth=1, dodge=True, order = o)        
        
    ax = sns.pointplot(x='strain', y='cell_volume', \
                       data = df_all_expt_mean_sort,\
                       capsize = 0.8, join=False, color='k', order = o)
        
    plt.ylabel(u'Mother cell volume (${\mu}m^3$)', fontsize=ft)
    plt.xlabel(None)
    ax.set(xticks=[])
    ax.tick_params('both', length=5, which='both')
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.ylim([0, 350])
    plt.legend([],[], frameon=False)
    # plt.savefig('201217_wt_cells_volume.svg')
    

#plot L vs diameter    
with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    sns.set_palette(cmap)
    
    sns.swarmplot(x='strain', y='len_ratio', data = df, linewidth=0.5,\
                  edgecolor='k', zorder=0, size=7, dodge=True)   
        
    ax = sns.stripplot(x='strain', y='len_ratio', \
                       data = df_means_sort[:4], size=15,\
                       color='grey', edgecolor='k', marker="s",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='len_ratio', \
                       data = df_means_sort[4:8], size=15,\
                       color='grey', edgecolor='k', marker="o",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='len_ratio', \
                       data = df_means_sort[8:12], size=15,\
                       color='grey', edgecolor='k', marker="^",\
                       linewidth=1, dodge=True, order = o)        
        
    ax = sns.pointplot(x='strain', y='len_ratio', data = df_means,\
                       capsize = 0.8, join=False, color='k')
        
    plt.ylabel('Cable length/Mother cell length ratio', fontsize=20)
    plt.xlabel(None)
    ax.set(xticks=[])
    ax.tick_params('both', length=5, which='both')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))            
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.ylim([0, 2.4])
    plt.legend([],[], frameon=False)
    # plt.savefig('201217_wt_cells_diameter_scaling.svg')
    

#plot L vs volume    
with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    sns.set_palette(cmap)
    
    sns.swarmplot(x='strain', y='vol_ratio', data = df, linewidth=0.5,\
                  edgecolor='k', zorder=0, size=7, dodge=True)   
        
    ax = sns.stripplot(x='strain', y='vol_ratio',\
                       data = df_means_sort[:4], size=15,\
                       color='grey', edgecolor='k', marker="s",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='vol_ratio', \
                       data = df_means_sort[4:8], size=15,\
                       color='grey', edgecolor='k', marker="o",\
                       linewidth=1, dodge=True, order = o)
        
    ax = sns.stripplot(x='strain', y='vol_ratio', \
                       data = df_means_sort[8:12], size=15,\
                       color='grey', edgecolor='k', marker="^",\
                       linewidth=1, dodge=True, order = o)        
        
    ax = sns.pointplot(x='strain', y='vol_ratio', data = df_means,\
                       capsize = 0.8, join=False, color='k')
        
    plt.ylabel('Cable length/Mother cell volume ratio', fontsize=20)
    plt.xlabel(None)
    ax.set(xticks=[])
    ax.tick_params('both', length=5, which='both')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))            
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.ylim([0, 0.3])
    plt.legend([],[], frameon=False)
    # plt.savefig('201217_wt_cells_volume_scaling.svg')
    

#plot distributions of L vs diameter 
# import joypy

# fig, axes = joypy.joyplot(df, column=['len_ratio'], overlap=1.0, by="strain",\
#                           ylim='own', x_range=(0,2), fill=True,\
#                           figsize=(8,10), legend=False, xlabels=True,\
#                           ylabels=False, \
#                           alpha=1.0, linewidth=2, linecolor='k')#, fade=True)
# plt.axvline(x=1.0, lw=2, color='k')    
# plt.rc("font", size=30)
# plt.xlabel('Cable length / Cell length ratio', fontsize=32, color='k', alpha=1)
# plt.ylabel('strain', fontsize=ft, color='grey', alpha=1)
# # plt.savefig('201012_wt_cells_length_ratio_distribution.svg')

# fig, axes = joypy.joyplot(df, column=['vol_ratio'], overlap=1.0, by="strain",\
#                           ylim='own', x_range=(0,0.25), fill=True,\
#                           figsize=(8,10), legend=False, xlabels=True,\
#                           ylabels=False, \
#                           alpha=1.0, linewidth=2, linecolor='k')#, fade=True)
# # plt.axvline(x=1.0, lw=2, color='k')    
# plt.rc("font", size=30)
# plt.xlabel('Cable length / Cell volume ratio', fontsize=32, color='k', alpha=1)
# plt.ylabel('strain', fontsize=ft, color='grey', alpha=1)
# plt.savefig('201218_wt_cells_volume_ratio_distribution.svg')     

#============================================================================
#do statistical tests on the data

#organize data by stain and expt number in a new dataframe
df_cable_stats = df_means[['D', 'L', 'L/D', 'len_ratio', 'vol_ratio']]

df_cell_stats = df_all_expt_mean[['cell_diameter', 'cell_volume']]

#the order is cdc28 uninduced, cdc28 induced, yBG12, yBG9    

#compare yBG12, yBG9
stat_12_9_cable, pval_12_9_cable = scipy.stats.ttest_ind(df_cable_stats[:3],\
                                                         df_cable_stats[3:6])
    
stat_12_9_cell, pval_12_9_cell = scipy.stats.ttest_ind(df_cell_stats[:3],\
                                                         df_cell_stats[3:6])
    

#compare yBG12, cdc28_uninduced
stat_12_u_cable, pval_12_u_cable = scipy.stats.ttest_ind(df_cable_stats[:3],\
                                                         df_cable_stats[6:9])
    
stat_12_u_cell, pval_12_u_cell = scipy.stats.ttest_ind(df_cell_stats[:3],\
                                                         df_cell_stats[6:9])
    
    

#compare yBG12, cdc28_induced
stat_12_i_cable, pval_12_i_cable = scipy.stats.ttest_ind(df_cable_stats[:3],\
                                                         df_cable_stats[9:12])

stat_12_i_cell, pval_12_i_cell = scipy.stats.ttest_ind(df_cell_stats[:3],\
                                                         df_cell_stats[9:12])
    

    
#compare yBG9, cdc28_uninduced
stat_9_u_cable, pval_9_u_cable = scipy.stats.ttest_ind(df_cable_stats[3:6],\
                                                         df_cable_stats[6:9])

stat_9_u_cell, pval_9_u_cell = scipy.stats.ttest_ind(df_cell_stats[3:6],\
                                                         df_cell_stats[6:9])

    
#compare yBG9, cdc28_induced
stat_9_i_cable, pval_9_i_cable = scipy.stats.ttest_ind(df_cable_stats[3:6],\
                                                         df_cable_stats[9:12])

stat_9_i_cell, pval_9_i_cell = scipy.stats.ttest_ind(df_cell_stats[3:6],\
                                                         df_cell_stats[9:12])

    
#compare cdc28_uninduced, cdc28_induced
stat_u_i_cable, pval_u_i_cable = scipy.stats.ttest_ind(df_cable_stats[6:9],\
                                                         df_cable_stats[9:12])

stat_u_i_cell, pval_u_i_cell = scipy.stats.ttest_ind(df_cell_stats[6:9],\
                                                         df_cell_stats[9:12])











