# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 12:33:20 2020

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
import seaborn as sns
import math
from math import sqrt
from scipy.spatial import distance

#import os to handle operating system
import os
#============================================================================= 
#Goal: Import individual cable trajectories and calculate distance, velocity, 
#total time, length, etc. Save these to a csv for each data set. 

#Note: Must change the frame rate (r) to match the imaging session!
#==============================================================================
#setup the data directory
datadir = "D:\\Goode_Lab\\Projects\\actin_cables\\data\\cable_trajectory_data\\test\\"

#initalize data frame to append all data 
df = pd.DataFrame()

df2 = pd.DataFrame()

df3 = pd.DataFrame()

r = 2.3 #frame interval; must set for each microscopy session

#==============================================================================

#Read in each trace from the directory
for f in glob.glob(datadir + '*' + '.csv'):
    
    #open each file
    df = pd.read_csv(f)
    #initialize an empty dataframe to write calculations to
    df2 = pd.DataFrame()
    #calculate the differences in x-coordinates
    df2['deltx'] = list(map(operator.sub, df['X'][:-1], df['X'][1:]))
    #calculate the differences in y-coordinates
    df2['delty'] = list(map(operator.sub, df['Y'][:-1], df['Y'][1:]))
    #use distance formula to calculate distance traveled per frame
    df2['dist'] = np.sqrt(df2['deltx']**2 + df2['delty']**2)
    #convert distance from pixels to microns
    df2['dist_um'] = df2['dist'] * 0.133
    #calculate cumulative distance travelled
    df2['cum_dist'] = np.cumsum(df2['dist_um'])
    #calculate the relative distance from the bud neck
    df2['rel_dist'] = df2['cum_dist'] / df2['cum_dist'].max()
    #calculate the absoulte difference from rear of mother cell
    df2['rear_dist'] = abs(df2['cum_dist'] - df2['cum_dist'].max())
    #calculate the length from the origin
    df2['neck_dist'] = abs(df2['cum_dist'] - df2['cum_dist'].min())
    #calculate velocity from distance and time variables
    df2['vel'] = df2['dist_um'] / r
    #write the frame number to the dataframe
    df2['frame'] = df['Index']
    #calculate cable lifetime from frame and interval (~2.3sec)
    df2['lifetime'] = df['Index'] * r
    #add file name for plotting of individual traces
    df2['cell'] = os.path.basename(f)
    
    #plot each trajectory 
    sns.lineplot(x=df2['lifetime'], y=df2['vel'], hue=df2['cell'], \
              palette='Set1', ci=95)       
    plt.ylim([0, 1])
    plt.xlim([0, 60])        

    plt.show()    
    
    #append data from all traces into a single dataframe
    df3 = df3.append(df2)
    
    
# output the data to csv
# outputdir = datadir + '\\summary_data\\'
# df3.to_csv(datadir + '\\200826_yBG12-Abp140Envy_cable_trajectory_analysis.csv',\
#             index=False)


   

