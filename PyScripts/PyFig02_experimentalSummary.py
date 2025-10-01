#%% Libraries 
"""
Fig01 
Experimental summary and representative patterns
"""
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
"""
importing changes to the RcParams
"""
from STIMKO_Options import *
matplotlib.rcParams.update(myRcParams())
#%% Data import
"""
Data for the histogram is taken from the STIMSUMMARY.csv file
Data for the time series is taken from the STIMSAMPLES.csv file
"""
oscNames = [
    'Single Release', 
    'Wide Spike', 
    'Burst', 
    'Narrow Spike', 
]
data = pd.read_csv('experimentalData/STIMSummary.csv', index_col=0)
STIMsamples = pd.read_csv('experimentalData/STIMsamples.csv', index_col=0)
data.columns = oscNames #Rewriting colnames to easy plotting
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 1, nrows=3, figure=fig)
subSpec = spec[1:3,0].subgridspec(2,2) #Separating the ax for the histogram with the time series
fig_width,fig_height = set_figsize(1, (1.5,1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
Plotting histogram
"""
ax = fig.add_subplot(spec[0,0])
data.plot.bar(stacked=True, ax=ax, rot=0, backend='matplotlib')
# Ax settings
ax.set_ylim(0, 140)
ax.set_yticks([0 ,50, 100, 140])
ax.set_xlabel('')
ax.set_ylabel('\# \n cells', rotation='horizontal', labelpad=0)
ax.set_title('A', loc='left')
ax.legend( loc='best', ncol=5, bbox_to_anchor=(0, -0.35, 1, 0.2), mode='expand')
"""
Plotting single release sample
"""
ylim = [0, 0.7]
yticks = [0, 0.7]
ax = fig.add_subplot(subSpec[0])
STIMsamples.plot(y='Single Release', ax=ax, color='C0')
ax.set_xlim([0, 450])
ax.set_xticks([])
ax.set_xlabel('')
ax.set_ylabel(r'F$_{340}$/F$_{380}$')
ax.set_ylim(ylim)
ax.set_yticks(yticks)
ax.set_title('B', loc='left')
"""
Plotting wide spike time series
"""
ax = fig.add_subplot(subSpec[1])
STIMsamples.plot(y='Wide Spike', ax=ax, color='C1')
ax.set_xlim([0, 450])
ax.set_xticks([])
ax.set_xlabel('')
ax.set_ylim(ylim)
ax.set_yticks([])
ax.set_title('C', loc='left')
"""
Plotting bursting time series
"""
ax = fig.add_subplot(subSpec[2])
STIMsamples.plot(y='Burst', ax=ax, color='C2')
ax.set_xlim([0, 450])
ax.set_xticks([0, 450])
ax.set_xlabel('Time (seconds)')
ax.set_ylabel(r'F$_{340}$/F$_{380}$')
ax.set_ylim(ylim)
ax.set_yticks(yticks)
ax.set_title('D', loc='left')
"""
Plotting narrow spike time series
"""
ax = fig.add_subplot(subSpec[3])
STIMsamples.plot(y='Narrow Spike', ax=ax, color='C3')
ax.set_xlim([0, 450])
ax.set_xticks([0, 450])
ax.set_xlabel('Time (seconds)')
ax.set_ylim(ylim)
ax.set_yticks([])
ax.set_title('E', loc='left')

filename = 'Fig2.pdf'
saveFigure(filename, fig)
# %%
