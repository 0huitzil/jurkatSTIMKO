#%% Libraries 
"""
Fig01 
Experimental time series
"""
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
"""
importing changes to the RcParams
"""
from STIMKO_Options import *
matplotlib.rcParams.update(myRcParams())
from STIMKO_Models import plotDataSample

#%% Data import
"""
Data is taken from the samples.csv files
"""
colors = ['tab:purple', 'tab:red', 'tab:blue']
WT = pd.read_csv('experimentalData/WTsamples.csv', index_col=0)
S1 = pd.read_csv('experimentalData/S1samples.csv', index_col=0)
S2 = pd.read_csv('experimentalData/S2samples.csv', index_col=0)
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 5, nrows=3, figure=fig)
fig_width,fig_height = set_figsize(1, (1.4,1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
PLotting WT samples
"""
for i, column in enumerate(WT.columns):
    ax = fig.add_subplot(spec[0, i])
    title = 'A' + str(i+1)
    ylim = [0,1]
    if i==0:
        yticks = [0, 0.9]
        ylabel = r'F$_{340}$/F$_{380}$'
    else:
        yticks = []
        ylabel = ''
    plotDataSample(df=WT, colname=column, ax=ax, color=colors[0], title=title, ylim=ylim, yticks=yticks, ylabel=ylabel)
"""
PLotting S1KO samples
"""
for i, column in enumerate(S2.columns):
    ax = fig.add_subplot(spec[1, i])
    title = 'B' + str(i+1)
    ylim = [0,0.6]
    if i==0:
        yticks = ylim
        ylabel = r'F$_{340}$/F$_{380}$'
    else:
        yticks = []
        ylabel = ''
    plotDataSample(df=S2, colname=column, ax=ax, color=colors[1], title=title, ylim=ylim, yticks=yticks, ylabel=ylabel)
"""
PLotting S2KO samples
"""
for i, column in enumerate(S1.columns):
    ax = fig.add_subplot(spec[2, i])
    title = 'C' + str(i+1)
    ylim = [0,0.7]
    if i==0:
        yticks = ylim
        ylabel = r'F$_{340}$/F$_{380}$'
    else:
        yticks = []
        ylabel = ''
    xticks = [0, 450]
    xlabel = 'Time (s)'
    plotDataSample(df=S1, colname=column, ax=ax, color=colors[2], title=title, ylim=ylim, xticks=xticks, xlabel=xlabel, yticks=yticks, ylabel=ylabel)

filename = 'Fig1.pdf'
saveFigure(filename, fig)
# %%
