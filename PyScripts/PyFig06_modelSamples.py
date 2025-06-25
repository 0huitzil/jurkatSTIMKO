#%%
"""
Fig06
Sample time series of the model
"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from auto import run, load, save, merge, relabel, cl, klb
"""
importing changes to the RcParams
"""
from myOptions import *
from PyModels import  sampleS2KO, sampleS1KO, sampleWT, parJurkatCell, plotTS
matplotlib.rcParams.update(myRcParams())
#%% Data collection
"""
Collecting data for the WT, S1 (S2KO), S2 (S1KO) versions of the model
"""
Vplc = [0.06, 0.1, 0.15, 0.2, 0.25]
par = parJurkatCell()
tf = 250
WT = []
for V in Vplc:
    par['Vplc']=0
    sample = sampleWT(par, V, tf)
    WT.append(sample)
S1 = []
for V in Vplc:
    par['Vplc']=0
    sample = sampleS2KO(par, V, tf)
    S1.append(sample)
S2 = []
for V in Vplc:
    par['Vplc']=0
    sample = sampleS1KO(par, V, tf)
    S2.append(sample)
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 5, nrows=3, figure=fig)
fig_width,fig_height = set_figsize(1, (1.6,1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
Plotting the WT time series
"""
for i, sample in enumerate(WT):
    ax = fig.add_subplot(spec[0, i])
    title = 'A' + str(i+1)
    ylim = [0,0.6]
    color = 'tab:purple'
    xtrace = sample[0]
    ytrace = sample[1]
    xlim = [0, tf]
    if i==0:
        yticks = ylim
        ylabel = '$c$'
    else:
        yticks = []
        ylabel = ''
    plotTS(ax, xtrace=xtrace, ytrace=ytrace, xlim=xlim, ylim=ylim, yticks=yticks, color=color,  ylabel=ylabel, subtitle=title, legend=False)
"""
Plotting the S2 (S1KO) time series
"""
for i, sample in enumerate(S2):
    ax = fig.add_subplot(spec[1, i])
    title = 'B' + str(i+1)
    ylim = [0,0.6]
    color = 'tab:red'
    xtrace = sample[0]
    ytrace = sample[1]
    xlim = [0, tf]
    if i==0:
        yticks = ylim
        ylabel = '$c$'
    else:
        yticks = []
        ylabel = ''
    plotTS(ax, xtrace=xtrace, ytrace=ytrace, xlim=xlim, ylim=ylim, yticks=yticks, color=color,  ylabel=ylabel, subtitle=title, legend=False)
"""
Plotting the S1 (S2KO) time series
"""
for i, sample in enumerate(S1):
    ax = fig.add_subplot(spec[2, i])
    title = 'C' + str(i+1)
    ylim = [0,0.6]
    color = 'tab:blue'
    xtrace = sample[0]
    ytrace = sample[1]
    xlim = [0, tf]
    xlabel = 'Time (s) \n\n' + r'V$_{plc} = $ ' + str(Vplc[i]) + r' $\mu$M'
    if i==0:
        yticks = ylim
        ylabel = '$c$'
    else:
        yticks = []
        ylabel = ''
    plotTS(ax, xtrace=xtrace, ytrace=ytrace, xlim=xlim, ylim=ylim, yticks=yticks, xticks=xlim, color=color,  ylabel=ylabel, xlabel=xlabel, subtitle=title, legend=False)

filename = 'Fig6.pdf'
saveFigure(filename, fig)