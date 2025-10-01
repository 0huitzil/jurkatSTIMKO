#%%
"""
This file contains all the scripts necessary to create fig #1
"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from auto import run, load, save, merge, relabel, cl, klb
"""
importing changes to the RcParams
"""
from STIMKO_Options import *
from STIMKO_Models import  sampleS2KO, sampleS1KO, sampleWT, parJurkatCell, plotTS, sampleS2KOBase, sampleS1KOBase, sampleWTBase, parBaseJurkatCell
matplotlib.rcParams.update(myRcParams())

#%%
"""
Data collection
Model traces
"""
Vplc = [0.09, 0.1, 0.12]
tf = 250
# Base model
WT = []
par = parBaseJurkatCell()
for V in Vplc:
    par['Vplc']=0
    sample = sampleWTBase(par, V, tf)
    WT.append(sample)
S2KO = []
par = parBaseJurkatCell()
for V in Vplc:
    par['Vplc']=0
    sample = sampleS2KOBase(par, V, tf)
    S2KO.append(sample)
S1KO = []
par = parBaseJurkatCell()
for V in Vplc:
    par['Vplc']=0
    sample = sampleS1KOBase(par, V, tf)
    S1KO.append(sample)
#%%
#% Figure creation
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 3, nrows=3, figure=fig)
fig_width,fig_height = set_figsize(1, (1.4,1), export=True) #184 is Beamer width
# fig.suptitle(sheetnames[k] + r'-Traces' + r'')
fig.set_size_inches([fig_width,fig_height])
xlabel = ''
for i, sample in enumerate(WT):
    ax = fig.add_subplot(spec[0, i])
    title = 'A' + str(i+1)
    ylim = [0,0.5]
    color = 'tab:purple'
    xtrace = sample[0]
    ytrace = sample[1]
    xlim = [0, tf]
    yticks = ylim
    ylabel = '$c$' 
    # label = r'$V_{sM}$ = '  + str(VsM[i])
    xticks = []
    xlabel = ''
    legendsize =8
    if i==0:
        yticks = ylim
        ylabel = '$c$' 
    else:
        yticks = []
        ylabel = ''
    plotTS(ax, xtrace=xtrace, ytrace=ytrace, xlim=xlim, ylim=ylim, yticks=yticks, xticks=xticks, color=color,  ylabel=ylabel, subtitle=title, legend=False, xlabel=xlabel, ylabelpad=-10, label='', legendsize=legendsize)

for i, sample in enumerate(S1KO):
    ax = fig.add_subplot(spec[1, i])
    title = 'B' + str(i+1)
    ylim = [0,0.5]
    color = 'tab:red'
    xtrace = sample[0]
    ytrace = sample[1]
    xlim = [0, tf]
    yticks = []
    ylabel = ''
    if i==0:
        yticks = ylim
        ylabel = '$c$' 
    else:
        yticks = []
        ylabel = ''
    plotTS(ax, xtrace=xtrace, ytrace=ytrace, xlim=xlim, ylim=ylim, yticks=yticks, xticks=xticks, color=color,  ylabel=ylabel, subtitle=title, legend=False, xlabel=xlabel)

for i, sample in enumerate(S2KO):
    ax = fig.add_subplot(spec[2, i])
    title = 'C' + str(i+1)
    ylim = [0,0.5]
    color = 'tab:blue'
    xtrace = sample[0]
    ytrace = sample[1]
    xlim = [0, tf]
    # xlabel = 'Time (s) \n\n' + r'$V_{\mathrm{PLC}} =$ ' + str(Vplc[i])
    yticks = []
    ylabel = ''
    xlabel = 'Time (s) \n\n' + r'$V_{\mathrm{PLC}} =$ ' + str(Vplc[i])
    xticks = xlim
    if i==0:
        yticks = ylim
        ylabel = '$c$' 
    else:
        yticks = []
        ylabel = ''
    plotTS(ax, xtrace=xtrace, ytrace=ytrace, xlim=xlim, ylim=ylim, yticks=yticks, xticks=xticks, color=color,  ylabel=ylabel, xlabel=xlabel, subtitle=title, legend=False)

filename = 'Fig_baseModelOscillations.pdf'
# saveFigure(filename, fig)
# %%
