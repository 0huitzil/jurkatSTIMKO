#%% Libraries
"""
Fig07
Bifurcation diagram and sample time series of the S2 (S1KO) model
"""
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

"""
The path to the AUTO library is setup in myOptions.py
Make sure to set the correct path there first before running these files
"""
from myOptions import auto_directory
sys.path.append(auto_directory)
from auto import *
from auto import run, load, save, merge, relabel, cl, klb, loadbd
"""
This command allows me to export the files directly to the 
sister LaTeX directory. Feel free to comment
"""
from pathlib import Path
parentPath = str(Path(os.getcwd()).parent)
sys.path.append(parentPath)
latexPath = Path(os.getcwd()).parent/'Latex'
"""
importing changes to the RcParams
"""
from myOptions import *
from PyModels import parJurkatCell, sampleS2Stable, plotTS
matplotlib.rcParams.update(myRcParams())

#%% Data Collection
"""
Data Vplc Diag 
Only run this is if the .S2 files have not been created yet
"""
# file = "AUTOJurkatS2Cell"
# model = load(file) 
# Vplc = [0.044, 0.07, 0.1, 0.16]
# eqS2 = run(
#     model, 
#     IPS=1, 
#     ICP=['Vplc'], 
#     NMX=20000,
#     DS=1e-2,
#     DSMAX=5e-2,
#     UZSTOP={'Vplc': 2}
# )
# #%
# cycleS2 = run(
#     eqS2('HB')[0], 
#     IPS=2,
#     ICP=['Vplc', 11], 
#     NMX=140000,
#     NTST=1500,
#     DS=1e-2,
#     DSMAX=3e-2,
#     DSMIN=1e-6,
#     SP=['LP0', 'UZ', 'PD0', 'TR0'],
#     UZSTOP={'Vplc': 1}, 
#     UZR={'Vplc': Vplc}, 
#     THL={1: 1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:0.1}
# )
# save(eqS2+cycleS2, 'S2')
# cl()
#%% Data import
"""
Loading S2 AUTO file
"""
eqS2 = loadbd('S2')[0]
cycleS2 = loadbd('S2')[1]
#%% Data import - time series
"""
Data Simulation - S2 cells
"""
Vplc = [0.044, 0.07, 0.1, 0.16]
par = parJurkatCell()
tf = 250
S2 = []
for V in Vplc:
    par['Vplc']=0
    sample = sampleS2Stable(par, V, tf)
    S2.append(sample)
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 4, nrows=5, figure=fig)
fig_width,fig_height = set_figsize(1, (1.9,1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
Bifurcation diagram
"""
ax = fig.add_subplot(spec[0:2,:])
# Extracting the data from the bifurcation diagram object to plot
eqCurve = eqS2
eqCurve = eqCurve[0::1]
eqStab = eqCurve.stability()
eqStab.insert(0,0)
cyCurve = cycleS2
cyCurve = cyCurve[0:109500] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
#Alphabet
alph = [chr(item) for item in range(ord("B"), ord("Z") + 1)]
#  Equilibrium line 
for i in range(1,len(eqStab)):
    if eqStab[i]< 0: #Stable curve is plotted as solid line
        ax.plot(
            eqCurve[max(np.abs(eqStab[i-1])-1,0):np.abs(eqStab[i])]['Vplc'], 
            eqCurve[max(np.abs(eqStab[i-1])-1,0):np.abs(eqStab[i])]['c'], 
            'k', 
            ls='solid', 
            alpha=1, 
            zorder=-1,
            # label='EQ points'
        )
    else: #unstable curve is plotted as dashed line
        ax.plot(
            eqCurve[max(np.abs(eqStab[i-1])-1,0):np.abs(eqStab[i])]['Vplc'], 
            eqCurve[max(np.abs(eqStab[i-1])-1,0):np.abs(eqStab[i])]['c'], 
            'k',
            ls='dashed',
            alpha=1,
            zorder=-1,
        )
#Limit cycle curve
for i in range(1,len(cycleStab)):
    if cycleStab[i]< 0:
        ax.plot( #Stable curve is plotted as solid line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:red', 
            ls='solid', 
            zorder=-1,
            # label='Limit cycle'
        )
    else:
        ax.plot( #unstable curve is plotted as dashed line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:red',
            ls='dotted',
        )
# Labeling of the HB points
for i in eqCurve('HB').getLabels(): 
    point = eqCurve.getLabel(i)
    if i==eqCurve('HB').getLabels()[0]:
        label = 'HB'
    else:
        label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'o', 
        color =  'tab:red', 
        label = label,
        # facecolors = 'none',
        zorder=1,
            )
# Labeling of the UZR points
j=0
for i in cyCurve('UZ').getLabels()[0:4]: #Labeling of the UZR points
    point = cyCurve.getLabel(i)
    if point['Vplc'] < 0.5:
        ax.scatter(
            point['Vplc'],
            max(point['c']),
            marker = 'o', 
            color =  'C7',
            s=10,
            zorder=3, #High zorder to ensure the point shows above the line
            )
        ax.annotate(
            alph[j],
            [point['Vplc'],max(point['c'])*0.9],
        )
        j=j+1
# Ax limits 
xlim = [0,0.18]
ylim = [-0.01, 0.6]
yticks = [0, 0.6]
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks(yticks)
ax.set_xticks(xlim)
ax.set_xlabel(r'$V_{\mathrm{PLC}}$',)
ax.set_ylabel(r'$c$', rotation='horizontal')
ax.set_title('A', loc='left')
ax.legend()

"""
Time series
"""
# Ax limits (common)
xlim = [0, 250]
xticks = []
yclim = [0,0.6]
ycelim = [330, 380]
ycmlim = [-10, 400]
ycmticks = [0, 400]
yticks = []
panels = [r'$c$', r'$c_e$', r'$c_m$']
color = 'tab:red'
"""
Long function ahead
I need to plot 12 time series, 4 for c, ce and cm 
This is done in this big, ugly loop
Someone smarter than me can do it more elegantly I am sure 
"""
for col, sample in enumerate(S2):
    t = sample[0] #Time array
    for row in range(3):
        ax = fig.add_subplot(spec[row+2, col]) #First two rows are for the big diag
        trace = sample[row+1] #c, ce or cm time series
        subtitle = alph[col] + str(row+1)
        # All these if's control the plotting of the axes, axlimits and ylabel 
        # A lot of individual cases. It is a big plot
        if row==0:
            ylim = yclim
            xticks = []
            xlabel = ''
            if col==0:
                yticks=ylim
            else:
                yticks = []
        if row==1:
            ylim = ycelim
            xticks = []
            xlabel = ''
            if col==0:
                yticks=ylim
            else:
                yticks = []
        if row==2:
            ylim = ycmlim
            xticks = xlim
            xlabel = 'Time (s) \n\n' + r'$V_{\mathrm{PLC}}$ = ' + str(Vplc[col]) + r' $\mu$M/s'
            if col==0:
                yticks=ycmticks
            else:
                yticks = []
        if col==0:
            ylabel = panels[row]
        else:
            ylabel=''
        plotTS(
            ax,
            xtrace = t, 
            ytrace = trace, 
            xlim = xlim, 
            xticks = xticks, 
            ylim = ylim, 
            yticks = yticks, 
            label = label, 
            color = color,  
            xlabel= xlabel, 
            ylabel = ylabel, 
            subtitle = subtitle, 
        )

filename = 'Fig8.pdf'
fig.savefig(latexPath/filename)
fig.savefig('Figures/' + filename)
# %%
