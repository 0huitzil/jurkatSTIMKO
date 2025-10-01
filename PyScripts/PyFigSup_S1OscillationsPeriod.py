#%% Libraries
"""
Fig07
Bifurcation diagram and sample time series of the S1 (S1KO) model
"""
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from pathlib import Path
parentPath = str(Path(os.getcwd()).parent)
sys.path.append(parentPath)
"""
The path to the AUTO library is setup in myOptions.py
Make sure to set the correct path there first before running these files
"""
from STIMKO_Options import auto_directory
sys.path.append(auto_directory)
from auto import *
from auto import run, load, save, merge, relabel, cl, klb, loadbd
"""
This command allows me to export the files directly to the 
sister LaTeX directory. Feel free to comment
"""
from STIMKO_Options import latexPath
"""
importing changes to the RcParams
"""
from STIMKO_Options import *
from STIMKO_Models import parJurkatCell, sampleS1Stable, plotTS
matplotlib.rcParams.update(myRcParams())

#%% Data Collection
"""
Data Vplc Diag 
Only run this is if the .S1a files have not been created yet
"""
# file = "AUTOJurkatS1Cell"
# model = load(file) 
# eqS1 = run(
#     model, 
#     IPS=1, 
#     ICP=['Vplc'], 
#     NMX=20000,
#     DS=1e-2,
#     DSMAX=5e-2,
#     UZSTOP={'Vplc': 2}
# )
# #%
# # Vplc = [0.06, 0.1, 0.2, 0.35]
# cycleS1 = run(
#     eqS1('HB')[0], 
#     IPS=2,
#     JAC=1,
#     ICP=['Vplc', 11], 
#     NMX=150000,
#     NTST=1000,
#     DS=1e-3,
#     DSMAX=3e-2,
#     DSMIN=1e-6,
#     SP=['LP0', 'UZ', 'PD', 'TR0'],
#     UZSTOP={'Vplc': 1.76}, 
#     # UZR={'Vplc': Vplc}, 
# )
# save(eqS1+cycleS1, 'S1a')
# cycleS1b = run(
#     eqS1('HB')[1], 
#     IPS=2,
#     JAC=1,
#     ICP=['Vplc', 11], 
#     NMX=40000,
#     NTST=1000,
#     DS=1e-3,
#     DSMAX=3e-2,
#     DSMIN=1e-6,
#     SP=['LP0', 'UZ', 'PD0', 'TR0'],
#     UZSTOP={'Vplc': 2}, 
#     # UZR={'Vplc': Vplc}, 
# )
# save(eqS1+cycleS1+cycleS1b, 'S1a')
# cycleS1c = run(
#     eqS1('HB')[2], 
#     IPS=2,
#     JAC=1,
#     ICP=['Vplc', 11], 
#     NMX=40000,
#     NTST=1000,
#     DS=1e-3,
#     DSMAX=3e-2,
#     DSMIN=1e-6,
#     SP=['LP0', 'UZ', 'PD0', 'TR0'],
#     UZSTOP={'Vplc': 2}, 
#     # UZR={'Vplc': Vplc}, 
# )
# save(eqS1+cycleS1+cycleS1b+cycleS1c, 'S1a')
# cl()
#%% Data import
"""
Loading S1 AUTO file
"""
eqS1 = loadbd('S1a')[0]
cycleS1 = loadbd('S1a')[1]
cycleS1b = loadbd('S1b')[0]
cycleS1c = loadbd('S1a')[3]
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 1, nrows=3, figure=fig)
fig_width,fig_height = set_figsize(1, (2,1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
Bifurcation diagram
"""
ax = fig.add_subplot(spec[0])
# Extracting the data from the bifurcation diagram object to plot
cyCurve = cycleS1
# cyCurve = cyCurve[0:109500] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
#Alphabet
alph = [chr(item) for item in range(ord("B"), ord("Z") + 1)]
#Limit cycle curve
for i in range(1,len(cycleStab)):
    if cycleStab[i]< 0:
        ax.plot( #Stable curve is plotted as solid line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['PERIOD'], 
            'tab:blue', 
            ls='solid', 
            zorder=-1,
            # label='Limit cycle'
        )
    else:
        ax.plot( #unstable curve is plotted as dashed line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['PERIOD'], 
            'tab:blue',
            ls='dotted',
        )
point = cyCurve[0]
ax.scatter(point['Vplc'],
        point['PERIOD'],
        marker = 'o', 
        facecolors = 'tab:blue',
        color =  'tab:blue',
        label = 'HB1', 
        zorder=1,
        )
xlim = [0,0.4]
ylim = [525]
yticks = [0, 0.8]
ax.set_xlim(xlim)
# ax.set_ylim(ylim)
ax.set_yticks([0, 175, 350, 525])
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
# ax.set_xlabel(r'$V_{\mathrm{PLC}}$',)
ax.set_ylabel(r'Period (s)', rotation='vertical')
ax.set_title('A', loc='left')
ax.legend()

# Second limit cycle curve
ax = fig.add_subplot(spec[1])
cyCurve = cycleS1b
cyCurve = cyCurve[0:45000] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
for i in range(1,len(cycleStab)):
    if cycleStab[i]< 0:
        ax.plot( #Stable curve is plotted as solid line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['PERIOD'], 
            'tab:blue', 
            ls='solid', 
            zorder=-1,
            # label='Limit cycle'
        )
    else:
        ax.plot( #unstable curve is plotted as dashed line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['PERIOD'], 
            'tab:blue',
            ls='dotted',
        )
# Labeling of the UZR points
point = cyCurve[0]
ax.scatter(point['Vplc'],
        point['PERIOD'],
        marker = 'o', 
        facecolors = 'tab:blue',
        color =  'tab:blue',
        label = 'HB2', 
        zorder=1,
        )
xlim = [0,0.4]
ylim = [0, 75]
yticks = [0, 0.8]
ax.set_xlim(xlim)
# ax.set_ylim(ylim)
ax.set_yticks([0, 25, 50, 75])
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
# ax.set_xlabel(r'$V_{\mathrm{PLC}}$',)
ax.set_ylabel(r'Period (s)', rotation='vertical')
ax.set_title('B', loc='left')
ax.legend()

# Third limit cycle curve
ax = fig.add_subplot(spec[2])
cyCurve = cycleS1c
cyCurve = cyCurve[0:109500] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
for i in range(1,len(cycleStab)):
    if cycleStab[i]< 0:
        ax.plot( #Stable curve is plotted as solid line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['PERIOD'], 
            'tab:blue', 
            ls='solid', 
            zorder=-1,
            # label='Limit cycle'
        )
    else:
        ax.plot( #unstable curve is plotted as dashed line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['PERIOD'], 
            'tab:blue',
            ls='dotted',
        )
# Labeling of the UZR points
point = cyCurve[0]
ax.scatter(point['Vplc'],
        point['PERIOD'],
        marker = 'o', 
        facecolors = 'tab:blue',
        color =  'tab:blue',
        label = 'HB3', 
        zorder=1,
        )
xlim = [0,0.4]
ylim = [-20, 1200]
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks([0, 400, 800, 1200])
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
ax.set_xlabel(r'$V_{\mathrm{PLC}}$',)
ax.set_ylabel(r'Period (s)', rotation='vertical')
ax.set_title('C', loc='left')
ax.legend()

filename = 'Fig_S1OscillationsPeriod.pdf'
saveFigure(filename, fig)
# %%
