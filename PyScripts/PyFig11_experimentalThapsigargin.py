#%% Libraries
"""
Fig11
Simulation of the Max SOCE (Thapsigargin) experiment
"""
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import pandas as pd
"""
importing changes to the RcParams
"""
from myOptions import *
from PyModels import parJurkatCell, jurkatWTCDICell, jurkatS2KOCDICell, jurkatS1KOCDICell, icsJurkatCell, thapsProtocol
matplotlib.rcParams.update(myRcParams())
#%% Data collection
"""
Data for the MAX SOCE assay is stored in the thapsMean and thapsSEM
Data for the ER refilling assay is stored in the CPAMean and CPAsem 
"""
thapsMean = pd.read_csv('experimentalData/thapsMean.csv', index_col=0)
thapssem = pd.read_csv('experimentalData/thapssem.csv', index_col=0)
# To reference the different experiments
sheetnames = [
    'Time', 
    'WT', 
    'STIM1 KO\#1', 
    'STIM1 KO\#2', 
    'STIM2 KO\#1', 
    'STIM2 KO\#2', 
]
# Short titles for the legend
shortnames = [
    'Time', 
    'WT', 
    'S1-KO1', 
    'S1-KO2', 
    'S2-KO1', 
    'S2-KO2', 
]
colors = ['tab:purple', 'C3', 'C1', 'C0', 'C9', 'C4']
"""
Data Simulation - WT, S1 and S2 cells
"""
ics = icsJurkatCell()
Tg = 0.5
tini, tf = [0, 3000]
tini, tthaps, tsoce, tguad, tfinal = [0, 100, 1500, 1700, 2000]

par = parJurkatCell()
model = jurkatWTCDICell
WT = thapsProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tguad, tfinal)

par = parJurkatCell()
model = jurkatS2KOCDICell
S1 = thapsProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tguad, tfinal)

par = parJurkatCell()
model = jurkatS1KOCDICell
S2 = thapsProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tguad, tfinal)
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 4, nrows=7, figure=fig)
fig_width,fig_height = set_figsize(1, (1.9, 1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
Max SOCE assay
Plotting individual time series, using SEM as the errorbar
"""
ax = fig.add_subplot(spec[0:3, 0:5])
for k in [5,4,3,2,1]:
    time = thapsMean['Time']
    mean = thapsMean[sheetnames[k]]
    sem = thapssem[sheetnames[k]]
    ax.errorbar(time, mean, yerr=sem, fmt='', color = colors[k-1], alpha=0.5, markersize = 1,)
    ax.plot(time, mean, marker='o', color = colors[k-1], label = shortnames[k], markersize = 1)
# Tg arrow
ax.annotate('Tg', [1.5, 0.72])
ax.annotate('', [2,0.5], [2, 0.7],  arrowprops=dict(arrowstyle="->"))
# Gd3+ arrow
ax.annotate(r'Gd$^{3+}$', [15, 1.03])
ax.annotate('', [15,0.8], [15, 1],  arrowprops=dict(arrowstyle="->"))
#Black rectangle at the top,goes black->white->black
ax.add_patch(Rectangle((0, 1.4), 3, 0.1, fc='k', ec = 'k'))
ax.annotate('0 mM', [3, 1.2])
ax.add_patch(Rectangle((3, 1.4), 9, 0.1, fc='None', ec = 'k'))
ax.annotate('2 mM', [12, 1.2])
ax.add_patch(Rectangle((12, 1.4), 15, 0.1, fc='k', ec = 'k'))
# Ax settings
ax.set_xlim(0, 20)
ax.set_ylim(0.1, 1.5)
ax.set_yticks([0.2, 1.4])
ax.set_xticks([0,4, 8, 12, 16, 20])
ax.set_xlabel('Time (Min)')
ax.set_ylabel(r'F$_{340}$/F$_{380}$')
ax.set_title('A', loc='left')
ax.legend(loc='best', ncol=1, bbox_to_anchor=(1, 0.6,)) #go fall in a ditch, matplotlib.legend
"""
Numerical simulation
Cytoplasmic time series
"""
ax = fig.add_subplot(spec[3:5, :])
ax.plot(WT[0], WT[1], 'tab:purple', label = 'WT')
ax.plot(S1[0], S1[1], 'tab:blue', label = 'S2-KO')
ax.plot(S2[0], S2[1], 'tab:red', label = 'S1-KO')
# These h' are used to fine tune the position of the black rectangles at the top
h1 = 1.1
h12 = 0.13
h2 = 1.1
# Tg arrow
ax.annotate('Tg', [tthaps-50, 0.35])
ax.annotate('', [tthaps,0.1], [tthaps, 0.3],  arrowprops=dict(arrowstyle="->"))
# Gd3+ arrow
ax.annotate('Gd$^{3+}$', [tguad, 0.9])
ax.annotate('', [tguad, 0.65], [tguad, 0.9],  arrowprops=dict(arrowstyle="->"))
#Black rectangle at the top,goes black->white->black
ax.add_patch(Rectangle((1, h1), tthaps, h12, fc='k', ec = 'k'))
ax.add_patch(Rectangle((tthaps+1, h1), tsoce-tthaps, h12, fc='None', ec = 'k'))
ax.add_patch(Rectangle((tsoce+1, h1), tfinal-tsoce-10, h12, fc='k', ec = 'k'))
# Ax limits 
ax.set_xlim(0, tfinal)
ax.set_ylim(0, 1.23)
ax.set_yticks([0, 1.1])
ax.set_xticks([])
# ax.set_xlabel('Time (seconds)')
ax.set_ylabel(r'$c$', rotation='horizontal', labelpad=0)
ax.set_title('B', loc='left')
ax.set_xticks([0,350, 700, 1050, 1400, 1700, 2000])
# ax.set_xlabel('Time (seconds)', labelpad=10)
# ax.legend(loc='best', ncol=3, bbox_to_anchor=(0.1, 0.8, 0.9, 0.5), mode='expand',prop={'size': 8})
ax.legend(loc='best', ncol=1, bbox_to_anchor=(1, 0.6,))
"""
Numerical simulation
Microdomain time series
"""
ax = fig.add_subplot(spec[5:7, :])
ax.plot(WT[0], WT[3], 'tab:purple', label = 'WT')
ax.plot(S1[0], S1[3], 'tab:blue', label = 'S2-KO')
ax.plot(S2[0], S2[3], 'tab:red', label = 'S1-KO')
# These h' are used to fine tune the position of the black rectangles at the top
h1 = 1100
h12 = 150
h2 = 925
#Black rectangle at the top,goes black->white->black
ax.add_patch(Rectangle((1, h1), tthaps, h12, fc='k', ec = 'k'))
ax.add_patch(Rectangle((tthaps+1, h1), tsoce-tthaps, h12, fc='None', ec = 'k'))
ax.add_patch(Rectangle((tsoce+1, h1), tfinal-tsoce-10, h12, fc='k', ec = 'k'))
# Ki line
ax.annotate(r'$K_i$', [700, 730])
ax.plot([S2[0][0], S2[0][-1]], [900, 900], ls = 'dashed', color = 'C7')
# Ax limits 
ax.set_xlim(0, tfinal)
ax.set_ylim(0, 1250)
ax.set_yticks([0, 1100])
ax.set_xticks([0,350, 700, 1050, 1400, 1700, 2000])
ax.set_xlabel('Time (seconds)', labelpad=0)
ax.set_ylabel(r'$c_m$', rotation='horizontal', labelpad=-5)
ax.set_title('C', loc='left')
ax.legend(loc='best', ncol=1, bbox_to_anchor=(1, 0.6,))

filename = 'Fig11.pdf'
saveFigure(filename, fig)