#%% Libraries
"""
Fig12
Simulation of the ER refilling (CPA) experiment
"""
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import pandas as pd
from numpy import exp
"""
importing changes to the RcParams
"""
from myOptions import *
from PyModels import parJurkatCell, jurkatWTCDICell, jurkatS2KOCDICell, jurkatS1KOCDICell, icsJurkatCell, cpaProtocol, jurkatS2KODynamicCell, jurkatWTDynamicCell, icsJurkatCellDynamicK1
matplotlib.rcParams.update(myRcParams())
#%% Data collection
"""
Data for the MAX SOCE assay is stored in the thapsMean and thapsSEM
Data for the ER refilling assay is stored in the CPAMean and CPAsem 
"""
CPAMean = pd.read_csv('experimentalData/CPAMean.csv', index_col=0)
CPAsem = pd.read_csv('experimentalData/CPAsem.csv', index_col=0)
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
Data Simulation - WT, S1 and S2 normal models
"""
ics = icsJurkatCell()
Tg = 1
tini, tf = [0, 3000]
tini, tthaps, tsoce,tfinal = [0, 200, 4000, 450]
VsDebuff = 0.5

par = parJurkatCell()
model = jurkatWTCDICell
WT = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)

par = parJurkatCell()
model = jurkatS2KOCDICell
S1 = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)

par = parJurkatCell()
model = jurkatS1KOCDICell
S2 = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)
"""
Data collection - model Dynamic K1
"""
ics = icsJurkatCellDynamicK1()
Tg = 1
tini, tf = [0, 3000]
tini, tthaps, tsoce,tfinal = [0, 200, 4000, 450]

par = parJurkatCell()
model = jurkatWTDynamicCell
par['a0']=0.025
WTd = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)

par = parJurkatCell()
model = jurkatS2KODynamicCell
par['a0']=0.025
S1d = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)

ics = icsJurkatCell()
par = parJurkatCell()
model = jurkatS1KOCDICell
par['a0']=0.01
S2d = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)
"""
Additional K_infty and tau_k curves as functions of ce
"""
K1 = par['K1']
Vk=par['Vk']
kn=par['kn']
taukn=par['taun']
Tk=par['Tkmax']
taun=par['taun']
nk=par['nk']
kRest=par['k1Rest']
a0=0.025
w1=par['w1']
n=par['n']
c_e = np.linspace(0, 900, 200)
kCurve = (Vk*kn**nk)/(kn**nk + c_e**nk) + kRest
TauK=Tk*(c_e**nk / (taukn**nk + c_e**nk))
TauExp=Tk*(S1d[2]**nk / (taukn**nk + S1d[2]**nk))
S1low = a0 + w1/(1+exp(n*(c_e-K1)))
S1high = a0 + w1/(1+exp(n*(c_e-560)))
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 2, nrows=4, figure=fig)
subSpec = spec[0:2,:].subgridspec(1,1)
fig_width,fig_height = set_figsize(1, (1.9,1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
ER Refilling assay
Plotting individual time series, using SEM as the errorbar
"""
ax = fig.add_subplot(subSpec[0, :])
for k in [5,4,3,2,1]:
    time = CPAMean['Time']
    mean = CPAMean[sheetnames[k]]
    sem = CPAsem[sheetnames[k]]
    # ax.errorbar(time, mean, yerr=sem, fmt='o', color = colors[k-1], alpha=0.5, markersize=1)
    ax.plot(time, mean, marker='o', color = colors[k-1], label = shortnames[k], markersize=4)
# CPA arrow
ax.annotate('CPA', [2.7, 560])
ax.annotate('', [3,450], [3, 550],  arrowprops=dict(arrowstyle="->"))
# Washes arrow
ax.annotate('Washes', [9.7, 560])
ax.annotate('', [10,450], [10, 550],  arrowprops=dict(arrowstyle="->"))
ax.annotate('', [11,450], [11, 550],  arrowprops=dict(arrowstyle="->"))
#Black rectangle at the top,goes black->white->black
ax.add_patch(Rectangle((0, 650), 3, 50, fc='k', ec = 'k'))
ax.annotate('0 mM', [3, 600])
ax.add_patch(Rectangle((3, 650), 9, 50, fc='None', ec = 'k'))
ax.annotate('2 mM', [12, 600])
ax.add_patch(Rectangle((12, 650), 15, 50, fc='k', ec = 'k')) 
# Ax settings
ax.set_xlim(0, 15)
ax.set_ylim(100, 700)
ax.set_yticks([ 200, 400, 600])
ax.set_xticks([0,3,6,9,12,15])
ax.set_xlabel('Time (Min)')
ax.set_ylabel('RFU')
ax.set_title('A', loc='left')
ax.legend(loc='best', ncol=1, bbox_to_anchor=(1, 0.2, 0.25, 0.5), mode='expand') #This legend will haunt me in my nightmares 
"""
Numerical simulation
ER time series
"""
ax = fig.add_subplot(spec[2, 0])
ax.plot(WT[0], WT[2], 'tab:purple', label = 'WT')
ax.plot(S1[0], S1[2], 'tab:blue', label = 'S2-KO')
ax.plot(S2[0], S2[2], 'tab:red', label = 'S1-KO')
# These h' are used to fine tune the position of the black rectangles at the top
h1 = 1200
h12 = 100
h2 = 950
h22 = 820
h3= 1100
# CPA arrow
ax.annotate('CPA', [tthaps, h2])
ax.annotate('', [tthaps,h22], [tthaps, h2],  arrowprops=dict(arrowstyle="->"))
# Washes arrow
ax.annotate('Wash', [tsoce, h2])
ax.annotate('', [tsoce,h22], [tsoce, h2],  arrowprops=dict(arrowstyle="->"))
#Black rectangle at the top,goes black->white->black
ax.add_patch(Rectangle((1, h1), tthaps, h12, fc='k', ec = 'k'))
ax.add_patch(Rectangle((tthaps+1, h1), tsoce-tthaps, h12, fc='None', ec = 'k'))
ax.add_patch(Rectangle((tsoce+1, h1), tthaps+tfinal, h12, fc='k', ec = 'k'))
# Ax limits 
xlim = [0,tfinal+tsoce]
ylim = [0, 1300]
yticks = [0, 1200]
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks(yticks)
ax.set_xticks(xlim)
ax.set_xlabel('Time (seconds)', labelpad=-10)
ax.set_ylabel(r'$c_e$', rotation='horizontal', labelpad=-10)
ax.set_title('B', loc='left')
ax.legend(loc='best', ncol=3, bbox_to_anchor=(0.1, 0.7, 1.1, 0.5), mode='expand',prop={'size': 8})
"""
Numerical simulation
ER time series - Dynamic K1 model
"""
ax = fig.add_subplot(spec[2,1])
ax.plot(WTd[0], WTd[2], 'tab:purple', label = 'WT')
ax.plot(S1d[0], S1d[2], 'tab:blue', label = 'S1')
ax.plot(S2d[0], S2d[2], 'tab:red', label = 'S2')
# Box 
h1 = 1200
h12 = 100
h2 = 950
h22 = 820
h3= 1100
ax.annotate('CPA', [tthaps, h2])
ax.annotate('', [tthaps,h22], [tthaps, h2],  arrowprops=dict(arrowstyle="->"))
ax.annotate('Wash', [tsoce, h2])
ax.annotate('', [tsoce,h22], [tsoce, h2],  arrowprops=dict(arrowstyle="->"))
ax.add_patch(Rectangle((1, h1), tthaps, h12, fc='k', ec = 'k'))
# ax.annotate('0 mM', [tthaps, h3])
ax.add_patch(Rectangle((tthaps+1, h1), tsoce-tthaps, h12, fc='None', ec = 'k'))
# ax.annotate('2 mM', [tsoce, h3])
ax.add_patch(Rectangle((tsoce+1, h1), tthaps+tfinal, h12, fc='k', ec = 'k'))
# Ax limits 
xlim = [0,tfinal+tsoce]
ylim = [0, 1300]
yticks = [0, 1200]
#Ax limits
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks(yticks)
ax.set_xticks(xlim)
ax.set_xlabel('Time (seconds)', labelpad=-10)
ax.set_title('C', loc='left')
"""
Numerical simulation
K1 vs ER Ca2+ phase plane
"""
ax = fig.add_subplot(spec[3, 0])
ax.plot(S1d[2], S1d[8], color = 'C0', lw=2, label='S1-KO')
# ax.plot(S1d[2], TauExp, 'C3', label = r'$\tau_k(c_e)$')
ax.plot(c_e, kCurve, 'C1', ls =  'dashed', lw=1, label = r'$K_\infty$')
ax.plot(c_e, TauK, 'C2', ls =  'dashed', lw=1, label = r'$\tau_k$')
xlim = [0,800]
ylim = [0, 1000]
#Ax limits
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks([0, 900])
ax.set_xticks(xlim)
ax.set_xlabel(r'$c_e$', labelpad=0)
ax.set_ylabel(r'$K_1$', labelpad=-10, rotation=0 )
ax.set_title('D', loc='left')
ax.legend(loc='best', ncol=3, bbox_to_anchor=(0.1, 0.7, 1, 0.5), mode='expand',prop={'size': 8})
"""
Numerical simulation
CRAC current vs ER Ca2+
"""
ax = fig.add_subplot(spec[3,1])
# STIM2 KO Cell
ax.plot(S1d[2], S1d[6], color = 'C0', lw=2, label='S1-KO')
ax.plot(c_e, S1low, 'C1', ls =  'dashed', lw=1, label = 'Regular \n' + r'$S_1^*$')
ax.plot(c_e, S1high, 'C2', ls =  'dashed', lw=1, label = 'Low ER \n' + r'$S_1^*$')
ax.annotate('', [500, 0.5], [220, 0.5],  arrowprops=dict(arrowstyle="->"))
xlim = [0,800]
ylim = [0, 1.15]
#Ax limits
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks([0, 1.1])
ax.set_xticks(xlim)
ax.set_xlabel(r'$c_e$', labelpad=0)
ax.set_ylabel(r'CRAC current', labelpad=-10)
ax.set_title('E', loc='left')
ax.legend(loc='best', ncol=3, bbox_to_anchor=(0.1, 0.7, 1, 0.5), mode='expand',prop={'size': 6})

filename = 'Fig12.pdf'
saveFigure(filename, fig)