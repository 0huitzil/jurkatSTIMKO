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
#% Data collection 
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

#% Data collection - Modified STIM weights
"""
Data simulation - new parameter values
"""
w1 = 0.4 # w1 = 1
w2 = 1.5 # w2 = 0.7
w12 = 3  # w12 = 3
w21 = 3  # w21 = 3
def updateWeights(w1, w2, w12, w21, par):
    par['w1'] = w1
    par['w2'] = w2
    par['w12'] = w12
    par['w21'] = w21
    return 0 
"""
Data Simulation - WT, S1 and S2 normal models
"""
ics = icsJurkatCell()
Tg = 1
tini, tf = [0, 3000]
tini, tthaps, tsoce,tfinal = [0, 200, 4000, 450]
VsDebuff = 0.5

par = parJurkatCell()
updateWeights(w1, w2, w12, w21, par)
model = jurkatWTCDICell
WT_2 = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)

par = parJurkatCell()
updateWeights(w1, w2, w12, w21, par)
model = jurkatS2KOCDICell
S1_2 = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)

par = parJurkatCell()
updateWeights(w1, w2, w12, w21, par)
model = jurkatS1KOCDICell
S2_2 = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)
"""
Data collection - model Dynamic K1
"""
ics = icsJurkatCellDynamicK1()
Tg = 1
tini, tf = [0, 3000]
tini, tthaps, tsoce,tfinal = [0, 200, 4000, 450]

par = parJurkatCell()
updateWeights(w1, w2, w12, w21, par)
model = jurkatWTDynamicCell
par['a0']=0.025
WTd_2 = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)

par = parJurkatCell()
updateWeights(w1, w2, w12, w21, par)
model = jurkatS2KODynamicCell
par['a0']=0.025
S1d_2 = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)

ics = icsJurkatCell()
par = parJurkatCell()
updateWeights(w1, w2, w12, w21, par)
model = jurkatS1KOCDICell
par['a0']=0.01
S2d_2 = cpaProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tfinal, VsDebuff)
#% Plotting
"""
Figure settings
"""
fig = plt.figure(tight_layout = True)
spec = gridspec.GridSpec(ncols = 2, nrows=2, figure=fig)
subSpec = spec[0:2,:].subgridspec(1,1)
fig_width,fig_height = set_figsize(1, (1.5,1), export=True) 
fig.set_size_inches([fig_width,fig_height])

"""
Numerical simulation
ER time series
"""
ax = fig.add_subplot(spec[0, 0])
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
ax.set_title('Base model - normal weights', loc='left')

"""
Numerical simulation
ER time series - Dynamic K1 model
"""
ax = fig.add_subplot(spec[0,1])
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
ax.set_title('Dynamic K1 model - normal weights', loc='left')


"""
Numerical simulation
ER time series  - modified weights
"""
ax = fig.add_subplot(spec[1, 0])
ax.plot(WT_2[0], WT_2[2], 'tab:purple', label = 'WT')
ax.plot(S1_2[0], S1_2[2], 'tab:blue', label = 'S2-KO')
ax.plot(S2_2[0], S2_2[2], 'tab:red', label = 'S1-KO')
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
ax.set_title('Base model - modified weights', loc='left')
ax.legend(loc='best', ncol=3, bbox_to_anchor=(0, 0.8, 1.2, 0.5), mode='expand',prop={'size': 8})
"""
Numerical simulation
ER time series - Dynamic K1 model - modified weights
"""
ax = fig.add_subplot(spec[1,1])
ax.plot(WTd_2[0], WTd_2[2], 'tab:purple', label = 'WT')
ax.plot(S1d_2[0], S1d_2[2], 'tab:blue', label = 'S1')
ax.plot(S2d_2[0], S2d_2[2], 'tab:red', label = 'S2')
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
ax.set_title('Dynamic K1 model - modified weights', loc='left')

filename = 'FigSup_STIMWeights.pdf'
saveFigure(filename, fig)
# %%
