#%% Libraries
"""
Fig07
Bifurcation diagram and sample time series of the S2 (S2KO) model
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
from STIMKO_Models import parJurkatCell, sampleS2Stable, plotTS
from STIMKO_AUTO import loadBifDiag
matplotlib.rcParams.update(myRcParams())
#%% Data import
"""
Loading S2 AUTO file. This file is now uploaded to the repository by default due to its size. Generate it locally first with the corresponding script in the AUTO folder. """
bd = loadBifDiag('S2a')
eqS2 = bd[0]
cycleS2 = bd[1]
cycleS2b = bd[2]
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 2, nrows=2, figure=fig)
fig_width,fig_height = set_figsize(1, (2,1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
Bifurcation diagram
"""
ax = fig.add_subplot(spec[0,:])
# Extracting the data from the bifurcation diagram object to plot
eqCurve = eqS2
eqCurve = eqCurve[0::1]
eqStab = eqCurve.stability()
eqStab.insert(0,0)
cyCurve = cycleS2
# cyCurve = cyCurve[0:109500] #don't worry about it
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
# Labeling of the HB points
j=1
xlab = [0.0349, 0.0767, 0.1010, 0.1737]
ylab = [0.059, 0.049, 0.07, 0.173]
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
    ax.annotate(
            'HB' + str(j),
            [xlab[j-1], ylab[j-1]],
        )
    j=j+1
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
# Labeling of the UZR points
j=0
for i in cyCurve('PD').getLabels(): #Labeling of the PD points
    point = cyCurve.getLabel(i)
    # if i==cyCurve('PD').getLabels()[0]:
    #     label = 'PD'
    # else:
    label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'D', 
        facecolors = 'none',
        color =  'tab:red',
        label = label, 
        zorder=1,
        )
# for i in cyCurve('LP').getLabels(): #Labeling of the PD points
#     point = cyCurve.getLabel(i)
#     if i==cyCurve('LP').getLabels()[0]:
#         label = 'SNPO'
#     else:
#         label = ''
#     ax.scatter(point['Vplc'],
#         max(point['c']),
#         marker = 'v', 
#         facecolors = 'none',
#         color =  'tab:red',
#         label = label, 
#         zorder=1,
#         )
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
# Second limit cycle curve
cyCurve = cycleS2b
cyCurve = cyCurve[0:109500] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
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
# Labeling of the UZR points
j=0
for i in cyCurve('PD').getLabels(): #Labeling of the PD points
    point = cyCurve.getLabel(i)
    if i==cyCurve('PD').getLabels()[0]:
        label = 'PD'
    else:
        label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'D', 
        facecolors = 'none',
        color =  'tab:red',
        label = label, 
        zorder=1,
        )
# for i in cyCurve('LP').getLabels(): #Labeling of the PD points
#     point = cyCurve.getLabel(i)
#     if i==cyCurve('LP').getLabels()[0]:
#         label = 'SNPO'
#     else:
#         label = ''
#     ax.scatter(point['Vplc'],
#         max(point['c']),
#         marker = 'v', 
#         facecolors = 'none',
#         color =  'tab:red',
#         label = label, 
#         zorder=1,
#         )
# Ax limits 
ax.add_patch(Rectangle((0.16, 0.3), 0.175-0.16, 0.6-0.3, fc='None', ec = 'k'))
ax.add_patch(Rectangle((0.05, 0.04), 0.06, 0.13-0.04, fc='None', ec = 'k'))
ax.annotate('B', (0.178, 0.6))
ax.annotate('C', (0.1115, 0.13))
xlim = [0,0.2]
ylim = [-0.001, 0.6]
yticks = [0, 0.6]
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks([0, 0.2, 0.4, 0.6])
ax.set_xticks([0,0.05, 0.1, 0.15, 0.2])
ax.set_xlabel(r'$V_{\mathrm{PLC}}$',)
ax.set_ylabel(r'$c$', rotation='horizontal')
ax.set_title('A', loc='left')
ax.legend(loc = 'upper left')

ax = fig.add_subplot(spec[2])
# Extracting the data from the bifurcation diagram object to plot
eqCurve = eqS2
eqCurve = eqCurve[0::1]
eqStab = eqCurve.stability()
eqStab.insert(0,0)
cyCurve = cycleS2
# cyCurve = cyCurve[0:109500] #don't worry about it
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
# Labeling of the UZR points
j=0
for i in cyCurve('PD').getLabels(): #Labeling of the PD points
    point = cyCurve.getLabel(i)
    # if i==cyCurve('PD').getLabels()[0]:
    #     label = 'PD'
    # else:
    label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'D', 
        facecolors = 'none',
        color =  'tab:red',
        label = label, 
        zorder=1,
        )
for i in cyCurve('LP').getLabels(): #Labeling of the PD points
    point = cyCurve.getLabel(i)
    if i==cyCurve('LP').getLabels()[0]:
        label = 'SNPO'
    else:
        label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'v', 
        facecolors = 'none',
        color =  'tab:red',
        label = label, 
        zorder=1,
        )
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
# Second limit cycle curve
cyCurve = cycleS2b
cyCurve = cyCurve[0:109500] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
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
# Labeling of the UZR points
j=0
for i in cyCurve('PD').getLabels(): #Labeling of the PD points
    point = cyCurve.getLabel(i)
    if i==cyCurve('PD').getLabels()[0]:
        label = 'PD'
    else:
        label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'D', 
        facecolors = 'none',
        color =  'tab:red',
        label = label, 
        zorder=1,
        )
# for i in cyCurve('LP').getLabels(): #Labeling of the PD points
#     point = cyCurve.getLabel(i)
#     if i==cyCurve('LP').getLabels()[0]:
#         label = 'SNPO'
#     else:
#         label = ''
#     ax.scatter(point['Vplc'],
#         max(point['c']),
#         marker = 'v', 
#         facecolors = 'none',
#         color =  'tab:red',
#         label = label, 
#         zorder=1,
#         )
# Ax limits 
xlim = [0.16,0.175]
ylim = [0.3, 0.6]
yticks = [0, 0.6]
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks([0.3, 0.4, 0.5, 0.6])
ax.set_xticks([0.16, 0.165, 0.17,0.175])
ax.set_xlabel(r'$V_{\mathrm{PLC}}$',)
ax.set_ylabel(r'$c$', rotation='horizontal')
ax.set_title('B', loc='left')
# ax.legend()


ax = fig.add_subplot(spec[3])
# Extracting the data from the bifurcation diagram object to plot
eqCurve = eqS2
eqCurve = eqCurve[0::1]
eqStab = eqCurve.stability()
eqStab.insert(0,0)
cyCurve = cycleS2
# cyCurve = cyCurve[0:109500] #don't worry about it
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
# Labeling of the UZR points
j=0
for i in cyCurve('PD').getLabels(): #Labeling of the PD points
    point = cyCurve.getLabel(i)
    # if i==cyCurve('PD').getLabels()[0]:
    #     label = 'PD'
    # else:
    label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'D', 
        facecolors = 'none',
        color =  'tab:red',
        label = label, 
        zorder=1,
        )
for i in cyCurve('LP').getLabels(): #Labeling of the PD points
    point = cyCurve.getLabel(i)
    if i==cyCurve('LP').getLabels()[0]:
        label = 'SNPO'
    else:
        label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'v', 
        facecolors = 'none',
        color =  'tab:red',
        label = label, 
        zorder=1,
        )
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
# Second limit cycle curve
cyCurve = cycleS2b
cyCurve = cyCurve[0:109500] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
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
# Labeling of the UZR points
j=0
for i in cyCurve('PD').getLabels(): #Labeling of the PD points
    point = cyCurve.getLabel(i)
    if i==cyCurve('PD').getLabels()[0]:
        label = 'PD'
    else:
        label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'D', 
        facecolors = 'none',
        color =  'tab:red',
        label = label, 
        zorder=1,
        )
# for i in cyCurve('LP').getLabels(): #Labeling of the PD points
#     point = cyCurve.getLabel(i)
#     if i==cyCurve('LP').getLabels()[0]:
#         label = 'SNPO'
#     else:
#         label = ''
#     ax.scatter(point['Vplc'],
#         max(point['c']),
#         marker = 'v', 
#         facecolors = 'none',
#         color =  'tab:red',
#         label = label, 
#         zorder=1,
#         )
# Ax limits 
xlim = [0.05,0.11]
ylim = [0.04, 0.13]
yticks = [0, 0.6]
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks([0.04, 0.07, 0.1, 0.13])
ax.set_xticks([0.05, 0.07, 0.09, 0.11])
ax.set_xlabel(r'$V_{\mathrm{PLC}}$',)
# ax.set_ylabel(r'$c$', rotation='horizontal')
ax.set_title('C', loc='left')
# ax.legend()

filename = 'Fig_S2OscillationsComplete.pdf'
fig.savefig(latexPath/filename)
saveFigure(filename, fig)# %%
