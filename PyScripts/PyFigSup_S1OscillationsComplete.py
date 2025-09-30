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
from myOptions import auto_directory
sys.path.append(auto_directory)
from auto import *
from auto import run, load, save, merge, relabel, cl, klb, loadbd
"""
This command allows me to export the files directly to the 
sister LaTeX directory. Feel free to comment
"""
from myOptions import latexPath
"""
importing changes to the RcParams
"""
from myOptions import *
from PyModels import parJurkatCell, sampleS1Stable, plotTS
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
#     NMX=80000,
#     NTST=1000,
#     DS=1e-3,
#     DSMAX=3e-2,
#     DSMIN=1e-6,
#     SP=['LP100', 'UZ', 'PD100', 'TR0'],
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
spec = gridspec.GridSpec(ncols = 2, nrows=2, figure=fig)
fig_width,fig_height = set_figsize(1, (2,1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
Bifurcation diagram
"""
ax = fig.add_subplot(spec[0, :])
# Extracting the data from the bifurcation diagram object to plot
eqCurve = eqS1
eqCurve = eqCurve[0::1]
eqStab = eqCurve.stability()
eqStab.insert(0,0)
cyCurve = cycleS1
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
xlab = [0.042, 0.1168, 0.1692, 2]
ylab = [0.046, 0.058, 0.148, 2]
for i in eqCurve('HB').getLabels(): 
    point = eqCurve.getLabel(i)
    if i==eqCurve('HB').getLabels()[0]:
        label = 'HB'
    else:
        label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'o', 
        color =  'tab:blue', 
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
            'tab:blue', 
            ls='solid', 
            zorder=-1,
            # label='Limit cycle'
        )
    else:
        ax.plot( #unstable curve is plotted as dashed line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:blue',
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
        color =  'tab:blue',
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
        color =  'tab:blue',
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
cyCurve = cycleS1b
cyCurve = cyCurve[0:40000] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
for i in range(1,len(cycleStab)):
    if cycleStab[i]< 0:
        ax.plot( #Stable curve is plotted as solid line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:blue', 
            ls='solid', 
            zorder=-1,
            # label='Limit cycle'
        )
    else:
        ax.plot( #unstable curve is plotted as dashed line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:blue',
            ls='dotted',
        )
# Labeling of the UZR points
# j=0
# for i in cyCurve('PD').getLabels(): #Labeling of the PD points
#     point = cyCurve.getLabel(i)
#     if i==cyCurve('PD').getLabels()[0]:
#         label = 'PD'
#     else:
#         label = ''
#     ax.scatter(point['Vplc'],
#         max(point['c']),
#         marker = 'D', 
#         facecolors = 'none',
#         color =  'tab:blue',
#         label = label, 
#         zorder=1,
#         )
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
#         color =  'tab:blue',
#         label = label, 
#         zorder=1,
#         )
# Third limit cycle curve
cyCurve = cycleS1c
cyCurve = cyCurve[0:109500] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
for i in range(1,len(cycleStab)):
    if cycleStab[i]< 0:
        ax.plot( #Stable curve is plotted as solid line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:blue', 
            ls='solid', 
            zorder=-1,
            # label='Limit cycle'
        )
    else:
        ax.plot( #unstable curve is plotted as dashed line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:blue',
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
        color =  'tab:blue',
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
        color =  'tab:blue',
        label = label, 
        zorder=1,
        )
# Ax limits 
ax.add_patch(Rectangle((0.28, 0.3), 0.34-0.28, 0.75-0.3, fc='None', ec = 'k'))
ax.add_patch(Rectangle((0.09, 0.05), 0.165-0.09, 0.5-0.05, fc='None', ec = 'k'))
ax.annotate('B', (0.345, 0.75))
ax.annotate('C', (0.17, 0.5))
xlim = [0,0.4]
ylim = [-0.01, 0.75]
yticks = [0, 0.8]
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks([0, 0.25, 0.5, 0.75])
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
ax.set_xlabel(r'$V_{\mathrm{PLC}}$',)
ax.set_ylabel(r'$c$', rotation='horizontal')
ax.set_title('A', loc='left')
ax.legend()

ax = fig.add_subplot(spec[2])
# Extracting the data from the bifurcation diagram object to plot
eqCurve = eqS1
eqCurve = eqCurve[0::1]
eqStab = eqCurve.stability()
eqStab.insert(0,0)
cyCurve = cycleS1
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
for i in eqCurve('HB').getLabels(): 
    point = eqCurve.getLabel(i)
    if i==eqCurve('HB').getLabels()[0]:
        label = 'HB'
    else:
        label = ''
    ax.scatter(point['Vplc'],
        max(point['c']),
        marker = 'o', 
        color =  'tab:blue', 
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
            'tab:blue', 
            ls='solid', 
            zorder=-1,
            # label='Limit cycle'
        )
    else:
        ax.plot( #unstable curve is plotted as dashed line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:blue',
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
        color =  'tab:blue',
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
        color =  'tab:blue',
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
# Ax limits 
xlim = [0.28,0.34]
ylim = [0.3, 0.75]
yticks = [0, 0.8]
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks([0.3, 0.45, 0.6, 0.75])
ax.set_xticks([0.28, 0.3, 0.32, 0.34])
ax.set_xlabel(r'$V_{\mathrm{PLC}}$',)
ax.set_ylabel(r'$c$', rotation='horizontal')
ax.set_title('B', loc='left')
# ax.legend()

ax = fig.add_subplot(spec[3])
# Extracting the data from the bifurcation diagram object to plot
eqCurve = eqS1
eqCurve = eqCurve[0::1]
eqStab = eqCurve.stability()
eqStab.insert(0,0)
cyCurve = cycleS1
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
        color =  'tab:blue', 
        label = label,
        # facecolors = 'none',
        zorder=1,
            )
# Second limit cycle curve
cyCurve = cycleS1b
cyCurve = cyCurve[0:45000] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
for i in range(1,len(cycleStab)):
    if cycleStab[i]< 0:
        ax.plot( #Stable curve is plotted as solid line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:blue', 
            ls='solid', 
            zorder=-1,
            # label='Limit cycle'
        )
    else:
        ax.plot( #unstable curve is plotted as dashed line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:blue',
            ls='dotted',
        )
# Labeling of the UZR points
# j=0
# for i in cyCurve('PD').getLabels(): #Labeling of the PD points
#     point = cyCurve.getLabel(i)
#     if i==cyCurve('PD').getLabels()[0]:
#         label = 'PD'
#     else:
#         label = ''
#     ax.scatter(point['Vplc'],
#         max(point['c']),
#         marker = 'D', 
#         facecolors = 'none',
#         color =  'tab:blue',
#         label = label, 
#         zorder=1,
#         )
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
#         color =  'tab:blue',
#         label = label, 
#         zorder=1,
#         )
# Third limit cycle curve
cyCurve = cycleS1c
cyCurve = cyCurve[0:109500] #don't worry about it
cycleStab = cyCurve.stability()
cycleStab.insert(0,0)
for i in range(1,len(cycleStab)):
    if cycleStab[i]< 0:
        ax.plot( #Stable curve is plotted as solid line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:blue', 
            ls='solid', 
            zorder=-1,
            # label='Limit cycle'
        )
    else:
        ax.plot( #unstable curve is plotted as dashed line
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['Vplc'], 
            cyCurve[np.abs(cycleStab[i-1]):np.abs(cycleStab[i])]['MAX c'], 
            'tab:blue',
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
        color =  'tab:blue',
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
        color =  'tab:blue',
        label = label, 
        zorder=1,
        )
# Ax limits 
xlim = [0.09,0.15]
ylim = [0.05, 0.45]
yticks = [0, 0.8]
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks([0.05, 0.2, 0.35, 0.5])
ax.set_xticks([0.09, 0.115, 0.14, 0.165])
ax.set_xlabel(r'$V_{\mathrm{PLC}}$',)
# ax.set_ylabel(r'$c$', rotation='horizontal')
ax.set_title('C', loc='left')
# ax.legend()

filename = 'Fig_S1OscillationsComplete.pdf'
fig.savefig(latexPath/filename)
saveFigure(filename, fig)# %%
