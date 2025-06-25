#%% Libraries 
"""
Fig01 
Secondary experimental assays
"""
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import pandas as pd
"""
importing changes to the RcParams
"""
from myOptions import *
matplotlib.rcParams.update(myRcParams())
#%% Data import
"""
Data for the MAX SOCE assay is stored in the thapsMean and thapsSEM
Data for the ER refilling assay is stored in the CPAMean and CPAsem 
"""
thapsMean = pd.read_csv('experimentalData/thapsMean.csv', index_col=0)
thapssem = pd.read_csv('experimentalData/thapssem.csv', index_col=0)
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
    'STIM1-KO1', 
    'STIM1-KO2', 
    'STIM2-KO1', 
    'STIM2-KO2', 
]
colors = ['tab:purple', 'C3', 'C1', 'C0', 'C9', 'C4']
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
"""
Now, one may ask, why does am I using an 18x2 grid for my 2x1 plot
The answer is simple: the legend box keep rendering outside of the figure
I know it doesn't when I set the grid this way. 
"""
spec = gridspec.GridSpec(ncols = 18, nrows=2, figure=fig)
subSpec1 = spec[0,0:-1].subgridspec(1,1)
subSpec2 = spec[1,0:-1].subgridspec(1,1)
fig_width,fig_height = set_figsize(1, (1.5, 1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
Max SOCE assay
Plotting individual time series, using SEM as the errorbar
"""
ax = fig.add_subplot(subSpec1[0])
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
# Ionomycin arrow
ax.annotate(r'Ionomycin', [20, 0.2])
ax.annotate('', [18,0.2], [19.8, 0.2],  arrowprops=dict(arrowstyle="->"))
#Black rectangle at the top,goes black->white->black
ax.add_patch(Rectangle((0, 1.4), 3, 0.1, fc='k', ec = 'k'))
ax.annotate('0 mM', [3, 1.25])
ax.add_patch(Rectangle((3, 1.4), 9, 0.1, fc='None', ec = 'k'))
ax.annotate('2 mM', [12, 1.25])
ax.add_patch(Rectangle((12, 1.4), 15, 0.1, fc='k', ec = 'k'))
# Ax settings
ax.set_xlim(0, 20)
ax.set_ylim(0.1, 1.5)
ax.set_yticks([0.2, 1.4])
ax.set_xticks([0,4, 8, 12, 16, 20])
ax.set_xlabel('Time (Min)')
ax.set_ylabel(r'F$_{340}$/F$_{380}$')
ax.set_title('A', loc='left')
ax.legend(loc='center right', ncol=1, bbox_to_anchor=(1.05, 0.2, 0.35, 0.5), mode='expand') #This legend will haunt me in my nightmares

"""
ER refilling assay
Plotting individual time series, using SEM as the errorbar
"""
ax = fig.add_subplot(subSpec2[0])
for k in [5,4,3,2,1]:
    time = CPAMean['Time']
    mean = CPAMean[sheetnames[k]]
    sem = CPAsem[sheetnames[k]]
    ax.errorbar(time, mean, yerr=sem, fmt='o', color = colors[k-1], alpha=0.5, markersize=1)
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
ax.set_title('B', loc='left')
ax.legend(loc='center right', ncol=1, bbox_to_anchor=(1, 0.2, 0.4, 0.5), mode='')#This legend will haunt me in my nightmares 

filename = 'Fig3.pdf'
saveFigure(filename, fig)
# %%
