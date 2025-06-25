#%% Libraries
"""
Fig06
Bifurcation diagram and sample time series of the WT model
"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from auto import run, load, save, merge, relabel, cl, klb, loadbd
"""
importing changes to the RcParams
"""
from myOptions import *
from PyModels import parJurkatCell, sampleWTStable, plotTS, saveFigure
matplotlib.rcParams.update(myRcParams())

#%% Data Collection
"""
Bifurcation diagram
Only run this is section if the .WT files have not been created yet
"""
# file = "AUTOJurkatWTCell"
# model = load(file) 
# Vplc = [0.1, 0.2, 0.3]
# eqWT = run(
#     model, 
#     IPS=1, 
#     ICP=['Vplc'], 
#     NMX=20000,
#     DS=1e-2,
#     DSMAX=5e-2,
#     UZSTOP={'Vplc': 2}, 
#     UZR={'Vplc': Vplc}
# )
# Vplc = [0.1, 0.2, 0.3]
# cycleWT = run(
#     eqWT('HB'), 
#     IPS=2,
#     # ISP=2, 
#     ICP=['Vplc', 11, 'MIN c'], 
#     NMX=5000,
#     NTST=500,
#     DS=1e-2,
#     DSMAX=1e-0,
#     SP=['LP2', 'UZ', 'PD2',],
#     UZSTOP={'Vplc': 2}, 
#     UZR={'Vplc': Vplc}
# )
# save(eqWT+cycleWT, 'WT')
# cl()
#%% Data import
"""
Loading WT AUTO file
"""
eqWT = loadbd('WT')[0]
cycleWT = loadbd('WT')[1]
#%% Data import - time series
"""
Data Simulation - WT cells
"""
Vplc = [0.1, 0.2, 0.3]
par = parJurkatCell()
tf = 250
WT = []
for V in Vplc:
    par['Vplc']=0
    sample = sampleWTStable(par, V, tf)
    WT.append(sample)
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 3, nrows=5, figure=fig)
fig_width,fig_height = set_figsize(1, (1.9,1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
Bifurcation diagram
"""
ax = fig.add_subplot(spec[0:2,:])
# Extracting the data from the bifurcation diagram object to plot
eqCurve = eqWT
eqCurve = eqCurve[0::1]
eqStab = eqCurve.stability()
eqStab.insert(0,0)
cyCurve = cycleWT
cyCurve = cyCurve[0::2]
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
    else:
        ax.plot( #unstable curve is plotted as dashed line
            eqCurve[max(np.abs(eqStab[i-1])-1,0):np.abs(eqStab[i])]['Vplc'], 
            eqCurve[max(np.abs(eqStab[i-1])-1,0):np.abs(eqStab[i])]['c'], 
            'k',
            ls='dashed',
            alpha=1, 
            zorder=-1,
        )
#Limit cycle curve
# This curve is completely stable, so no need to do the solid-dashed separation
ax.plot(cyCurve['Vplc'], cyCurve['MAX c'], 'tab:purple')
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
        color =  'tab:purple', 
        label = label,
        zorder=1,
            )
# Labeling of the UZR points
j=0
for i in cyCurve('UZ').getLabels(): 
    point = cyCurve.getLabel(i)
    if point['Vplc'] < 0.5:
        ax.scatter(
            point['Vplc'],
            max(point['c']),
            marker = 'o', 
            color =  'C7' ,
            s=10,
            zorder=3, #High zorder to ensure the point shows above the line
            )
        ax.annotate(
            alph[j],
            [point['Vplc'],max(point['c'])-0.07],
        )
        j=j+1
# Ax limits 
xlim = [0,0.4]
ylim = [0, 0.5]
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks(ylim)
ax.set_xticks(xlim)
ax.set_xlabel(r'$V_{\mathrm{PLC}}$')
ax.set_ylabel(r'$c$')
ax.set_title('A', loc='left')
# ax.legend(ncol=4)
ax.legend( loc='best', ncol=4, bbox_to_anchor=(0, 0.8, 0.2, 0.2), mode='expand')

"""
Time series
"""
# Ax limits (common)
xlim = [0, 250]
xticks = []
yclim = [0, 0.5]
ycelim = [800, 860]
ycmlim = [-10, 750]
ycmticks = [0, 750]
yticks = []
panels = [r'$c$', r'$c_e$', r'$c_m$']
color = 'tab:purple'
xlabel = r'Time (seconds)'
"""
Long function ahead
I need to plot 12 time series, 4 for c, ce and cm 
This is done in this big, ugly loop
Someone smarter than me can do it more elegantly I am sure 
"""
for col, sample in enumerate(WT):
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

filename = 'Fig7.pdf'
saveFigure(filename, fig)
# %%
