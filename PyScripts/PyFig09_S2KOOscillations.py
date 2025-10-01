#%% Libraries
"""
Fig08
Bifurcation diagram and sample time series of the S2KO model
"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from auto import run, load, save, merge, relabel, cl, klb, loadbd
"""
importing changes to the RcParams
"""
from STIMKO_Options import *
from STIMKO_Models import parJurkatCell, sampleS2KOStable, plotTS
from STIMKO_AUTO import loadBifDiag
matplotlib.rcParams.update(myRcParams())

#%% Data import
"""
Loading S2KO AUTO file. If it does not exist, obtain it with the getS1KOBifDiag function on the STIMKO_AUTO script. 
"""
bd = loadBifDiag('S2KO')
eqS2KO = bd[0]
cycleS2KO = bd[1]
#%% Data import - time series
"""
Data Simulation - S2KO cells
"""
Vplc = [0.06, 0.1, 0.2, 0.35]
par = parJurkatCell()
tf = 250
S2KO = []
for V in Vplc:
    par['Vplc']=0
    sample = sampleS2KOStable(par, V, tf)
    S2KO.append(sample)
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 4, nrows=5, figure=fig)
fig_width,fig_height = set_figsize(1, (1.9,1), export=True) 
fig.set_size_inches([fig_width,fig_height])
"""
Ax settings Vplc Diag 
"""
ax = fig.add_subplot(spec[0:2,:])
# Extracting the data from the bifurcation diagram object to plot
eqCurve = eqS2KO
eqCurve = eqCurve[0::1]
eqStab = eqCurve.stability()
eqStab.insert(0,0)
cyCurve = cycleS2KO
cyCurve = cyCurve[0::1]
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
            color =  'C7',
            s=10,
            zorder=3, #High zorder to ensure the point shows above the line
            )
        if point['Vplc'] < 0.3: #Sneaky avoiding plotting an extra UZR around Vplc=2
            ax.annotate(
                alph[j],
                [point['Vplc'],max(point['c'])-0.1],
            )
        else:
            ax.annotate(
                alph[j],
                [point['Vplc'],max(point['c'])+0.05],
            )
        j=j+1
# Ax limits 
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
ax.legend(ncol=1)

"""
Time series
"""
# Ax limits (common)
xlim = [0, 250]
xticks = []
yclim = [0,0.8]
ycelim = [180, 240]
ycmlim = [-10, 500]
ycmticks = [0, 500]
yticks = []
panels = [r'$c$', r'$c_e$', r'$c_m$']
color = 'tab:blue'
"""
Long function ahead
I need to plot 12 time series, 4 for c, ce and cm 
This is done in this big, ugly loop
Someone smarter than me can do it more elegantly I am sure 
"""
for col, sample in enumerate(S2KO):
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

filename = 'Fig9.pdf'
saveFigure(filename, fig)
# %%
