#%% Libraries
"""
Fig06
Bifurcation diagram and sample time series of the WT model
"""
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy import exp
"""
importing changes to the RcParams
"""
from STIMKO_Options import *
from STIMKO_Models import parJurkatCell
matplotlib.rcParams.update(myRcParams())
#%% STIM expressions
"""
Equations for the STIM binding model 
Details of the derivation of these expressions are found in the PySym_STIMKOmodel and PySym_STIMWTmodel
"""

def phi1(K_1, K_2, K_12, K_21, c_e, n):
    return exp(n*(-K_1 + c_e))

def phi2(K_1, K_2, K_12, K_21, c_e, n):
    return exp(n*(-K_2 + c_e))

def phi3(K_1, K_2, K_12, K_21, c_e, n):
    return exp(n*(-K_1 + K_21 + c_e))

def phi4(K_1, K_2, K_12, K_21, c_e, n):
    return exp(K_12*n)

def S1a(K_1, K_2, K_12, K_21, c_e, n):
    phi_1=phi1(K_1, K_2, K_12, K_21, c_e, n)
    phi_2=phi2(K_1, K_2, K_12, K_21, c_e, n)
    phi_3=phi3(K_1, K_2, K_12, K_21, c_e, n)
    phi_4=phi4(K_1, K_2, K_12, K_21, c_e, n)
    return (phi_2 + 1)*(np.sqrt((phi_1 + 1)*(phi_2 + 1)*(phi_1*phi_2 + phi_1 + phi_2 + 4*phi_3 + 4*phi_4 + 1))/(2*(phi_2*phi_3 + phi_2*phi_4 + phi_3 + phi_4)) - (phi_1 + 1)/(2*(phi_3 + phi_4)))/(phi_1 + 1)

def S2a(K_1, K_2, K_12, K_21, c_e, n):
    phi_1=phi1(K_1, K_2, K_12, K_21, c_e, n)
    phi_2=phi2(K_1, K_2, K_12, K_21, c_e, n)
    phi_3=phi3(K_1, K_2, K_12, K_21, c_e, n)
    phi_4=phi4(K_1, K_2, K_12, K_21, c_e, n)
    return np.sqrt((phi_1 + 1)*(phi_2 + 1)*(phi_1*phi_2 + phi_1 + phi_2 + 4*phi_3 + 4*phi_4 + 1))/(2*(phi_2*phi_3 + phi_2*phi_4 + phi_3 + phi_4)) - (phi_1 + 1)/(2*(phi_3 + phi_4))

def S21(K_1, K_2, K_12, K_21, c_e, n):
    phi_1=phi1(K_1, K_2, K_12, K_21, c_e, n)
    phi_2=phi2(K_1, K_2, K_12, K_21, c_e, n)
    phi_3=phi3(K_1, K_2, K_12, K_21, c_e, n)
    phi_4=phi4(K_1, K_2, K_12, K_21, c_e, n)
    S1A=S1a(K_1, K_2, K_12, K_21, c_e, n)
    S2A=S2a(K_1, K_2, K_12, K_21, c_e, n)
    return S1A*S2A*phi_3

def S12(K_1, K_2, K_12, K_21, c_e, n):
    phi_1=phi1(K_1, K_2, K_12, K_21, c_e, n)
    phi_2=phi2(K_1, K_2, K_12, K_21, c_e, n)
    phi_3=phi3(K_1, K_2, K_12, K_21, c_e, n)
    phi_4=phi4(K_1, K_2, K_12, K_21, c_e, n)
    S1A=S1a(K_1, K_2, K_12, K_21, c_e, n)
    S2A=S2a(K_1, K_2, K_12, K_21, c_e, n)
    return S1A*S2A*phi_4

#%% Data collection - CRAC expressions
"""
Parameter values
"""
par = parJurkatCell()
n = par['n']
n_ce = 200
c_e = np.linspace(0, 1500, n_ce)
K_2 = par['K2']
K_1 = par['K1']
K_12 = par['K12']
K_21 = par['K21']
w1 = par['w1']
w2 = par['w2']
w12 = par['w12']
w21 = par['w21']
"""
Expression values - WT model
"""
s1a = S1a(K_1, K_2, K_12, K_21, c_e, n)
s2a = S2a(K_1, K_2, K_12, K_21, c_e, n)
s12 = S12(K_1, K_2, K_12, K_21, c_e, n)
s21 = S21(K_1, K_2, K_12, K_21, c_e, n)
Jcrac = w1*s1a + w2*s2a + w21*s12 + w12*s21

"""
Expression values - KO model
"""
S1 = w1/(1+exp(n*(c_e-K_1)))
S2 = w2/(1+exp(n*(c_e-K_2)))

#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 2, nrows=2, figure=fig)
fig_width,fig_height = set_figsize(1, (1.4,1), export=True)
fig.set_size_inches([fig_width,fig_height])
"""
s_infty functions 
WT, S1 (S2KO), S2 (S1KO)
"""
ax = fig.add_subplot(spec[0,:])
ax.plot(c_e, Jcrac, 'tab:purple', ls =  'solid', lw=2, label = r'WT')
ax.plot(c_e, S1, 'tab:blue', ls =  'solid', lw=2, label = r'S2-KO')
ax.plot(c_e, S2, 'tab:red', ls =  'solid', lw=2, label = r'S1-KO')
#Ax limits
ax.set_ylim([0,2.5])
ax.set_yticks([0,2])
ax.set_xlim([0,1000])
ax.set_xticks([0,250, 500, 750, 1000])
ax.set_title('A', loc ='left')
ax.set_ylabel('CRAC \n current', rotation=0, labelpad=20)
ax.set_xlabel(r'ER Ca$^{2+}$ concentration ($\mu$M)')
ax.legend(loc='best', ncol=1, bbox_to_anchor=(0.85, 0.8))
"""
s_infty, S12 and S21 in the WT model
"""
ax = fig.add_subplot(spec[1,0])
ax.plot(c_e, s12, 'tab:orange', ls =  'solid', lw=2, label = r'$S_{12}^*$')
ax.plot(c_e, s21, 'tab:green', ls =  'solid', lw=2, label = r'$S_{21}^*$')
ax.plot(c_e, Jcrac, 'tab:purple', ls =  'solid', lw=2, label = r'$s_\infty$')
#Ax limits
ax.set_ylim([0,2.6])
ax.set_yticks([0,2])
ax.set_xlim([0,1000])
ax.set_xticks([0,250, 500, 750, 1000])
ax.set_title('B', loc ='left')
ax.set_ylabel('Activated \n STIM', rotation=0, labelpad=20)
ax.set_xlabel(r'ER Ca$^{2+}$ concentration ($\mu$M)')
ax.legend(loc='best', ncol=2, bbox_to_anchor=(0.1, 0.8, 0.8, 0.2), mode='expand')

"""
S12, S21, S1a and S2a in the WT model
"""
ax = fig.add_subplot(spec[1,1])
ax.plot(c_e, s12, 'tab:orange', ls =  'solid', lw=2, label = r'$S_{12}^*$')
ax.plot(c_e, s21, 'tab:green', ls =  'solid', lw=2, label = r'$S_{21}^*$')
ax.plot(c_e, s1a, 'tab:blue', ls =  'solid', lw=2, label = r'$S_{1}^*$' )
ax.plot(c_e, s2a, 'tab:red', ls =  'solid', lw=2, label = r'$S_{2}^*$')
#Ax limits
ax.set_ylim([-0.01, 1.4])
ax.set_yticks([0,1])
ax.set_xlim([0,1000])
ax.set_xticks([0,250, 500, 750, 1000])
ax.set_title('C', loc ='left')
ax.set_xlabel(r'ER Ca$^{2+}$ concentration ($\mu$M)')
ax.legend(loc='best', ncol=2, bbox_to_anchor=(0.1, 0.9, 0.8, 0.2), mode='expand')
filename = 'Fig10.pdf'
saveFigure(filename, fig)