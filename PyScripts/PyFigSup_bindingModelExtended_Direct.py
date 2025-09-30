#%% Libraries
"""
Fig10
Bifurcation diagram and sample time series of the WT model
"""
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy import exp, sqrt
from scipy.optimize import newton, brentq
"""
importing changes to the RcParams
"""
from myOptions import *
from PyModels import parExtendedJurkatCell, parJurkatCell
# from PySym_STIMWTModelExtended import S1aPoly
matplotlib.rcParams.update(myRcParams())
#%% STIM expressions
"""
Equations for the STIM binding model 
Details of the derivation of these expressions are found in the PySym_STIMKOmodelExtended and PySym_STIMWTmodelExtended
"""

def phi1(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n):
    return exp(K_1*n)*exp(-c_e*n)

def phi2(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n):
    return exp(K_2*n)*exp(-c_e*n)

def phi3(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n):
    return exp(-K_1*n)*exp(K_21*n)*exp(c_e*n)

def phi4(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n):
    return exp(K_12*n)

def phi5(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n):
    return exp(K_22*n)

def phi6(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n):
    return exp(K_11*n)

def phi7(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n):
    return exp(K_21*n)

def S1star(S1, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    # F holds the implicit equation on (c_e, S1) used to calculate S1
    F = S1**2*phi_2**2*(phi_1*phi_4 + phi_7)**2*(2*S1**2*phi_6 + S1*phi_1 + S1 - 1) + S1**2*phi_2**2*(phi_1*phi_4 + phi_7)**2 + S1*phi_2**2*(phi_1*phi_4 + phi_7)*(2*S1**2*phi_6 + S1*phi_1 + S1 - 1) + S1*phi_2*(phi_1*phi_4 + phi_7)*(2*S1**2*phi_6 + S1*phi_1 + S1 - 1) - 2*phi_5*(2*S1**2*phi_6 + S1*phi_1 + S1 - 1)**2
    return F 

def S1starPrime(S1, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    # F holds the derivative of S1aS"tar
    Fprime = S1**2*phi_2**2*(phi_1*phi_4 + phi_7)**2*(4*S1*phi_6 + phi_1 + 1) + 2*S1*phi_2**2*(phi_1*phi_4 + phi_7)**2*(2*S1**2*phi_6 + S1*phi_1 + S1 - 1) + 2*S1*phi_2**2*(phi_1*phi_4 + phi_7)**2 + S1*phi_2**2*(phi_1*phi_4 + phi_7)*(4*S1*phi_6 + phi_1 + 1) + S1*phi_2*(phi_1*phi_4 + phi_7)*(4*S1*phi_6 + phi_1 + 1) + phi_2**2*(phi_1*phi_4 + phi_7)*(2*S1**2*phi_6 + S1*phi_1 + S1 - 1) + phi_2*(phi_1*phi_4 + phi_7)*(2*S1**2*phi_6 + S1*phi_1 + S1 - 1) - 2*phi_5*(8*S1*phi_6 + 2*phi_1 + 2)*(2*S1**2*phi_6 + S1*phi_1 + S1 - 1)
    return Fprime

def S1starPrime2(S1, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    # F holds the second derivative of S1aStar
    Fprime = 4*S1**2*phi_2**2*phi_6*(phi_1*phi_4 + phi_7)**2 + 4*S1*phi_2**2*phi_6*(phi_1*phi_4 + phi_7) + 4*S1*phi_2**2*(phi_1*phi_4 + phi_7)**2*(4*S1*phi_6 + phi_1 + 1) + 4*S1*phi_2*phi_6*(phi_1*phi_4 + phi_7) + 2*phi_2**2*(phi_1*phi_4 + phi_7)**2*(2*S1**2*phi_6 + S1*phi_1 + S1 - 1) + 2*phi_2**2*(phi_1*phi_4 + phi_7)**2 + 2*phi_2**2*(phi_1*phi_4 + phi_7)*(4*S1*phi_6 + phi_1 + 1) + 2*phi_2*(phi_1*phi_4 + phi_7)*(4*S1*phi_6 + phi_1 + 1) - 16*phi_5*phi_6*(2*S1**2*phi_6 + S1*phi_1 + S1 - 1) - 2*phi_5*(4*S1*phi_6 + phi_1 + 1)*(8*S1*phi_6 + 2*phi_1 + 2)
    return Fprime

def S1(phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    # We used the Scipy optimization package to numerically calculate 
    # the root of S1astar
    S1a = newton(S1star, x0=0.001, args=(phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n) , fprime=S1starPrime, fprime2=S1starPrime2, maxiter=100, tol=1e-5, rtol=1e-5)
    # S1a = brentq(S1star, a=0.1, b=20, args=(phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n), xtol=1e-30, maxiter=300)
    # if S1a < 1e-10:
    #     S1a=0
    return S1a


S1 = np.vectorize(S1)

def S2(S1, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    # return -(2*s1a**2*phi_6 + s1a*phi_1 + s1a - 1)/(s1a*(phi_1*phi_7 + phi_4))
    return -(2*S1**2*phi_6 + S1*phi_1 + S1 - 1)/(S1*phi_2*(phi_1*phi_4 + phi_7))

def S21(s1, s2, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    return s1*s2*phi_2*phi_7

def S12(s1, s2, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    return s1*s2*phi_1*phi_2*phi_4

def S11(s1, s2, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    return s1**2*phi_6

def S22(s1, s2, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    return  s2**2*phi_5

def S1a(s1, s2, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    return s1*phi_1

def S2a(s1, s2, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n):
    return  s2*phi_2

#%% Data collection - CRAC expressions
"""
Parameter values
"""
par = parExtendedJurkatCell()
n = par['n']
n_ce = 1000
c_e = np.linspace(1, 1000, n_ce)
# c_e = 150
K_2 = par['K2']
K_1 = par['K1']
K_22 = par['K22']
K_11 = par['K11']
K_12 = par['K12']
K_21 = par['K21']
w1 = par['w1']
w2 = par['w2']
w12 = par['w12']    
w21 = par['w21']
w11 = par['w11']
w22 = par['w22']
#%% Data collection - WT model
"""
Expression values - WT model
"""
# Phi expressions 
phi_1=phi1(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n)
phi_2=phi2(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n)
phi_3=phi3(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n)
phi_4=phi4(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n)
phi_5=phi5(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n)
phi_6=phi6(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n)
phi_7=phi7(K_1, K_2, K_12, K_21, K_11, K_22, c_e, n)
# STIM expressions
s1 = S1(phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n)
s2 = S2(s1, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n)
s1a = S1a(s1, s2, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n)
s2a = S2a(s1, s2, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n)
s12 = S12(s1a, s2a, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n)
s21 = S21(s1a, s2a, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n)
s22 = S22(s1a, s2a, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n)
s11 = S11(s1a, s2a, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, c_e, n)
Jcrac = w1*s1a + w2*s2a + w21*s12 + w12*s21 + w11*s11 + w22*s22
#%% Data collection - STIM1-KO model
"""
Expression values STIM1-KO model
"""
S2aKO = -(phi_2 + 1)/(4*phi_5) + sqrt(phi_2**2 + 2*phi_2 + 8*phi_5 + 1)/(4*phi_5)
S22KO = S2aKO**2*phi_5
S1KO_CRAC = w2*S2aKO + w22*S22KO
#%% Data collection - STIM2-KO model
"""
Expression values STIM2-KO model
"""
S1aKO = -(phi_1 + 1)/(4*phi_6) + sqrt(phi_1**2 + 2*phi_1 + 8*phi_6 + 1)/(4*phi_6)
S11KO = S1aKO**2*phi_6
S2KO_CRAC = w1*S1aKO + w11*S11KO

#%% Data collection - Base model
"""
Expression values (base model)
"""
parBase = parJurkatCell()
n1 = parBase['n']
K_1 = parBase['K1']
K_2 = parBase['K2']
w1 = parBase['w1']
w2 = parBase['w2']
S2KOBase = w1/(1+exp(n1*(c_e-K_1)))
S1KOBase = w2/(1+exp(n1*(c_e-K_2)))
WTBase = 3/(1+exp(n1*(c_e-800)))
#%% Plotting
"""
Figure settings
"""
fig = plt.figure(constrained_layout = True)
spec = gridspec.GridSpec(ncols = 1, nrows=3, figure=fig)
fig_width,fig_height = set_figsize(1, (1.85,1), export=True)
fig.set_size_inches([fig_width,fig_height])

"""
S12, S21, S1a and S2a in the WT model
"""
ax = fig.add_subplot(spec[0,0])
ax.plot(c_e, s1, 'tab:blue', ls =  'dashed', lw=2, label = r'$S_{1}$' )
ax.plot(c_e, s2, 'tab:red', ls =  'dashed', lw=2, label = r'$S_{2}$')
# ax.plot(c_e, s1a, 'tab:blue', ls =  'solid', lw=2, label = r'$S_{1}^*$' )
# ax.plot(c_e, s2a, 'tab:red', ls =  'solid', lw=2, label = r'$S_{2}^*$')
ax.plot(c_e, s12, 'tab:orange', ls =  'solid', lw=2, label = r'$S_{12}^*$')
ax.plot(c_e, s21, 'tab:green', ls =  'solid', lw=2, label = r'$S_{21}^*$')
# ax.plot(c_e, s11, 'tab:pink', ls =  'solid', lw=2, label = r'$S_{11}^*$' )
# ax.plot(c_e, s22, 'tab:cyan', ls =  'solid', lw=2, label = r'$S_{22}^*$')
ax.plot(c_e, Jcrac, 'tab:purple', ls =  'solid', lw=2, label = r'WT CRAC')
ax.plot(c_e, WTBase, 'tab:purple', ls =  'dashed', lw=2, label = r'WT CRAC (base)')
#Ax limits
ax.set_ylim([-0.01, 4])
ax.set_yticks([0, 1, 3])
ax.set_xlim([0,1000])
# ax.set_xticks([0,250, 400,])
ax.set_title('A', loc ='left')
ax.set_xlabel(r'ER Ca$^{2+}$ concentration ($\mu$M)')
ax.legend(loc='best', ncol=2, bbox_to_anchor=(0.1, 0.9, 0.8, 0.2), mode='expand')

"""
S12, S21, S1a and S2a in the S2-KO model
"""
ax = fig.add_subplot(spec[1,0])
ax.plot(c_e, S1aKO, 'tab:olive', ls =  'solid', lw=2, label = r'$S_{1}^*$')
ax.plot(c_e, S11KO, 'tab:cyan', ls =  'solid', lw=2, label = r'$S_{11}^*$')
ax.plot(c_e, S2KO_CRAC, 'tab:blue', ls =  'solid', lw=2, label = r'S2KO CRAC')
ax.plot(c_e, S2KOBase, 'tab:blue', ls =  'dashed', lw=2, label = r'S2KO CRAC (base)')
#Ax limits
ax.set_ylim([-0.01, 1.4])
ax.set_yticks([0,1])
ax.set_xlim([0,1000])
ax.set_xticks([0,250, 500, 750, 1000])
ax.set_title('B', loc ='left')
ax.set_xlabel(r'ER Ca$^{2+}$ concentration ($\mu$M)')
ax.legend(loc='best', ncol=2, bbox_to_anchor=(0.1, 0.9, 0.8, 0.2), mode='expand')
"""
S12, S21, S1a and S2a in the S1-KO model
"""
ax = fig.add_subplot(spec[2,0])
ax.plot(c_e, S2aKO, 'tab:brown', ls =  'solid', lw=2, label = r'$S_{2}^*$')
ax.plot(c_e, S22KO, 'tab:pink', ls =  'solid', lw=2, label = r'$S_{22}^*$')
ax.plot(c_e, S1KO_CRAC, 'tab:red', ls =  'solid', lw=2, label = r'S1KO CRAC')
ax.plot(c_e, S1KOBase, 'tab:red', ls =  'dashed', lw=2, label = r'S1KO CRAC (base)')
#Ax limits
ax.set_ylim([-0.01, 1.4])
ax.set_yticks([0,1])
ax.set_xlim([0,1000])
ax.set_xticks([0,250, 500, 750, 1000])
ax.set_title('C', loc ='left')
ax.set_xlabel(r'ER Ca$^{2+}$ concentration ($\mu$M)')
ax.legend(loc='best', ncol=2, bbox_to_anchor=(0.1, 0.9, 0.8, 0.2), mode='expand')

filename = 'Fig10Extended_noCoop.pdf'
saveFigure(filename, fig)
# %%
