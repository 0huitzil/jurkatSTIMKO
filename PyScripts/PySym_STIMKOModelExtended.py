#%%
"""
"""
from sympy import *
from sympy.abc import x, y, z, a, b
from sympy.plotting import plot

init_printing() 

#%% 
"""
Parameters 
"""
S1=Symbol('S1')
S1a=Symbol('S1a')
S2=Symbol('S2')
S2a=Symbol('S2a')
S11=Symbol('S11')
S12=Symbol('S12')
S21=Symbol('S21')
S22=Symbol('S22')
#Total STIM
S1t=1
S2t=1
# Reaction rates 
k1=Symbol('K_1')
k2=Symbol('K_2')
k11=Symbol('K_11')
k22=Symbol('K_22')
k12=Symbol('K_12')
k21=Symbol('K_21')
#Calcium
ce=Symbol('c_e')
n=Symbol('n')
f1=exp(ce*n)
#%%
"""
STIM-2 KO experiments 
"""
eq1 = S1 - exp(-n*k1)*f1*S1a
eq2 = S1a**2 - exp(-n*k11)*S11
#%%
"""
Total STIM
"""
eq3 = S1t - (S1 + S1a + 2*S11)
#%%
"""
Solving eq3
"""
s1 = solveset(eq1, S1).args[0]
s11 = solveset(eq2, S11).args[0]
eq3 = eq3.subs(S1, s1).subs(S11, s11)
#%%
"""
EQ3 can be simplifed by introducing the following phi expressions
"""
phi1 = Symbol('phi_1')
phi2 = Symbol('phi_2')
phi3 = Symbol('phi_3')
phi4 = Symbol('phi_4')
phi5 = Symbol('phi_5')
phi6 = Symbol('phi_6')
phi7 = Symbol('phi_7')
# Phi substitution
ph1 = exp(n*(ce - k1)).expand()
ph2 = exp(n*(ce - k2)).expand()
ph3 = exp(n*(ce + k21 - k1)).expand() #might not be necessary?
ph4 = exp(n*(k12)).expand()
ph5 = exp(n*(k22)).expand()
ph6 = exp(n*(k11)).expand()
ph7 = exp(n*(k21)).expand()

eq3 = eq3.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
#%%
"""
EQ3 can be solved symbolically, which gives up the species expressions
"""
S1aStar = solveset(eq3, S1a).args[1]
S11Star = s11.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
#%%
"""
The expressions for all the species are printed (with the phi expressions substitued in) to be used in the figure scripts.
"""
print("\n\nSTIM2-KO expressions \n\n")
print('S1aKO = ' + str(S1aStar))
print('S11KO = ' + str(S11Star))

#%%
"""
STIM-1 KO experiments 
"""
eq1 = S2 - exp(-n*k2)*f1*S2a
eq2 = S2a**2 - exp(-n*k22)*S22
#%%
"""
Total STIM
"""
eq3 = S2t - (S2 + S2a + 2*S22)
#%%
"""
STIM dimers as functions of ce (using S2a and S2a as pivots)
"""
s2 = solveset(eq1, S2).args[0]
s22 = solveset(eq2, S22).args[0]
eq3 = eq3.subs(S2, s2).subs(S22, s22)
#%%
"""
EQ3 can be simplifed by introducing the following phi expressions
"""
phi1 = Symbol('phi_1')
phi2 = Symbol('phi_2')
phi3 = Symbol('phi_3')
phi4 = Symbol('phi_4')
phi5 = Symbol('phi_5')
phi6 = Symbol('phi_6')
phi7 = Symbol('phi_7')
# Phi substitution
ph1 = exp(n*(ce - k1)).expand()
ph2 = exp(n*(ce - k2)).expand()
ph3 = exp(n*(ce + k21 - k1)).expand() #might not be necessary?
ph4 = exp(n*(k12)).expand()
ph5 = exp(n*(k22)).expand()
ph6 = exp(n*(k11)).expand()
ph7 = exp(n*(k21)).expand()

eq3 = eq3.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
#%%
"""
EQ3 can be solved symbolically, which gives up the species expressions
"""
S2aStar = solveset(eq3, S2a).args[1]
S22Star = s22.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
#%%
"""
The expressions for all the species are printed (with the phi expressions substitued in) to be used in the figure scripts.
"""
print("\n\nSTIM1-KO expressions \n\n")
print('S2aKO = ' + str(S2aStar))
print('S22KO = ' + str(S22Star))

# %%
