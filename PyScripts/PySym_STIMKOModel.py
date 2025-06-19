#%%
"""
This file is used to create the s_infty expressions for the KO models
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
S1a=Symbol('S1^*')
S2=Symbol('S2')
S2a=Symbol('S2^*')
S12=Symbol('S12^*')
S21=Symbol('S21^*')
#Total STIM
S1t=1
S2t=1
# Reaction rates 
k1=Symbol('K_1')
k2=Symbol('K_2')
k12=Symbol('K_12')
k21=Symbol('K_21')
#Calcium
ce=Symbol('c_e')
n=Symbol('n')
f=exp(ce*n)

#%% 
"""
STIM-2 KO experiments 
"""

eq1 = S1 - exp(-k1*n)*f*S1a
eq2 = S1t -(S1 + S1a)

s1 = solveset(eq1, S1).args[0]
eq2 = eq2.subs(S1, s1)
s1a = solveset(eq2, S1a).args[0]

print_latex(s1a.subs(ce, 'x')) # This last line is only used to export the expression to desmos.

#%% 
"""
STIM-1 KO experiments 
"""
eq1 = S2 - exp(-k2*n)*f*S2a
eq2 = S2t -(S2 + S2a)

s2 = solveset(eq1, S2).args[0]
eq2 = eq2.subs(S2, s2)
s2a = solveset(eq2, S2a).args[0]

print_latex(s2a.subs(ce, 'x')) # This last line is only used to export the expression to desmos.