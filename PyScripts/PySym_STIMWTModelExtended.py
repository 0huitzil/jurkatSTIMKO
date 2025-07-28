#%%
"""
"""
import numpy as np
from sympy import *
from sympy.abc import x, y, z, a, b
from sympy.plotting import plot
from PyModels import parExtendedJurkatCell
par = parExtendedJurkatCell()
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
n=Symbol('n')
# Reaction rates - numerical
# k1=par['K1']
# k2=par['K2']
# k11=par['K11']
# k22=par['K22']
# k12=par['K12']
# k21=par['K21']
# n=par['n']
#Calcium
ce=Symbol('c_e')

f1=exp(ce*n)
#%%
"""
Chemical reactions at equilibrium
"""
eq1 = S1 - exp(-n*k1)*f1*S1a
eq2 = S2 - exp(-n*k2)*f1*S2a
eq3 = S1a**2 - exp(-n*k11)*S11
eq4 = S2a**2 - exp(-n*k22)*S22
eq5 = S2a*S1 - exp(-n*k21)*S21
eq6 = S1a*S2a - exp(-n*k12)*S12
#%%
"""
Total STIM
"""
eq7 = S1t - (S1 + S1a + S12 + S21 + 2*S11)
eq8 = S2t - (S2 + S2a + S12 + S21 + 2*S22)

#%%
"""
STIM dimers as functions of ce (using S1a and S2a as pivots)
The system of equations can be reduced to a smaller system of two equations on S1a and S2a, given by eq7 and eq8
"""
s1 = solveset(eq1, S1).args[0]
s2 = solveset(eq2, S2).args[0]

s11 = solveset(eq3, S11).args[0]
s22 = solveset(eq4, S22).args[0]

s21 = solveset(eq5, S21).args[0]
s12 = solveset(eq6, S12).args[0]

s21 = s21.subs(S1, s1).subs(S2, s2)
s12 = s12.subs(S1, s1).subs(S2, s2)



# eq5 = eq5.subs(S1, s1).subs(S2, s2)
# eq6 = eq6.subs(S1, s1).subs(S2, s2)



# eq3 = eq3.subs(S1, s1)
# eq4 = eq4.subs(S2, s2)






eq7 = eq7.subs(S1, s1).subs(S2, s2).subs(S12, s12).subs(S21, s21).subs(S11, s11).subs(S22, s22)
eq8 = eq8.subs(S1, s1).subs(S2, s2).subs(S12, s12).subs(S21, s21).subs(S11, s11).subs(S22, s22)
#%%
"""
Equations eq7 and eq8 can be simplifed by introducing the following phi expressions, the updated equations are named p1a and p2a, respectively
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
# Rewritting the expressions as functions of S1a and S2a
p1a = eq7.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
p2a = eq8.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
#%%
"""
p2a can be solved in terms of S2a (S2astar). We substitute this solution into p1a, which results in a quartic polynomial in S1a (S1astar). The roots of this polynomial hold the required S1a value. 
"""
# Solve eq7 for S2a
S2aStar = solveset(p1a, S2a).args[0]
# Substitute solution in eq8, leaving us with one quartic polynomial in one variable
S1aStar = p2a.subs(S2a, S2aStar).simplify()
# Only need the numerator 
S1aStar = S1aStar.args[2]
#Polynomial version?
S1aPoly = S1aStar.subs(phi4, ph4).subs(phi5, ph5).subs(phi6, ph6).subs(phi7, ph7)
S1aPoly = S1aPoly.subs(ce, 891)
S1aPoly = S1aPoly.as_poly(S1a)
np.roots(S1aPoly.all_coeffs())
#%%
"""
The expressions for all the species are printed (with the phi expressions substitued in) to be used in the figure scripts. S1a is obtained numerically from implicit function, S2a is a function of S1a, and the four reminaing species are functions of S1a and S2a
"""
S1aPrint = S1aStar
S2aPrint = S2aStar
S11Print = s11.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
S22Print = s22.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
S12Print = s12.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
S21Print = s21.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
S1Print = s1.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
S2Print = s2.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4).subs(ph5, phi5).subs(ph6, phi6).subs(ph7, phi7)
print('S1a = ' + str(S1aPrint))
print('S2a = ' + str(S2aPrint))
print('S11 = ' + str(S11Print))
print('S22 = ' + str(S22Print))
print('S12 = ' + str(S12Print))
print('S21 = ' + str(S21Print))
print('S1 = ' + str(S1Print))
print('S2 = ' + str(S2Print))
# %%
