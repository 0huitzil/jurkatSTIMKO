#%%
"""
This file is used to create the s_infty expressions for the WT model
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
n1=Symbol('n_1')
n2=Symbol('n_2')
f1=exp(ce*n1)
f2=exp(ce*n2)
#%%
"""
Chemical reactions at equilibrium (Regular coefficients)
"""
eq1 = S1 - k1*f1*S1a
eq2 = S2 - k2*f2*S2a
eq3 = S2a*S1 - k21*S21
eq4 = S1a*S2a - k12*S12
#%%
"""
Chemical reactions at equilibrium (Exponential coefficients)
"""
eq1 = S1 - exp(-n1*k1)*f1*S1a
eq2 = S2 - exp(-n2*k2)*f2*S2a
eq3 = S2a*S1 - exp(-n1*k21)*S21
eq4 = S1a*S2a - exp(-n1*k12)*S12
#%%
"""
Total STIM
"""
eq5 = S1t - (S1 + S1a + S12 + S21 )
eq6 = S2t - (S2 + S2a + S12 + S21 )

#%%
"""
STIM dimers as functions of ce (using S1a and S2a as pivots)
"""
s1 = solveset(eq1, S1).args[0].simplify()
s2 = solveset(eq2, S2).args[0].simplify()


eq3 = eq3.subs(S1, s1).subs(S2, s2)
eq4 = eq4.subs(S1, s1).subs(S2, s2)

s21 = solveset(eq3, S21).args[0].simplify()
s12 = solveset(eq4, S12).args[0].simplify()

eq5 = eq5.subs(S1, s1).subs(S2, s2).subs(S12, s12).subs(S21, s21)
eq6 = eq6.subs(S1, s1).subs(S2, s2).subs(S12, s12).subs(S21, s21)
#%%
#%%
"""
Polynomials
"""
phi1 = Symbol('\phi_1')
phi2 = Symbol('\phi_2')
phi3 = Symbol('\phi_3')
phi4 = Symbol('\phi_4')
# Phi substitution
ph1 = exp(n1*(ce - k1))
ph2 = exp(n2*(ce - k2))
ph3 = exp(n1*(ce + k21 - k1))
ph4 = exp(n1*(k12))
# Rewritting the expressions as polynomials of S1a and S2a
p1a = eq5.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4)
p2a = eq6.subs(ph1, phi1).subs(ph2, phi2).subs(ph3, phi3).subs(ph4, phi4)
p1 = p1a.as_poly(S1a, S2a)
p2 = p2a.as_poly(S1a, S2a)

s= solve_poly_system([p1, p2])
s = s[1]
#%%
"""
Explicit forms
"""
S1Astar = s[0]
S2Astar = s[1]
S21star = s21.subs(ph3, phi3).subs(ph4, phi4)
S12star = s12.subs(ph3, phi3).subs(ph4, phi4)
S1star = s1.subs(ph1, phi1).subs(ph2, phi2)
S2star = s2.subs(ph1, phi1).subs(ph2, phi2)

#%% Python printing
# 
print('phi1=')
print(ph1)
print('phi2=')
print(ph2)
print('phi3=')
print(ph3)
print('phi4=')
print(ph4)
print('S1a=')
print(S1Astar)
print('S2a=')
print(S2Astar)
print('S1=')
print(S1star)
print('S2=')
print(S2star)
print('S21=')
print(S21star)
print('S12=')
print(S12star)
# %% Desmos printing 
"""
Unused, but left here just in case. This code changes the variable names to 'x' and 'y' in order to plot them in Desmos. 
"""
# Desmos_S1AStar = S1Astar.subs(ce, 'x').subs(phi1, 'phi_1(x)').subs(phi2, 'phi_2(x)').subs(phi3, 'phi_3(x)').subs(phi4, 'phi_4(x)')
# Desmos_S2AStar = S2Astar.subs(ce, 'x').subs(phi1, 'phi_1(x)').subs(phi2, 'phi_2(x)').subs(phi3, 'phi_3(x)').subs(phi4, 'phi_4(x)')
# Desmos_S12AStar = S12star.subs(ce, 'x').subs(phi1, 'phi_1(x)').subs(phi2, 'phi_2(x)').subs(phi3, 'phi_3(x)').subs(phi4, 'phi_4(x)').subs(S1a, 'S_1a(x)').subs(S2a, 'S_2a(x)')
# Desmos_S21AStar = S21star.subs(ce, 'x').subs(phi1, 'phi_1(x)').subs(phi2, 'phi_2(x)').subs(phi3, 'phi_3(x)').subs(phi4, 'phi_4(x)').subs(S1a, 'S_1a(x)').subs(S2a, 'S_2a(x)')
# print_latex(Desmos_S1AStar)
# print_latex(Desmos_S2AStar)
# print_latex(Desmos_S12AStar)
# print_latex(Desmos_S21AStar)

# print_latex(s1.subs(ce, 'x').subs(S1a, 'S_1a(x)').subs(S2a, 'S_2a(x)'))
# print_latex(s2.subs(ce, 'x').subs(S1a, 'S_1a(x)').subs(S2a, 'S_2a(x)'))


# print_latex(s21.subs(ce, 'x').subs(S1a, 'S_1a(x)').subs(S2a, 'S_2a(x)'))
# print_latex(s12.subs(ce, 'x').subs(S1a, 'S_1a(x)').subs(S2a, 'S_2a(x)'))

# print_latex(ph1.subs(ce, 'x').subs(S1a, 'S_1a(x)').subs(S2a, 'S_2a(x)'))
# print_latex(ph2.subs(ce, 'x').subs(S1a, 'S_1a(x)').subs(S2a, 'S_2a(x)'))


# print_latex(ph3.subs(ce, 'x').subs(S1a, 'S_1a(x)').subs(S2a, 'S_2a(x)'))
# print_latex(ph4.subs(ce, 'x').subs(S1a, 'S_1a(x)').subs(S2a, 'S_2a(x)'))



