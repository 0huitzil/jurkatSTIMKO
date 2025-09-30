#%%
"""
"""
from sympy import *
from sympy.abc import x, y, z, a, b
from sympy.plotting import plot
from PyModels import parExtendedJurkatCell
init_printing() 

#%% 
"""
Parameters 
"""
def latexPar():
    par = {
    # 'Vpm': [r'V_{pm}', '\muMsAp' ],
    # 'Kpm': [r'K_{pm}', '\muMAp' ],
    # 'Ts': [r'\tau_s', '(s)'  ],
    # 'a0': [r'\alpha_0', '\muMsAp'  ],
    'n': [r'n', ' ' ],
    'K1': [r'K_1', '\muMAp'  ],
    'K11': [r'K_{11}', '\muMAp' ],
    'K2': [r'K_2', '\muMAp' ],
    'K22': [r'K_{22}', '\muMAp' ],
    'K12': [r'K_{12}', '\muMAp' ],
    'K21': [r'K_{21}', '\muMAp' ],
    'w1': [r'w_1', '\muMsAp' ],
    'w11': [r'w_{11}', '\muMsAp' ],
    'w2': [r'w_2', '\muMsAp' ],
    'w22': [r'w_{22}', '\muMsAp' ],
    'w12': [r'w_{12}', '\muMsAp' ],
    'w21': [r'w_{21}', '\muMsAp' ],
    # 'Vdeg': [r'L_\degr', '\muMsAp' ],
    # 'Kdeg': [r'K_\degr', '\muMAp' ],
    # # 'Vplc': [r'V_\plc', '\muMsAp' ],
    # 'Tp': [r'\tau_p', '(s)' ],
    # 'gammaE': [r'\gamma_E', ' ' ],
    # 'gammaM': [r'\gamma_M', ' ' ],
    # 'delta': [r'\delta', ' ' ],
    # 'Kf': [r'K_f', '\muMsAp' ],
    # 'kbeta': [r'k_\beta', ' ' ],
    # 'Kp': [r'K_p', '\muMAp' ],
    # 'Kc': [r'K_c', '\muMAp' ],
    # 'Kh': [r'K_h', '\muMAp' ],
    # 'Tmax': [r'\tau_{max}', '(s)' ],
    # 'Kt': [r'K_\tau', '\muMAp' ],
    # 'VsC': [r'V_{sC}', '\muMsAp' ],
    # 'VsM': [r'V_{sM}', '\muMsAp' ],
    # 'KbarC': [r'\bar{K}_C', ' ' ],
    # 'KbarM': [r'\bar{K}_M', ' ' ],
    # 'KsM': [r'K_{sM}', '\muMAp' ],
    # 'KsC': [r'K_{sC}', '\muMAp' ],
    # 'Tg': [r'T_g', ' ' ],
    # 'kdiff': [r'k_{diff}', '\muMsTwoAp' ],
    # 'Ki': [r'K_i', '\muMAp' ],
    # 'si': [r's_i', ' ' ],
    # 'Ti': [r'T_i', '(s)' ],
    # 'Vk': [r'V_k', ' ' ],
    # 'k1Rest': [r'k1_{Rest}', '\muMAp' ],
    # 'kn': [r'k_n',  '\muMAp'],
    # 'nk': [r'n_k', ' ' ],
    # 'Tkmax': [r'\tau_{kmax}', '(s)' ],
    # 'taun': [r'\tau_n', '(s)']
    }
    return par 
#%%
pars1 = []
pars2 = []
pars3 = []
vals1 = []
vals2 = []
vals3 = []
unit1 = []
unit2 = []
unit3 = []
parTex = latexPar()
parVals = parExtendedJurkatCell()
headers = ['Par', 'Value', '(Units)']*3
for k, p in list(enumerate(parTex)):
    if k < 5:
        pars1.append('$' + parTex[p][0] + '$')
        vals1.append(str(parVals[p]))
        unit1.append(' $' + parTex[p][1] + '$')
    elif k < 10:
        pars2.append('$' + parTex[p][0]+ '$')
        vals2.append(str(parVals[p]))
        unit2.append(' $' + parTex[p][1] + '$')
    else:
        pars3.append('$' + parTex[p][0] + '$')
        vals3.append(str(parVals[p]))
        unit3.append(' $' + parTex[p][1] + '$')
# %%

def theader(hcols):
  "Given a list of column headers, return the header lines of a LaTeX table."

  # Figure out the maximum number of lines in the header
  hcols = [ x.splitlines() for x in hcols ]
  maxlines = max(len(x) for x in hcols)
  hcols = [ [' ']*(maxlines - len(x)) + x for x in hcols ]

  # Assemble the rows of the header
  rows = []
  for i in range(maxlines):
    row = [ str(c[i]) for c in hcols ]
    rows.append(' & '.join(row) + r' \\')

  return '\n'.join(rows) + '\n\\midrule'

def tbody(*cols, group=0):
  "Given the columns, return the body of a LaTeX table."

  # Add blanks at the ends of short columns
  lens = [ len(c) for c in cols ]
  nrows = max(lens)
  cols = [ c + [' ']*(nrows - len(c)) for c in cols ]

  # Assemble the rows of the table
  rows = []
  for i in range(nrows):
    if (group > 0) and (i > 0) and (i % group == 0):
      rows.append(r'\addlinespace')
    row = [ str(c[i]) for c in cols ]
    rows.append(' & '.join(row) + r' \\')

  return '\n'.join(rows)
#%%
print(theader(headers))
print(tbody(pars1, vals1, unit1, pars2, vals2, unit2, pars3, vals3, unit3))
# %%
