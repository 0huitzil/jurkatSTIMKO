#%%
"""
"""
from sympy import *
from sympy.abc import x, y, z, a, b
from sympy.plotting import plot
from PyModels import parJurkatCell
init_printing() 

#%% 
"""
Parameters 
"""
def latexPar():
    par = {
    'Vpm': r'V_{pm}',
    'Kpm': r'K_{pm}',
    'Ts': r'\tau_s',
    'a0': r'\alpha_0',
    'n': r'n',
    'K1': r'K_1',
    'K2': r'K_2',
    'K12': r'K_{12}',
    'K21': r'K_{21}',
    'w1': r'w_1',
    'w2': r'w_2',
    'w12': r'w_{12}',
    'w21': r'w_{21}',
    'Vdeg': r'V_\degr',
    'Kdeg': r'K_\degr',
    'Vplc': r'V_\plc',
    'Tp': r'\tau_p',
    'gammaE': r'\gamma_E',
    'gammaM': r'\gamma_M',
    'delta': r'\delta',
    'Kf': r'K_f',
    'kbeta': r'k_\beta',
    'Kp': r'K_p',
    'Kc': r'K_c',
    'Kh': r'K_h',
    'Tmax': r'\tau_{max}',
    'Kt': r'K_t',
    'VsC': r'V_{sC}',
    'VsM': r'V_{sM}',
    'KbarC': r'\bar{K}_C',
    'KbarM': r'\bar{K}_M',
    'KsM': r'K_{sM}',
    'KsC': r'K_{sC}',
    'Tg': r'T_g',
    'kdiff': r'k_{diff}',
    'Ki': r'K_i',
    'si': r's_i',
    'Ti': r'T_i',
    'Vk': r'V_k',
    'k1Rest': r'k1_{Rest}',
    'kn': r'k_n',
    'nk': r'n_k',
    'Tkmax': r'\tau_{kmax}',
    'taun': r'\tau_n'}
    return par 
#%%
pars1 = []
pars2 = []
pars3 = []
vals1 = []
vals2 = []
vals3 = []
parTex = latexPar()
parVals = parJurkatCell()
headers = ['Par', 'Value']*3
for k, p in list(enumerate(parTex)):
    if k < 15:
        pars1.append('$' + parTex[p] + '$')
        vals1.append(parVals[p])
    elif k < 30:
        pars2.append('$' + parTex[p]+ '$')
        vals2.append(parVals[p])
    else:
        pars3.append('$' + parTex[p]+ '$')
        vals3.append(parVals[p])
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
print(tbody(pars1, vals1, pars2, vals2, pars3, vals3))