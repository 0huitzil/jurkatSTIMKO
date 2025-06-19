#%%
"""
This file is used to create the state derivatives for the f90 file. 
However, these derivates cannot be copied 1 to 1 to the f90 file 
Due to fortran size limitations, the line breaks (&) have to be inserted manually
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
def pythonPar():
    """
    Dictionary of named parameters
    Maybe I should pass this function to a separate PySym_Models file
     
    Parameters 
    ----------
    none

    Returns
    ----------
    none
    """
    par = {
    'Vpm': 'Vpm',
    'Kpm': 'Kpm',
    's1': 's1',
    'Vsoce': 'Vsoce',
    'Ts': 'Ts',
    'Ke': 'Ke',
    'a0': 'a0',
    'n': 'n',
    'K1': 'K1',
    'K2': 'K2',
    'K12': 'K12',
    'K21': 'K21',
    'w1': 'w1',
    'w2': 'w2',
    'w12': 'w12',
    'w21': 'w21',
    'Vdeg': 'Vdeg',
    'Kdeg': 'Kdeg',
    'Vplc': 'Vplc',
    'Tp': 'Tp',
    'gammaE': 'gammaE',
    'gammaM': 'gammaM',
    'delta': 'delta',
    'Ct': 'Ct',
    'Kf': 'Kf',
    'kbeta': 'kbeta',
    'Kp': 'Kp',
    'Kc': 'Kc',
    'Kh': 'Kh',
    'Tmax': 'Tmax',
    'Kt': 'Kt',
    'VsC': 'VsC',
    'VsM': 'VsM',
    'KbarC': 'KbarC',
    'KbarM': 'KbarM',
    'KsM': 'KsM',
    'KsC': 'KsC',
    'Tg': 'Tg',
    'kdiff': 'kdiff',
    'Ki': 'Ki',
    'si': 'si',
    'Ti': 'Ti',
    'Vk': 'Vk',
    'k1Rest': 'k1Rest',
    'kn': 'kn',
    'nk': 'nk',
    'Tkmax': 'Tkmax',
    'taun': 'taun'}
    return par 
def pythonVar():
    """
    Dictionary of named variables
     
    Parameters 
    ----------
    none

    Returns
    ----------
    none
    """
    par = {
        'c': 'c',
        'ce': 'ce',
        'cm': 'cm',
        'h': 'h',
        'p': 'p',
        's': 's',
        'i': 'i',
    }
    return par
par = pythonPar()
"""
Now beings the very elegant process of initializing every symbol individually
I swear this particular implementation was useful at some point
"""
Vpm=Symbol(par["Vpm"],)
Kpm=Symbol(par["Kpm"],)
s1=Symbol(par["s1"],)
Vsoce=Symbol(par["Vsoce"],)
Ts=Symbol(par["Ts"],)
Ke=Symbol(par["Ke"],)
a0=Symbol(par["a0"],)
n=Symbol(par["n"],)
K1=Symbol(par["K1"],)
K2=Symbol(par["K2"],)
K12=Symbol(par["K12"],)
K21=Symbol(par["K21"],)
w1=Symbol(par["w1"],)
w2=Symbol(par["w2"],)
w12=Symbol(par["w12"],)
w21=Symbol(par["w21"],)
Vdeg=Symbol(par["Vdeg"],)
Kdeg=Symbol(par["Kdeg"],)
Vplc=Symbol(par["Vplc"],)
Tp=Symbol(par["Tp"],)
gammaE=Symbol(par["gammaE"],)
gammaM=Symbol(par["gammaM"],)
delta=Symbol(par["delta"],)
Ct=Symbol(par["Ct"],)
Kf=Symbol(par["Kf"],)
kbeta=Symbol(par["kbeta"],)
Kp=Symbol(par["Kp"],)
Kc=Symbol(par["Kc"],)
Kh=Symbol(par["Kh"],)
Tmax=Symbol(par["Tmax"],)
Kt=Symbol(par["Kt"],)
VsC=Symbol(par["VsC"],)
VsM=Symbol(par["VsM"],)
KbarC=Symbol(par["KbarC"],)
KbarM=Symbol(par["KbarM"],)
KsM=Symbol(par["KsM"],)
KsC=Symbol(par["KsC"],)
Tg=Symbol(par["Tg"],)
kdiff=Symbol(par["kdiff"],)
Ki=Symbol(par["Ki"],)
si=Symbol(par["si"],)
Ti=Symbol(par["Ti"],)
Vk=Symbol(par["Vk"],)
k1Rest=Symbol(par["k1Rest"],)
kn=Symbol(par["kn"],)
nk=Symbol(par["nk"],)
Tkmax=Symbol(par["Tkmax"],)
taun=Symbol(par["taun"],)
"""
Initializing variables
"""
par = pythonVar()
c = Symbol(par['c'],)
ce = Symbol(par['ce'],)
cm = Symbol(par['cm'],)
h = Symbol(par['h'],)
p = Symbol(par['p'],)
s = Symbol(par['s'],)
i = Symbol(par['i'],)
"""
Equations
"""
#======================================
#IPR flux
ma = c**4 / (Kc**4 + c**4)
B = p**2 / (Kp**2 + p**2)
A = 1 - p**2 / (Kp**2 + p**2)
ha = Kh**4 / (Kh**4 + c**4)
alpha = A*(1-ma*ha)
beta = B*ma*h
P0 = beta/(beta + kbeta*(beta + alpha))   
Jipr=Kf*P0*(ce-c)
#======================================
#h fluxes
Jh=(Kh**4 / (Kh**4 + c**4)) - h
Th=Tmax*(Kt**4 / (Kt**4 + c**4))
#======================================
#SERCA flux
JsercaC=(VsC*(c**2*(1-Tg) - KbarC*ce**2))/(c**2 + KsC**2)
JsercaM=(VsM*(cm**2*(1-Tg) - KbarM*ce**2))/(cm**2 + KsM**2)
#======================================
#IP3 flux
IPdeg=Vdeg*c**2/(c**2+Kdeg**2)
#======================================
#PMCA flux
Jpm = (Vpm*c**2)/(c**2 + Kpm**2)
#======================================
#CRAC flux
phi1=exp(n*(-K1 + ce))
phi2=exp(n*(-K2 + ce))
phi3=exp(n*(-K1 + K21 + ce))
phi4=exp(K12*n)
num1=(phi1 + 1)*(phi2 + 1)*(phi1*phi2 + phi1 + phi2 + 4*phi3 + 4*phi4 + 1)
num2=2*(phi2*phi3 + phi2*phi4 + phi3 + phi4)
num3= (phi1 + 1)/(2*(phi3 + phi4))
S2a=sqrt(num1)/num2 - num3
S1a=((phi2+1)/(phi1+1))*S2a
S21=S1a*S2a*phi3
S12=S1a*S2a*phi4
JcracWT = w1*S1a + w2*S2a + w12*S12 + w21*S21
JcracS1 = w1/(1+exp(n*(ce-K1)))
JcracS2 = w2/(1+exp(n*(ce-K2)))
#======================================
#Linear diffusion
Jdiff = kdiff*(cm-c)
#======================================
#CDI 
cdiInfty = 1/(1+exp(si*(cm-Ki)))
#======================================
#Dynamic K1
kInfty = (Vk*kn**8)/(kn**8 + ce**8) + k1Rest
TauK=Tkmax*(ce**8 / (taun**8 + ce**8))
#======================================
#%% WT cell
"""
Create the equations for the WT cell 
The state derivatives are stored in the DxF object
"""
# Equations WT cell
dc=(Jipr - JsercaC + Jdiff) + delta*(-Jpm)
dce=gammaE*(JsercaC + JsercaM - Jipr)
dcm=gammaM*(delta*s - JsercaM - Jdiff)
dh=Jh/Th
dp=(Vplc - IPdeg*p)/Tp
ds=(JcracWT - s)/Ts
F = Matrix([dc, dce, dcm, dh, dp, ds])
x = Matrix([c, ce, cm, h, p, s])
DxF = F.jacobian(x)

"""
This loop prints every state derivative using the AUTO DFDU format, to be copied in the f90 file
"""
for row in range(6):
    for col in range(6):
        der = DxF[row,col]
        print('DFDU(' + str(row+1) + ',' +  str(col+1) + ')= ' + str(der))

"""
Parameter derivatives for WT Cell, using the AUTO DFDP format
The parameter derivatives are only calculated for the Vplc parameter (par 21)
"""
x = Matrix([Vplc])
DxF = F.jacobian(x)
for row in range(6):
    der = DxF[row,0]
    print('DFDP(' + str(row+1) + ',' +  str(21) + ')= ' + str(der))
#%% S1 cell
# Equations S1 cell
"""
Create the equations for the S1 (S2KO) cell 
The state derivatives are stored in the DxF object
"""
dc=(Jipr - JsercaC + Jdiff) + delta*(-Jpm)
dce=gammaE*(JsercaC + JsercaM - Jipr)
dcm=gammaM*(delta*s - JsercaM - Jdiff)
dh=Jh/Th
dp=(Vplc - IPdeg*p)/Tp
ds=(JcracS1 - s)/Ts
di=(cdiInfty - i)/Ti
F = Matrix([dc, dce, dcm, dh, dp, ds])
x = Matrix([c, ce, cm, h, p, s])
DxF = F.jacobian(x)
"""
State derivatives for S1 (S2KO) Cell
"""
for row in range(6):
    for col in range(6):
        der = DxF[row,col]
        print('DFDU(' + str(row+1) + ',' +  str(col+1) + ')= ' + str(der))
"""
Parameter derivatives for S1 (S2KO) Cell, using the AUTO DFDP format
The parameter derivatives are only calculated for the Vplc parameter (par 21)
"""
x = Matrix([Vplc])
DxF = F.jacobian(x)
for row in range(6):
    der = DxF[row,0]
    print('DFDP(' + str(row+1) + ',' +  str(21) + ')= ' + str(der))
#%%S1 cell
# Equations S2 cell
dc=(Jipr - JsercaC + Jdiff) + delta*(-Jpm)
dce=gammaE*(JsercaC + JsercaM - Jipr)
dcm=gammaM*(delta*s - JsercaM - Jdiff)
dh=Jh/Th
dp=(Vplc - IPdeg*p)/Tp
ds=(JcracS2 - s)/Ts
di=(cdiInfty - i)/Ti
F = Matrix([dc, dce, dcm, dh, dp, ds])
x = Matrix([c, ce, cm, h, p, s])
DxF = F.jacobian(x)
"""
State derivatives for S2 (S1KO) Cell
"""
for row in range(6):
    for col in range(6):
        der = DxF[row,col]
        print('DFDU(' + str(row+1) + ',' +  str(col+1) + ')= ' + str(der))
"""
Parameter derivatives for S2 (S1KO) Cell, using the AUTO DFDP format
The parameter derivatives are only calculated for the Vplc parameter (par 21)
"""
x = Matrix([Vplc])
DxF = F.jacobian(x)
for row in range(6):
    der = DxF[row,0]
    print('DFDP(' + str(row+1) + ',' +  str(21) + ')= ' + str(der))
