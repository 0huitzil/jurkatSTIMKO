#%% 
# region Libraries
import numpy as np
from scipy.integrate import solve_ivp
from myOptions import *
# endregion
#%% S2KO model versions
"""
S2KO model versions
"""
def jurkatS2KOCell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Model without CDI
    Parameters
    ----------
    t: double
            time 
    y: vector 
            state variables
    par: dict 
            dictionary with parameters

    Returns 
    ----------
    field: vector 
            field equations
    """

    """
    Declare variables
    """
    c, ce, cm, h, p, s, i = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    JsercaC= getSERCAC(y, par)
    JsercaM= getSERCAM(y, par)
    Jcrac = getS2KOCRAC(y, par)
    Jpm = getPMCA(y, par)
    Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    Jcdi = getCDI(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC + Jdiff) - delta*Jpm
    dce=gammaE*(JsercaC + JsercaM - Jipr)
    dcm=gammaM*(delta*s*i - JsercaM - Jdiff)
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di=0
    field = [dc, dce, dcm, dh, dp, ds, di]
    """
    Return
    """
    return field 

def jurkatS2KOCDICell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Model with CDI enabled
    Parameters
    ----------
    t: double
            time 
    y: vector 
            state variables
    par: dict 
            dictionary with parameters

    Returns 
    ----------
    field: vector 
            field equations
    """
    
    """
    Declare variables
    """
    c, ce, cm, h, p, s, i = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    JsercaC= getSERCAC(y, par)
    JsercaM= getSERCAM(y, par)
    Jcrac = getS2KOCRAC(y, par)
    Jpm = getPMCA(y, par)
    Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    Jcdi = getCDI(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC + Jdiff) - delta*Jpm
    dce=gammaE*(JsercaC + JsercaM - Jipr)
    dcm=gammaM*(delta*s*i - JsercaM - Jdiff)
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di=(Jcdi - i)/Ti
    field = [dc, dce, dcm, dh, dp, ds, di]
    """
    Return
    """
    return field 

def jurkatBaseS2KOCell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Base model without microdomain or CDI
    Easiest way to run this version of the model without changing 
    the flux functions is just to set the cm, i equations to zero
    Parameters
    ----------
    t: double
            time 
    y: vector 
            state variables
    par: dict 
            dictionary with parameters

    Returns 
    ----------
    field: vector 
            field equations
    """

    """
    Declare variables
    """
    c, ce, cm, h, p, s, i = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    JsercaC= getSERCAC(y, par)
    # JsercaM= getSERCAM(y, par)
    Jcrac = getS2KOCRAC(y, par)
    Jpm = getPMCA(y, par)
    # Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    # Jcdi = getCDI(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC) + delta*(s - Jpm)
    dce=gammaE*(JsercaC  - Jipr)
    dcm= 0
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di= 0
    field = [dc, dce, dcm, dh, dp, ds, di]
    """
    Return
    """
    return field 

def jurkatS2KODynamicCell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Model with dynamic K1 enabled
    Parameters
    ----------
    t: double
            time 
    y: vector 
            state variables
    par: dict 
            dictionary with parameters

    Returns 
    ----------
    field: vector 
            field equations
    """

    """
    Declare variables
    """
    c, ce, cm, h, p, s, i, K1 = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    Tauk = getTauK(y, par)
    JsercaC= getSERCAC(y, par)
    JsercaM= getSERCAM(y, par)
    Jcrac = getS2KOCRAC(y, par)
    Jpm = getPMCA(y, par)
    Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    Jcdi = getCDI(y, par)
    kCurve = getKcurve(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC + Jdiff) - delta*Jpm
    dce=gammaE*(JsercaC + JsercaM - Jipr)
    dcm=gammaM*(delta*s*i - JsercaM - Jdiff)
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di=(Jcdi - i)/Ti
    dK1=(kCurve-K1)/Tauk 
    field = [dc, dce, dcm, dh, dp, ds, di, dK1]
    """
    Return
    """
    return field 
#%% S1KO model versions
"""
S1KO model versions
There is not S1KODynamicCell
"""
def jurkatS1KOCell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Model without CDI
    Parameters
    ----------
    t: double, time 
    y: vector with the state variables

    Returns 
    ----------
    field: vector with the field equations
    """
    """
    Declare variables
    """
    c, ce, cm, h, p, s, i = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    JsercaC= getSERCAC(y, par)
    JsercaM= getSERCAM(y, par)
    Jcrac = getS1KOCRAC(y, par)
    Jpm = getPMCA(y, par)
    Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    Jcdi = getCDI(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC + Jdiff) - delta*Jpm
    dce=gammaE*(JsercaC + JsercaM - Jipr)
    dcm=gammaM*(delta*s*i - JsercaM - Jdiff)
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di=0
    field = [dc, dce, dcm, dh, dp, ds, di]
    """
    Return
    """
    return field 

def jurkatS1KOCDICell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Model with CDI
    Parameters
    ----------
    t: double, time 
    y: vector with the state variables

    Returns 
    ----------
    field: vector with the field equations
    """
    """
    Declare variables
    """
    c, ce, cm, h, p, s, i = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    JsercaC= getSERCAC(y, par)
    JsercaM= getSERCAM(y, par)
    Jcrac = getS1KOCRAC(y, par)
    Jpm = getPMCA(y, par)
    Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    Jcdi = getCDI(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC + Jdiff) - delta*Jpm
    dce=gammaE*(JsercaC + JsercaM - Jipr)
    dcm=gammaM*(delta*s*i - JsercaM - Jdiff)
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di=(Jcdi - i)/Ti
    field = [dc, dce, dcm, dh, dp, ds, di]
    """
    Return
    """
    return field 

def jurkatBaseS1KOCell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Base model without microdomain or CDI
    Easiest way to run this version of the model without changing 
    the flux functions is just to set the cm, i equations to zero
    Parameters
    ----------
    t: double
            time 
    y: vector 
            state variables
    par: dict 
            dictionary with parameters

    Returns 
    ----------
    field: vector 
            field equations
    """
    """
    Declare variables
    """
    c, ce, cm, h, p, s, i = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    JsercaC= getSERCAC(y, par)
    # JsercaM= getSERCAM(y, par)
    Jcrac = getS1KOCRAC(y, par)
    Jpm = getPMCA(y, par)
    # Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    # Jcdi = getCDI(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC) + delta*(s - Jpm)
    dce=gammaE*(JsercaC  - Jipr)
    dcm= 0
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di= 0
    field = [dc, dce, dcm, dh, dp, ds, di]
    """
    Return
    """
    return field 
#%% WT model versions
"""
WT model versions
"""

def jurkatWTCell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Model without CDI
    Parameters
    ----------
    t: double
            time 
    y: vector 
            state variables
    par: dict 
            dictionary with parameters

    Returns 
    ----------
    field: vector 
            field equations
    """
    
    """
    Declare variables
    """
    c, ce, cm, h, p, s, i = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    JsercaC= getSERCAC(y, par)
    JsercaM= getSERCAM(y, par)
    Jcrac = getWTCRAC(y, par)
    Jpm = getPMCA(y, par)
    Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    Jcdi = getCDI(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC + Jdiff) - delta*Jpm
    dce=gammaE*(JsercaC + JsercaM - Jipr)
    dcm=gammaM*(delta*s - JsercaM - Jdiff)
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di=0
    field = [dc, dce, dcm, dh, dp, ds, di]
    """
    Return
    """
    return field 

def jurkatWTCDICell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Model with CDI
    Parameters
    ----------
    t: double
            time 
    y: vector 
            state variables
    par: dict 
            dictionary with parameters

    Returns 
    ----------
    field: vector 
            field equations
    """

    """
    Declare variables
    """
    c, ce, cm, h, p, s, i = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    JsercaC= getSERCAC(y, par)
    JsercaM= getSERCAM(y, par)
    Jcrac = getWTCRAC(y, par)
    Jpm = getPMCA(y, par)
    Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    Jcdi = getCDI(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC + Jdiff) - delta*Jpm
    dce=gammaE*(JsercaC + JsercaM - Jipr)
    dcm=gammaM*(delta*s*i - JsercaM - Jdiff)
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di=(Jcdi - i)/Ti
    field = [dc, dce, dcm, dh, dp, ds, di]
    """
    Return
    """
    return field 

def jurkatBaseWTCell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Base model without microdomain or CDI
    Easiest way to run this version of the model without changing 
    the flux functions is just to set the cm, i equations to zero
    Parameters
    ----------
    t: double
            time 
    y: vector 
            state variables
    par: dict 
            dictionary with parameters

    Returns 
    ----------
    field: vector 
            field equations
    """
    
    """
    Declare variables
    """
    c, ce, cm, h, p, s, i = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    JsercaC= getSERCAC(y, par)
    # JsercaM= getSERCAM(y, par)
    Jcrac = getWTCRAC(y, par)
    Jpm = getPMCA(y, par)
    Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    # Jcdi = getCDI(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC) + delta*(s -Jpm)
    dce=gammaE*(JsercaC - Jipr)
    dcm=0
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di=0
    field = [dc, dce, dcm, dh, dp, ds, di]
    """
    Return
    """
    return field 

def jurkatWTDynamicCell(t, y, par):
    """
    Create the model equations using Scipy for integration (solve_ivp)
    Model with dynamic K1 enabled
    Parameters
    ----------
    t: double
            time 
    y: vector 
            state variables
    par: dict 
            dictionary with parameters

    Returns 
    ----------
    field: vector 
            field equations
    """

    """
    Declare variables
    """
    c, ce, cm, h, p, s, i, K1 = y
    """
    Declare parameters 
    """
    #Volume
    gammaE=par['gammaE']
    gammaM=par['gammaM']
    delta = par['delta']
    Ts = par['Ts']
    Tp = par['Tp']
    Ti = par['Ti']
    """ 
    Declare fluxes 
    """
    Jipr = getIPR(y, par)
    Jh = getJh(y, par)
    Th = getTh(y, par)
    Tauk = getTauK(y, par)
    JsercaC= getSERCAC(y, par)
    JsercaM= getSERCAM(y, par)
    Jcrac = getWTCRAC(y, par)
    Jpm = getPMCA(y, par)
    Jdiff= getDiff(y, par)
    Jdegr = getDegr(y, par)
    Jcdi = getCDI(y, par)
    kCurve = getKcurve(y, par)
    """
    Declare field 
    """
    dc=(Jipr - JsercaC + Jdiff) - delta*Jpm
    dce=gammaE*(JsercaC + JsercaM - Jipr)
    dcm=gammaM*(delta*s*i - JsercaM - Jdiff)
    dh=Jh/Th
    dp=(Jdegr)/Tp
    ds=(Jcrac - s)/Ts
    di=(Jcdi - i)/Ti
    dK1=(kCurve-K1)/Tauk 
    field = [dc, dce, dcm, dh, dp, ds, di, dK1]
    """
    Return
    """
    return field 

#%% Parameters
def parJurkatCell():
    """
    This is the main parameter set used in the jurkatCell model

    Parameters
    ----------

    Returns 
    ----------
    par: dict, set of named parameters
    """
    par = {
        # PMCA 
        'Vpm': 1.0, 
        'Kpm': 0.2, 
        # SOCE 
        's1': 0.15, 
        'Vsoce': 3, 
        'Ts': 4, 
        'Ke': 800, 
        'a0': 0,
        # Advanced SOCE 
        'n': 0.15,
        'K1': 200,
        "K2": 350,
        'K12': 400,
        'K21': 470,
        'w1': 1,
        'w2': 0.7,
        'w12': 2,
        'w21': 2,
        # IP3 
        'Vdeg': 6, 
        'Kdeg': 0.5,
        'Vplc': 0, 
        'Tp': 2.9, 
        # 'p': 0, 
        # Volume
        'gammaE': 5.5,
        'gammaM': 50,
        'delta': 3,
        # Total Ca, conserved for closed cell problems 
        'Ct': 147,
        # IPR 
        'Kf': 1.6,
        'kbeta': 0.4,
        'Kp': 10,
        'Kc': 0.16,
        # H
        'Kh': 0.168,
        'Tmax': 10,
        'Kt': 0.095,
        # SERCA 
        'VsC': 1,
        'VsM': 0.2,
        'KbarC': 2.2e-8,
        'KbarM': 2.2e-8,
        'KsM': 0.2,
        'KsC': 0.2,
        'Tg': 0,
        # Diffusion
        'kdiff': 0.004,
        # CDI
        'Ki': 800, 
        'si': 0.2,
        'Ti': 1, 
        # Dynamic K1
        'Vk': 670,
        'k1Rest': 200,
        'kn': 100,
        'nk': 8,
        'Tkmax': 800,
        'taun': 60,
    }
    return par

def parBaseJurkatCell():
    """
    This is the parameter set used in the base version of the model
    e.g. the one from the previous paper
    
    Parameters
    ----------
    
    Returns 
    ----------
    par: dict, set of named parameters
    """
    par = {
        # PMCA 
        'Vpm': 3.0, 
        'Kpm': 0.2, 
        # SOCE 
        's1': 0.2, 
        'Vsoce': 3, 
        'Ts': 15, 
        'Ke': 800, 
        'a0': 0,
        # Advanced SOCE 
        'n': 0.3,
        'K1': 200,
        "K2": 350,
        'K12': 400,
        'K21': 470,
        'w1': 1,
        'w2': 0.7,
        'w12': 3,
        'w21': 3,
        # IP3 
        'Vdeg': 6, 
        'Kdeg': 0.5,
        'Vplc': 0, 
        'Tp': 2, 
        # 'p': 0, 
        # Volume
        'gammaE': 5.5,
        'gammaM': 50,
        'delta': 2,
        # Total Ca, conserved for closed cell problems 
        'Ct': 147,
        # IPR 
        'Kf': 1.6,
        'kbeta': 0.4,
        'Kp': 10,
        'Kc': 0.16,
        # H
        'Kh': 0.168,
        'Tmax': 7.5,
        'Kt': 0.095,
        # SERCA 
        'VsC': 2.0,
        'VsM': 0.2,
        'KbarC': 1e-8,
        'KbarM': 1e-8,
        'KsM': 0.19,
        'KsC': 0.19,
        'Tg': 0,
        # Diffusion
        'kdiff': 0.004,
        # CDI
        'Ki': 800, 
        'si': 0.2,
        'Ti': 1, 
        # Dynamic K1
        'Vk': 670,
        'k1Rest': 200,
        'kn': 100,
        'nk': 8,
        'Tkmax': 800,
        'taun': 60,
    }
    return par


#%% Initial conditions
"""
All the sets of initial conditions
"""
def icsJurkatCell():
    """
    These are the coordinates of the resting equilibrium point 
    at rest (Vplc = 0)
    
    Parameters
    ----------

    Returns 
    ----------
    par: dict, set of named parameters
    """
    dic = {
        #Initial conditions 
        'c': 1.00495E-01,
        'ce': 6.14590E+02, 
        'cm': 5.05751E+01,
        'h': 8.86495E-01,
        'p': 0,
        's': 2.01584E-01,
        'i': 1, 
    }
    return dic


def icsBaseJurkatCell():
    """
    These are the coordinates of the resting equilibrium point 
    at rest (Vplc = 0) of the base model 

    Parameters
    ----------

    Returns 
    ----------
    par: dict, set of named parameters
    """
    dic = {
        #Initial conditions 
        'c': 8.09050E-02,
        'ce': 8.09050E+02, 
        'cm': 100,
        'h': 9.48960E-01,
        'p': 0,
        's': 4.21884E-01,
        'i': 1, 
    }
    return dic

def icsJurkatCellDynamicK1():
    """
    These are the coordinates of the resting equilibrium point 
    at rest (Vplc = 0) of the models with dynamic K1

    Parameters
    ----------
    Returns 
    ----------
    par: dict, set of named parameters
    """
    dic = {
        #Initial conditions 
        'c': 1.00495E-01,
        'ce': 6.14590E+02, 
        'cm': 5.05751E+01,
        'h': 8.86495E-01,
        'p': 0,
        's': 2.01584E-01,
        'i': 1, 
        'K1': 200
    }
    return dic


#%% Fluxes
"""
Each function corresponds to an individual flux
"""
def getPMCA(y, par):     
    """
    The PMCA flux

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jpm: double
            The PMCA flux

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    #PM 
    Vpm=par['Vpm']
    Kpm=par['Kpm']
    # PM
    Jpm = (Vpm*c**2)/(c**2 + Kpm**2)
    return Jpm
    
def getIPR(y, par): 
    """
    The IPR flux

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jpir: double
            The PMCA flux

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    # IPR 
    Kf=par['Kf']
    kbeta=par['kbeta']
    Kp=par['Kp']
    Kc=par['Kc']
    Kh=par['Kh']
    gammaE=par['gammaE']    
    ma = c**4 / (Kc**4 + c**4)
    B = p**2 / (Kp**2 + p**2)
    A = 1 - p**2 / (Kp**2 + p**2)
    ha = Kh**4 / (Kh**4 + c**4)
    alpha = A*(1-ma*ha)
    beta = B*ma*h
    P0 = beta/(beta + kbeta*(beta + alpha))   
    Jipr=Kf*P0*(ce-c)
    return Jipr

def getP0(y, par): 
    """
    The open probability of the IPR

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jpm: double
            The PMCA flux

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    # IPR 
    Kf=par['Kf']
    kbeta=par['kbeta']
    Kp=par['Kp']
    Kc=par['Kc']
    Kh=par['Kh']
    gammaE=par['gammaE']    
    ma = c**4 / (Kc**4 + c**4)
    B = p**2 / (Kp**2 + p**2)
    A = 1 - p**2 / (Kp**2 + p**2)
    ha = Kh**4 / (Kh**4 + c**4)
    alpha = A*(1-ma*ha)
    beta = B*ma*h
    P0 = beta/(beta + kbeta*(beta + alpha))   
    return P0

def getSERCAC(y, par): 
    """
    The bulk cytoplasm SERCA flux

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jserca: double
            The SERCA flux

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    # SERCA 
    Vs=par['VsC']
    Kbar=par['KbarC']
    Ks=par['KsC']
    gammaE=par['gammaE']    
    Tg=par['Tg']
    # kleakc=par['kleakc']
    Jserca=(Vs*((1-Tg)*c**2 - Kbar*ce**2))/(c**2 + Ks**2)
    # Jserca=(Vs*((1-Tg)*c**2))/(c**2 + Ks**2) - kleakc*(ce-c) #unidirectional SERCA, depreciated
    return Jserca

def getSERCAM(y, par): 
    """
    The microdomain SERCA flux

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jserca: double
            The SERCA flux

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    # SERCA 
    Vs=par['VsM']
    Kbar=par['KbarM']
    Ks=par['KsM']
    gammaE=par['gammaE']    
    Tg=par['Tg']
    # kleakm=par['kleakm']
    Jserca=(Vs*((1-Tg)*cm**2 - Kbar*ce**2))/(cm**2 + Ks**2) 
    # Jserca=(Vs*((1-Tg)*cm**2))/(cm**2 + Ks**2) - kleakm*(cm-ce)
    return Jserca

def getJh(y, par): 
    """
    The h_infty term

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jserca: double
            The SERCA flux

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    #h_infty
    Kh=par['Kh']
    Jh=(Kh**4 / (Kh**4 + c**4)) - h
    return Jh

def getTh(y, par): 
    """
    The timescale for the h variable

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jpm: double
            The PMCA flux

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    #Tauh
    Kt=par['Kt']
    Tmax=par['Tmax']
    Th=Tmax*(Kt**4 / (Kt**4 + c**4))
    return Th

def getS2KOCRAC(y, par): 
    """
    The s_infty term for the S2KO version of the model

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jpm: double
            The PMCA flux

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
        Ke=par['K1']
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
        Ke=K1
    # Ke=par['Ke']
    n=par['n']
    w1=par['w1']
    a0 = par['a0']
    Jcrac = a0 + w1/(1+np.exp(n*(ce-Ke)))
    return Jcrac

def getS1KOCRAC(y, par): 
    """
    The s_infty term for the S1KO version of the model

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jpm: double
            The PMCA flux

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    n=par['n']
    w2=par['w2']
    a0 = par['a0']
    Ke=par['K2']
    Jcrac = a0 + w2/(1+np.exp(n*(ce-Ke)))
    return Jcrac

def getWTCRAC(y, par): 
    """
    The s_infty term for the WT version of the model
    The derivation for the expressions used in this functions is in PySym_STIMWTModel.py

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jpm: double
            The PMCA flux

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
        K1 = par['K1']
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    # Paramters
    K2 = par['K2']
    K12 = par['K12']
    K21 = par['K21']
    n = par['n']
    w1 = par['w1']
    w2 = par['w2']
    w12 = par['w12']
    w21 = par['w21']
    a0 = par['a0'] #a0 = 0 in the non-dynamic K1 version of the model. 
    phi1=np.exp(n*(-K1 + ce))
    phi2=np.exp(n*(-K2 + ce))
    phi3=np.exp(n*(-K1 + K21 + ce))
    phi4=np.exp(K12*n)
    # Active species
    S1a=(phi2 + 1)*(np.sqrt((phi1 + 1)*(phi2 + 1)*(phi1*phi2 + phi1 + phi2 + 4*phi3 + 4*phi4 + 1))/(2*(phi2*phi3 + phi2*phi4 + phi3 + phi4)) - (phi1 + 1)/(2*(phi3 + phi4)))/(phi1 + 1)
    S2a=np.sqrt((phi1 + 1)*(phi2 + 1)*(phi1*phi2 + phi1 + phi2 + 4*phi3 + 4*phi4 + 1))/(2*(phi2*phi3 + phi2*phi4 + phi3 + phi4)) - (phi1 + 1)/(2*(phi3 + phi4))
    S21=S1a*S2a*phi3
    S12=S1a*S2a*phi4
    # Flux
    Jcrac = a0 + w1*S1a + w2*S2a + w12*S12 + w21*S21
    return Jcrac

def getCDI(y, par): 
    """
    The i_infty term

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jcdi: double
            The s+infty term

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
        K1 = par['K1']
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    #Parameters
    Ki=par['Ki']
    si=par['si']
    #Value
    Jcdi =1/(1+np.exp(si*(cm-Ki)))
    return Jcdi

def getKcurve(y, par): 
    """
    The K_infty term

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    kCurve: double
            The K_infty term

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
        K1 = par['K1']
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    #Parameters
    Vk=par['Vk']
    kn=par['kn']
    nk=par['nk']
    k1Rest=par['k1Rest']
    kCurve = (Vk*kn**nk)/(kn**nk + ce**nk) + k1Rest
    return kCurve

def getTauK(y, par): 
    """
    The K1 timescale

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    TauK: double
            The K1 timescale

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
        K1 = par['K1']
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    #Parameters
    Tkmax=par['Tkmax']
    taun=par['taun']
    nk=par['nk']
    TauK=Tkmax*(ce**nk / (taun**nk + ce**nk))
    return TauK


def getDiff(y, par): 
    """
    The diffusion term

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    Jdiff: double
            The diffusion term

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
        K1 = par['K1']
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    #Parameters
    kdiff=par['kdiff']
    Jdiff=kdiff*(cm-c)
    return Jdiff

def getDegr(y, par): 
    """
    A combination of the Vplc production and Ca2+-based degradation terms for IP3

    Parameters
    ----------
    y: vector
            with the state variables
    par: dic
            parameters 

    Returns 
    ----------
    L: double
            The overall change in IP3 concentration

    """
    #Variables 
    if len(y) == 7: #model without dynamic K1 
        c, ce, cm, h, p, s, i = y
        K1 = par['K1']
    if len(y) == 8: #model with dynamic K1
        c, ce, cm, h, p, s, i, K1 = y
    #Parameters
    Vdeg=par['Vdeg']
    Kdeg=par['Kdeg']
    Vplc=par['Vplc']
    IPdeg=Vdeg*c**2/(c**2+Kdeg**2)
    L=Vplc - IPdeg*p
    return L
#%% Numerical integration functions
def getIVP(model, par, ics, tini=0, tf=100): 
    """
    Simple integration using the solve_ivp routine from Scipy 

    Parameters
    ----------
    model: func
            model to integrate
    par: dict 
            set of named parameters 
    ics: dict
            set of named initial conditions
    tini: double
            start of integration time
    tf: double
            end of integration time 

    Returns 
    ----------
    data: array, first entry is t value, rest of values are the model variables
    """
    data = solve_ivp(model, [tini, tini+tf],list(ics.values()), method = 'Radau', max_step = 500, rtol = 1e-7, atol = 1e-9, args = (par,))
    data = np.vstack((data.t, data.y))
    return data

def getnewICS(data): 
    """
    Get the initial conditions of the open cell model using the 
    last data point of a former run

    Parameters
    ----------
    data: array
            a numerical integration run, like the one generated by getIVP
    Returns 
    ----------
    ics: vector 
            the new initial conditions
    """
    ics = data[1:][:,-1]
    if len(ics) == 8: #Dynamic K1 
        ics = {
        'c': ics[0],
        'ce': ics[1],
        'cm': ics[2],
        'h': ics[3], 
        'p': ics[4],
        's': ics[5],
        'i': ics[6],
        'K1': ics[7]
        }
    if len(ics) == 7: #Non-dynmanic K1
        ics = { 
        'c': ics[0],
        'ce': ics[1],
        'cm': ics[2],
        'h': ics[3], 
        'p': ics[4],
        's': ics[5],
        'i': ics[6]
        }
    return ics

def thapsProtocol(par, model, ics, Tg, tini, tf, tthaps, tsoce, tguad, tfinal): 
    """
    The numerical integration simulating the maximal SOCE experimental assay (aka the Thapsigargin protocol), comprised of four numerical integrations 

    Zero run: integrate for tf-tini seconds at rest to find the steady state
    
    First run: forward time integration of our three models for tthaps seconds

    Second run: Remove SOCE, add thapsigargin, integrate for tsoce-tthaps seconds

    Third run: Add back SOCE, integrate for tguad-tsoce seconds

    Fourth run: Add Gd3+, integrate for tfinal-tguad

    Parameters
    ----------
    par: dic
            Parameter values
    model: func
            model to integrate 
    ics: dict 
            initial conditions
    Tg: double
            amount of SERCA blockage due to Tg treatment, 1=no effect, 0=transport disabled in the cyto->ER direction
    tini: double 
            starting time, usually t=0
    tf: double
            time of integration of the zero run used to find the steady state
    tthaps: double
            time of thapsigargin treatment and removal of extracellular Ca2+
    tsoce: double 
            time of readdition of extracellular Ca2+
    tguad: double 
            time of Gd3+ addition
    tfinal: double
            time at which the run ends

    Returns 
    ----------
    thapsRun: array
            The array contains the time series for all phase variables across the four runs. 

    """

    """
    Zero run
    """
    # Setting original parameters
    w1 = par['w1']
    w2 = par['w2']  
    w12 = par['w12']
    w21 = par['w21']
    # Finding resting state
    newICS = getIVP(model, par, ics, tini=tini, tf=tf)
    ics = getnewICS(newICS)
    """
    First run
    """
    # First run, starting at SS
    tini, tf = [0, tthaps]
    thaps0 = getIVP(model, par, ics, tini=tini, tf=tf)
    """
    Second run 
    """
    # Second run, Removing SOCE, adding thapsigargin
    ics = getnewICS(thaps0)
    par['Tg']=Tg
    par['w1']=0
    par['w2']=0
    par['w12']=0
    par['w21']=0
    tini, tf = [tthaps, tsoce-tthaps]
    thaps1 = getIVP(model, par, ics, tini=tini, tf=tf)
    """
    Third run
    """
    # Third run, adding SOCE back
    par['w1']=w1
    par['w2']=w2
    par['w12']=w12
    par['w21']=w21
    ics = getnewICS(thaps1)
    tini, tf = [tsoce, tguad-tsoce]
    thaps2 = getIVP(model, par, ics, tini=tini, tf=tf)
    """
    Fourth run 
    """
    # Fourth run, Adding Gd3+
    par['w1']=0
    par['w2']=0
    par['w12']=0
    par['w21']=0
    ics = getnewICS(thaps2)
    tini, tf = [tguad, tfinal-tguad]
    thaps3 = getIVP(model, par, ics, tini=tini, tf=tf)
    """
    Joining the four runs
    """
    thapsRun = np.concatenate([thaps0, thaps1, thaps2, thaps3], axis=1)
    return thapsRun

def cpaProtocol(par, model, ics, CPA, tini, tf, tthaps, tsoce, tfinal, VsDebuff): 
    """
    The numerical integration simulating the ER refilling assay (aka the CPA protocol), comprised of three numerical integrations 

    Zero run: integrate for tf-tini seconds at rest to find the steady state
    
    First run: forward time integration of our three models for tthaps seconds

    Second run: Remove SOCE, add CPA, integrate for tsoce-tthaps seconds

    Third run: Add back SOCE, remove CPA, integrate for tfinal-tsoce seconds

    Parameters
    ----------
    par: dic
            Parameter values
    model: func
            model to integrate 
    ics: dict 
            initial conditions
    CPA: double
            amount of SERCA blockage due to CPA treatment, 1=no effect, 0=transport disabled in the cyto->ER direction
    tini: double 
            starting time, usually t=0
    tf: double
            time of integration of the zero run used to find the steady state
    tthaps: double
            time of thapsigargin treatment and removal of extracellular Ca2+
    tsoce: double 
            time of readdition of extracellular Ca2+
    tfinal: double
            time at which the run ends
    VsDebuff: double 
            penalty applied to the SERCA pumps after removal of CPA. The experimental group told us that the cells post-CPA treatment were 'unhappy' and that they believe the SERCA pumps didnt recover properly after the CPA treatment. 1=no effect, 0=SERCA pumps disabled on both directions 

    Returns 
    ----------
    thapsRun: array
            The array contains the time series for all phase variables across the three runs. 
    """
    # Setting original parameters
    w1 = par['w1']
    w2 = par['w2']  
    w12 = par['w12']
    w21 = par['w21']
    # Finding resting state
    newICS = getIVP(model, par, ics, tini=tini, tf=tf)
    ics = getnewICS(newICS)
    # First run, starting at SS
    tini, tf = [0, tthaps]
    thaps0 = getIVP(model, par, ics, tini=tini, tf=tf)
    # Second run, Removing SOCE, adding thapsigargin
    ics = getnewICS(thaps0)
    par['Tg']=CPA
    par['w1']=0
    par['w2']=0
    par['w12']=0
    par['w21']=0
    tini, tf = [tthaps, tsoce-tthaps]
    thaps1 = getIVP(model, par, ics, tini=tini, tf=tf)
    # Third run, adding SOCE back, removing CPA, debuffing SERCA because of CPA treatment
    par['w1']=w1
    par['w2']=w2
    par['w12']=w12
    par['w21']=w21
    par['Tg']=0
    par['VsC']=par['VsC']*VsDebuff
    par['VsM']=par['VsM']*VsDebuff
    ics = getnewICS(thaps1)
    tini, tf = [tsoce, tfinal]
    thaps2 = getIVP(model, par, ics, tini=tini, tf=tf)
    thapsRun = np.concatenate([thaps0, thaps1, thaps2], axis=1)
    return thapsRun

#%% Base model samples

def sampleWTBase(par, Vplc, tf):
    """
    Integration of the base WT model after stimulation (represented by Vplc) for tf seconds. this run always uses the initial condition corresponding to the equiilibrium point at rest (Vplc = 0)

    Parameters
    ----------
    par: dict 
            set of named parameters 
    Vplc: double    
            amount of stimulation
    tf: double 
            total time of integration
        
    Returns 
    ----------
    data: array
            first entry is t value, rest of values are the model variables
    """
    ics = icsBaseJurkatCell()
    model = jurkatBaseWTCell
    tini=0
    #Resting State 
    ics['ce']=800
    newICS = getIVP(model, par, ics, tini=tini, tf=tf)
    ics = getnewICS(newICS)
    #Stimulation
    par['Vplc']=Vplc
    sample = getIVP(model, par, ics, tini=tini, tf=tf)
    return sample 

def sampleS2KOBase(par, Vplc, tf):
    """
    Integration of the base S2KO model after stimulation (represented by Vplc) for tf seconds. this run always uses the initial condition corresponding to the equiilibrium point at rest (Vplc = 0)

    Parameters
    ----------
    par: dict 
            set of named parameters 
    Vplc: double    
            amount of stimulation
    tf: double 
            total time of integration

    Returns 
    ----------
    data: array
            first entry is t value, rest of values are the model variables
    """
    ics = icsBaseJurkatCell()
    model = jurkatBaseS2KOCell
    tini=0
    par['Ke']=par['K1']
    par['Vsoce']=par['w1']
    #Resting State 
    ics['ce']=par['Ke']
    newICS = getIVP(model, par, ics, tini=tini, tf=tf)
    ics = getnewICS(newICS)
    #Stimulation
    par['Vplc']=Vplc
    sample = getIVP(model, par, ics, tini=tini, tf=tf)
    return sample 

def sampleS1KOBase(par, Vplc, tf):
    """
    Integration of the base S1KO model after stimulation (represented by Vplc) for tf seconds. this run always uses the initial condition corresponding to the equiilibrium point at rest (Vplc = 0)

    Parameters
    ----------
    par: dict 
            set of named parameters 
    Vplc: double    
            amount of stimulation
    tf: double 
            total time of integration

    Returns 
    ----------
    data: array
            first entry is t value, rest of values are the model variables
    """
    ics = icsBaseJurkatCell()
    model = jurkatBaseS1KOCell
    tini=0
    par['Ke']=par['K2']
    par['Vsoce']=par['w2']
    #Resting State 
    ics['ce']=par['Ke']
    newICS = getIVP(model, par, ics, tini=tini, tf=tf)
    ics = getnewICS(newICS)
    #Stimulation
    par['Vplc']=Vplc
    sample = getIVP(model, par, ics, tini=tini, tf=tf)
    return sample 
#%% Jurkat model samples
def sampleWT(par, Vplc, tf):
    """
    Integration of the base WT after stimulation (represented by Vplc) for tf seconds. this run always uses the initial condition corresponding to the equiilibrium point at rest (Vplc = 0)

    Parameters
    ----------
    par: dict 
            set of named parameters 
    Vplc: double    
            amount of stimulation
    tf: double 
            total time of integration

    Returns 
    ----------
    data: array
            first entry is t value, rest of values are the model variables
    """
    ics = icsJurkatCell()
    model = jurkatWTCell
    tini=0
    #Resting State 
    ics['ce']=800
    newICS = getIVP(model, par, ics, tini=tini, tf=tf)
    ics = getnewICS(newICS)
    #Stimulation
    par['Vplc']=Vplc
    sample = getIVP(model, par, ics, tini=tini, tf=tf)
    return sample 

def sampleS2KO(par, Vplc, tf):
    """
    Integration of the S2KO model after stimulation (represented by Vplc) for tf seconds. this run always uses the initial condition corresponding to the equiilibrium point at rest (Vplc = 0)

    Parameters
    ----------
    par: dict 
            set of named parameters 
    Vplc: double    
            amount of stimulation
    tf: double 
            total time of integration

    Returns 
    ----------
    data: array
            first entry is t value, rest of values are the model variables
    """
    ics = icsJurkatCell()
    model = jurkatS2KOCell
    tini=0
    par['Ke']=par['K1']
    par['Vsoce']=par['w1']
    #Resting State 
    ics['ce']=par['Ke']
    newICS = getIVP(model, par, ics, tini=tini, tf=tf)
    ics = getnewICS(newICS)
    #Stimulation
    par['Vplc']=Vplc
    sample = getIVP(model, par, ics, tini=tini, tf=tf)
    return sample 

def sampleS1KO(par, Vplc, tf):
    """
    Integration of the S1KO model after stimulation (represented by Vplc) for tf seconds. this run always uses the initial condition corresponding to the equiilibrium point at rest (Vplc = 0)

    Parameters
    ----------
    par: dict 
            set of named parameters 
    Vplc: double    
            amount of stimulation
    tf: double 
            total time of integration

    Returns 
    ----------
    data: array
            first entry is t value, rest of values are the model variables
    """
    ics = icsJurkatCell()
    model = jurkatS1KOCell
    tini=0
    par['Ke']=par['K2']
    par['Vsoce']=par['w2']
    #Resting State 
    ics['ce']=par['Ke']
    newICS = getIVP(model, par, ics, tini=tini, tf=tf)
    ics = getnewICS(newICS)
    #Stimulation
    par['Vplc']=Vplc
    sample = getIVP(model, par, ics, tini=tini, tf=tf)
    return sample 
#%% Stable Jurkat model samples
def sampleWTStable(par, Vplc, tf: int, tf2=500):
    """
    Integration of the WT model after stimulation (represented by Vplc) for tf seconds. 
    This functions does an initial integration for tf2 seconds to reach the stable periodic orbit firs. 

    Parameters
    ----------
    par: dict 
            set of named parameters 
    Vplc: double    
            amount of stimulation
    tf: double 
            total time of integration
    tf2: double
            Zero-run time of integation

    Returns 
    ----------
    data: array
            first entry is t value, rest of values are the model variables
    """
    ics = icsJurkatCell()
    model = jurkatWTCell
    tini=0
    #Resting State 
    ics['ce']=800
    newICS = getIVP(model, par, ics, tini=tini, tf=tf2)
    ics = getnewICS(newICS)
    #Stimulation
    par['Vplc']=Vplc
    newICS = getIVP(model, par, ics, tini=tini, tf=tf2)
    ics = getnewICS(newICS)
    sample = getIVP(model, par, ics, tini=tini, tf=tf)
    return sample 

def sampleS2KOStable(par, Vplc, tf, tf2=500):
    """
    Integration of the S2KO model after stimulation (represented by Vplc) for tf seconds. 
    This functions does an initial integration for tf2 seconds to reach the stable periodic orbit firs. 

    Parameters
    ----------
    par: dict 
            set of named parameters 
    Vplc: double    
            amount of stimulation
    tf: double 
            total time of integration
    tf2: double
            Zero-run time of integation

    Returns 
    ----------
    data: array
            first entry is t value, rest of values are the model variables
    """
    ics = icsJurkatCell()
    model = jurkatS2KOCell        
    tini=0
    #These next lines ensure that the basic CRAC function uses the STIM1 parameters
    #this might not be necessary because I defined the S2KOCRAC separately, but I guess it can't hurt to leave it here
    par['Ke']=par['K1']
    par['Vsoce']=par['w1']
    #Resting State 
    ics['ce']=par['Ke']
    newICS = getIVP(model, par, ics, tini=tini, tf=tf2)
    ics = getnewICS(newICS)
    #Stimulation
    par['Vplc']=Vplc
    newICS = getIVP(model, par, ics, tini=tini, tf=tf2)
    ics = getnewICS(newICS)
    sample = getIVP(model, par, ics, tini=tini, tf=tf)
    return sample 

def sampleS1KOStable(par, Vplc, tf, tf2=500):
    """
    Integration of the S1KO model after stimulation (represented by Vplc) for tf seconds. 
    This functions does an initial integration for tf2 seconds to reach the stable periodic orbit firs. 

    Parameters
    ----------
    par: dict 
            set of named parameters 
    Vplc: double    
            amount of stimulation
    tf: double 
            total time of integration
    tf2: double
            Zero-run time of integation

    Returns 
    ----------
    data: array
            first entry is t value, rest of values are the model variables
    """
    ics = icsJurkatCell()
    model = jurkatS1KOCell
    tini=0
    #These next lines ensure that the basic CRAC function uses the STIM2 parameters
    #this might not be necessary because I defined the S1KOCRAC separately, but I guess it can't hurt to leave it here
    par['Ke']=par['K2']
    par['Vsoce']=par['w2']
    #Resting State 
    ics['ce']=par['Ke']
    newICS = getIVP(model, par, ics, tini=tini, tf=tf2)
    ics = getnewICS(newICS)
    #Stimulation
    par['Vplc']=Vplc
    newICS = getIVP(model, par, ics, tini=tini, tf=tf2)
    ics = getnewICS(newICS)
    sample = getIVP(model, par, ics, tini=tini, tf=tf)
    return sample 

#%% Time series plotting functions
def plotTS(ax, xtrace, ytrace, xlim, ylim,color, xticks=[], yticks=[],  label='',  xlabel='', ylabel='', subtitle='', legend=False, ylabelpad = -10, legendsize=8):
    """
    Plots one individual time series in a given axes
    I think there is a more elegant way to pass all these extra arguments to matlpotlib
    But I have not found it yet

    Parameters
    ----------
    ax: axes
            The axes to plot in
    xtrace: array
            usually the time variable 
    ytrace: array
            usually a time series 
    other: mutliple types
            extra plotting arguments for matplotlib

    Returns 
    ----------
    none
    """
    ax.plot(xtrace, ytrace, color = color, label = label)
    # Ax limits 
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_ylabel(ylabel, rotation='horizontal', labelpad=ylabelpad, ma='center')
    ax.set_xlabel(xlabel, labelpad=0)
    ax.set_title(subtitle, loc='left')
    if legend:
        ax.legend(prop={'size': legendsize})
    return 0 

def plotDataSample(df, colname, ax, color, title, ylim, xticks=[], yticks=[], xlabel='', ylabel = ''):
    """
    Plots one individual experimental time series in a given axes
    The experimental time series are stored in a pandas dataframe instead
    I think there is a more elegant way to pass all these extra arguments to matlpotlib
    But I have not found it yet

    Parameters
    ----------
    df: dataframe 
            Dataframe with the experimental time series
    colname: str
            Name of the sample to plot
    ax: axes
            The axes to plot in
    other: mutliple types
            extra plotting arguments for matplotlib

    Returns 
    ----------
    none
    """
    df.plot(y=colname, ax=ax, color=color, legend=False)
    ax.set_xlim([0, 450])
    ax.set_xticks(xticks)
    ax.set_xlabel(xlabel)
    ax.set_ylim(ylim)
    ax.set_yticks(yticks)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc='left')
    return 0 

def plotAUTOBranch(branch, xVar, yVar, ax, color = 'k'):
    eqStab = branch.stability()
    eqStab.insert(0,0)
    for i in range(1,len(eqStab)):
        if eqStab[i]< 0: #Stable curve is plotted as solid line
            ax.plot(
                branch[max(np.abs(eqStab[i-1])-1,0):np.abs(eqStab[i])][xVar], 
                branch[max(np.abs(eqStab[i-1])-1,0):np.abs(eqStab[i])][yVar], 
                color, 
                ls='solid', 
                alpha=1, 
                zorder=-1,
            )
        else:
            ax.plot( #unstable curve is plotted as dashed line
                branch[max(np.abs(eqStab[i-1])-1,0):np.abs(eqStab[i])][xVar], 
                branch[max(np.abs(eqStab[i-1])-1,0):np.abs(eqStab[i])][yVar], 
                color,
                ls='dashed',
                alpha=1,
                zorder=-1,
            )

#%% AUTO f90 writing files
"""
Miscellaneous functions designed to facilitate writing the AUTO f90 files 
"""
"""
getICS - get initial conditions in a printed format for the specified parameters
getParNum - return format gammaE=PAR(1) - f90 file
getParValues return format gammaE=1 - f90 file
getNumPar - return format PAR(1)=gammaE - f90 file
getParnames - return format "gammaE": 1 - constants file
"""

def printICS(ics): 
    """
    initial conditions in the f90 format

    Parameters 
    ----------
    ics: dict 
            set of initial conditions 

    Returns
    ----------
    none
    """
    i=1 
    for key, value in ics.items():
        print(r'U(' + str(i) + r')' + r'=' + str(value))
        i=i+1
    
def getParNum(par):
    """
    Prints parameters in the format parName=AUTOparName (gammaE=PAR(1))

    Parameters 
    ----------
    par: dict 
            set of parameters 

    Returns
    ----------
    none
    """
    i=1 
    for key, value in par.items():
        if type(value) != str: 
            print(str(key) + r'=PAR(' + str(i) + r')' )
            i=i+1
            if i==11:
                print(r'PERIOD=PAR(11)')
                i=i+1
            if i==14:
                print(r'TIME=PAR(14)')
                i=i+1

def getParValues(par):
    """
    Prints parameters in the format parName=parValue (gammaE=1)

    Parameters 
    ----------
    par: dict 
            set of parameters 

    Returns
    ----------
    none
    """
    i=1 
    for key, value in par.items():
        if type(value) != str: 
            print(str(key) + r'=' + str(value) + r'' )
            i=i+1
            if i==11:
                i=i+1
            if i==14:
                i=i+1

def getNumPar(par):
    """
    Prints parameters in the format AUTOparName=parName (PAR(1)=gammaE)

    Parameters 
    ----------
    par: dict 
            set of parameters 

    Returns
    ----------
    none
    """
    i=1 
    for key, value in par.items():
        if type(value) != str: 
            print(r'PAR(' + str(i) + r')=' + str(key) )
            i=i+1
            if i==11:
                print(r'PAR(11)=0.0')
                i=i+1
            if i==14:
                print(r'PAR(14)=0.0')
                i=i+1

def getParnames(par): 
    """
    Prints parameters in the format parName: parValue ("gammaE": 1 )

    Parameters 
    ----------
    par: dict 
            set of parameters 

    Returns
    ----------
    none
    """

    i=1
    for key, value in par.items():
        if type(value) != str: 
            print("" + str(i) + ": '" + key +"'," )
            i+=1
            if i==11:
                print("" + r'11' + ": '" + r'PERIOD' +"'," )
                i=i+1
            if i==14:
                print("" + r'14' + ": '" + r'TIME' +"'," )
                i=i+1

# %%
