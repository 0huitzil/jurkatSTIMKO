#%% Libraries
"""
Contains functions which generate the bifurcation diagram files for the S1KO, S2KO and WT models, in case they do not exist in the directory already
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from auto.AUTOCommands import run, load, save, clean, klb, loadbd
from pathlib import Path
import os
"""
This command allows me to export the files directly to the 
sister LaTeX directory. Feel free to comment
"""
from myOptions import latexPath
"""
importing changes to the RcParams
"""
# from myOptions import *
# from PyModels import parJurkatCell, sampleS1Stable, plotTS
# matplotlib.rcParams.update(myRcParams())

#%% Data Collection
"""
Data Vplc Diag 
Only run this is if the .S1a files have not been created yet
"""
# def getWTBifDiag():
#     file = "AUTOJurkatWTCell"
#     cts = 'AUTOJurkatCell'
#     model = load(file, constants = cts) 
#     Vplc = [0.1, 0.2, 0.3]
#     eqWT = run(
#         model, 
#         IPS=1, 
#         ICP=['Vplc'], 
#         NMX=20000,
#         DS=1e-2,
#         DSMAX=5e-2,
#         UZSTOP={'Vplc': 2}, 
#         UZR={'Vplc': Vplc}
#     )
#     Vplc = [0.1, 0.2, 0.3]
#     cycleWT = run(
#         eqWT('HB'), 
#         IPS=2,
#         # ISP=2, 
#         ICP=['Vplc', 11, 'MIN c'], 
#         NMX=5000,
#         NTST=500,
#         DS=1e-2,
#         DSMAX=1e-0,
#         SP=['LP2', 'UZ', 'PD2',],
#         UZSTOP={'Vplc': 2}, 
#         UZR={'Vplc': Vplc}
#     )
#     save(eqWT+cycleWT, 'WT')
#     clean()

# file = "AUTOJurkatS1Cell"
# model = load(file) 
# eqS1 = run(
#     model, 
#     IPS=1, 
#     ICP=['Vplc'], 
#     NMX=20000,
#     DS=1e-2,
#     DSMAX=5e-2,
#     UZSTOP={'Vplc': 2}
# )
# #%
# # Vplc = [0.06, 0.1, 0.2, 0.35]
# cycleS1 = run(
#     eqS1('HB')[0], 
#     IPS=2,
#     JAC=1,
#     ICP=['Vplc', 11], 
#     NMX=150000,
#     NTST=1000,
#     DS=1e-3,
#     DSMAX=3e-2,
#     DSMIN=1e-6,
#     SP=['LP0', 'UZ', 'PD', 'TR0'],
#     UZSTOP={'Vplc': 1.76}, 
#     # UZR={'Vplc': Vplc}, 
# )
# save(eqS1+cycleS1, 'S1a')
# cycleS1b = run(
#     eqS1('HB')[1], 
#     IPS=2,
#     JAC=1,
#     ICP=['Vplc', 11], 
#     NMX=80000,
#     NTST=1000,
#     DS=1e-3,
#     DSMAX=3e-2,
#     DSMIN=1e-6,
#     SP=['LP100', 'UZ', 'PD100', 'TR0'],
#     UZSTOP={'Vplc': 2}, 
#     # UZR={'Vplc': Vplc}, 
# )
# save(eqS1+cycleS1+cycleS1b, 'S1a')
# cycleS1c = run(
#     eqS1('HB')[2], 
#     IPS=2,
#     JAC=1,
#     ICP=['Vplc', 11], 
#     NMX=40000,
#     NTST=1000,
#     DS=1e-3,
#     DSMAX=3e-2,
#     DSMIN=1e-6,
#     SP=['LP0', 'UZ', 'PD0', 'TR0'],
#     UZSTOP={'Vplc': 2}, 
#     # UZR={'Vplc': Vplc}, 
# )
# save(eqS1+cycleS1+cycleS1b+cycleS1c, 'S1a')
# cl()

#%% 

def loadBifDiag(name=str):
    """
    Load a bifurcation diagram object from a subfolder.
    Must move into the subfolder, load the object and move back into the original directory
    """
    try:
        mainFolder = Path.cwd()
        autoFolder = mainFolder / 'AUTO'
        os.chdir(autoFolder)
        bd = loadbd(name)
        # os.chdir(mainFolder)
        return bd
    except:
        print ('Bifurcation diagram files not found in AUTO folder')
    