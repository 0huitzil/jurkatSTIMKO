#%%
"""
Contains functions which generate the bifurcation diagram files for the WT model
"""
from auto.AUTOCommands import run, load, save, clean

file = "AUTOJurkatWTCell"
cts = 'AUTOJurkatCell'
model = load(file, constants = cts) 
Vplc = [0.1, 0.2, 0.3]
eqWT = run(
    model, 
    IPS=1, 
    ICP=['Vplc'], 
    NMX=20000,
    DS=1e-2,
    DSMAX=5e-2,
    UZSTOP={'Vplc': 2}, 
    UZR={'Vplc': Vplc}
)
Vplc = [0.1, 0.2, 0.3]
cycleWT = run(
    eqWT('HB'), 
    IPS=2,
    # ISP=2, 
    ICP=['Vplc', 11, 'MIN c'], 
    NMX=5000,
    NTST=500,
    DS=1e-2,
    DSMAX=1e-0,
    SP=['LP2', 'UZ', 'PD2',],
    UZSTOP={'Vplc': 2}, 
    UZR={'Vplc': Vplc}
)
save(eqWT+cycleWT, 'WT')
clean()