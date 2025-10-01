#%%
"""
Contains functions which generate the bifurcation diagram files for the S2KO model
"""
from auto.AUTOCommands import run, load, save, clean

file = "AUTOJurkatS2KOCell"
cts = 'AUTOJurkatCell'
model = load(file, constants = cts) 
eqS2KO = run(
    model, 
    IPS=1, 
    ICP=['Vplc'], 
    NMX=20000,
    DS=1e-2,
    DSMAX=5e-2,
    UZSTOP={'Vplc': 2}
)
#%
Vplc = [0.06, 0.1, 0.2, 0.35]
cycleS2KO = run(
    eqS2KO('HB')[0], 
    IPS=2,
    JAC=1,
    ICP=['Vplc', 11], 
    NMX=150000,
    NTST=1000,
    DS=1e-3,
    DSMAX=1e-1,
    DSMIN=1e-6,
    SP=['LP0', 'UZ', 'PD0', 'TR0'],
    UZSTOP={'Vplc': 1}, 
    UZR={'Vplc': Vplc}, 
)
save(eqS2KO+cycleS2KO, 'S2KO')
clean()