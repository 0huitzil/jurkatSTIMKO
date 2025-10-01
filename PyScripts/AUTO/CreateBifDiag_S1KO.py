#%%
"""
Contains functions which generate the bifurcation diagram files for the S1KO model
"""
from auto.AUTOCommands import run, load, save, clean

file = "AUTOJurkatS1KOCell"
cts = 'AUTOJurkatCell'
model = load(file, constants = cts) 
Vplc = [0.044, 0.07, 0.1, 0.16]
eqS1KO = run(
    model, 
    IPS=1, 
    ICP=['Vplc'], 
    NMX=20000,
    DS=1e-2,
    DSMAX=5e-2,
    UZSTOP={'Vplc': 2}
)
cycleS1KO = run(
    eqS1KO('HB')[0], 
    IPS=2,
    JAC=1,
    ICP=['Vplc', 11], 
    NMX=140000,
    NPR=5000,
    NTST=1200,
    DS=1e-3,
    DSMAX=5e-2,
    DSMIN=1e-6,
    SP=['LP0', 'UZ7', 'PD0', 'TR0'],
    UZSTOP={'Vplc': 1}, 
    UZR={'Vplc': Vplc}, 
)
save(eqS1KO+cycleS1KO, 'S1KO')
clean()