#%%
"""
Contains functions which generate the bifurcation diagram files for the S1KO model
"""
from auto.AUTOCommands import run, load, save, clean

file = "AUTOJurkatS1Cell"
model = load(file) 
eqS1 = run(
    model, 
    IPS=1, 
    ICP=['Vplc'], 
    NMX=20000,
    DS=1e-2,
    DSMAX=5e-2,
    UZSTOP={'Vplc': 2}
)
#%
# Vplc = [0.06, 0.1, 0.2, 0.35]
cycleS1 = run(
    eqS1('HB')[0], 
    IPS=2,
    JAC=1,
    ICP=['Vplc', 11], 
    NMX=150000,
    NTST=1000,
    DS=1e-3,
    DSMAX=3e-2,
    DSMIN=1e-6,
    SP=['LP0', 'UZ', 'PD', 'TR0'],
    UZSTOP={'Vplc': 1.76}, 
    # UZR={'Vplc': Vplc}, 
)
save(eqS1+cycleS1, 'S1a')
cycleS1b = run(
    eqS1('HB')[1], 
    IPS=2,
    JAC=1,
    ICP=['Vplc', 11], 
    NMX=80000,
    NTST=1000,
    DS=1e-3,
    DSMAX=3e-2,
    DSMIN=1e-6,
    SP=['LP100', 'UZ', 'PD100', 'TR0'],
    UZSTOP={'Vplc': 2}, 
    # UZR={'Vplc': Vplc}, 
)
save(eqS1+cycleS1+cycleS1b, 'S1a')
cycleS1c = run(
    eqS1('HB')[2], 
    IPS=2,
    JAC=1,
    ICP=['Vplc', 11], 
    NMX=40000,
    NTST=1000,
    DS=1e-3,
    DSMAX=3e-2,
    DSMIN=1e-6,
    SP=['LP0', 'UZ', 'PD0', 'TR0'],
    UZSTOP={'Vplc': 2}, 
    # UZR={'Vplc': Vplc}, 
)
save(eqS1+cycleS1+cycleS1b+cycleS1c, 'S1a')
clean()