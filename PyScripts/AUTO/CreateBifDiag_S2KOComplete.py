#%%
"""
Contains functions which generate the bifurcation diagram files for the S2KO model
"""
from auto.AUTOCommands import run, load, save, clean

file = "AUTOJurkatS2Cell"
model = load(file) 
eqS2 = run(
    model, 
    IPS=1, 
    ICP=['Vplc'], 
    NMX=20000,
    DS=1e-2,
    DSMAX=5e-2,
    UZSTOP={'Vplc': 0.5}
)
#%
# Vplc = [0.06, 0.1, 0.2, 0.35]
cycleS2 = run(
    eqS2('HB')[0], 
    IPS=2,
    JAC=1,
    ICP=['Vplc', 11], 
    NMX=150000,
    NTST=1000,
    DS=1e-3,
    DSMAX=3e-2,
    DSMIN=1e-6,
    SP=['LP0', 'UZ', 'PD14', 'TR0'],
    UZSTOP={'Vplc': 1.76}, 
    # UZR={'Vplc': Vplc}, 
)
save(eqS2+cycleS2, 'S2a') #done up to this point 
cycleS2c = run(
    eqS2('HB')[2], 
    IPS=2,
    JAC=1,
    ICP=['Vplc', 11], 
    NMX=10000,
    NTST=1000,
    DS=1e-3,
    DSMAX=5e-2,
    DSMIN=1e-6,
    SP=['LP2', 'UZ', 'PD10', 'TR0'],
    UZSTOP={'Vplc': 2}, 
    # UZR={'Vplc': Vplc}, 
)
save(eqS2+cycleS2+cycleS2c, 'S2a')
clean()