#%% Libraries
"""
Ccontains some utility functions to deal with AUTO
"""
import os
from pathlib import Path
from auto.AUTOCommands import run, load, save, clean, klb, loadbd


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
        os.chdir(mainFolder)
        return bd
    except:
        print ('Bifurcation diagram files not found in AUTO folder')
    