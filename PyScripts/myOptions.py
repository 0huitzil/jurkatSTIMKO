import matplotlib
import numpy as np
from cycler import cycler

# auto_directory = "/home/huitzil/auto/07p/python" 

def myRcParams():
    """
    Personal chanes to the matplotlib RcParams matplotlib file
    """
    options = {
    # 'backend': 'PS',
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
     'font.serif': ['Computer Modern Roman'],
    'text.usetex': True,
    'axes.grid': False, 
    # 'axes.spines.left' : False, 
    'axes.spines.right' : False, 
    'axes.spines.top' : False, 
    # 'axes.spines.bottom' : False, 
    'lines.linewidth' : 1,
    'font.size': 8, 
    "axes.prop_cycle":cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']), #Th default matplotlib colors, AUTO somehow changes it? 
    'mathtext.default':  'regular' 
    }
    return options


def set_figsize(fraction=1, subplots=(1, 1), export=False, width = 370.0):
    """
    Set figure dimensions to sit nicely in our document.

    Parameters
    ----------
    width_pt: float
            Document width in points
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    export: Boolean
            Determines is the fig should be sized to fit inside a latex document (True), 
            or sized in a more appropiate way to be seen in python (False)
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    if export:
        fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])
    else:
        fig_height_in = fig_width_in * golden_ratio

    return (fig_width_in, fig_height_in)
myATol = 1e-12
myRTol = 1e-6

# def calculate_xticks(ax, ticks, round_to=0.1, center=False):
#     """ 
#     I made a function based on the solution above that makes sure that the labels don't end up to be something with a lot of decimals:
#     From: https://stackoverflow.com/questions/20243683/matplotlib-align-twinx-tick-marks
#     """
#     upperbound = np.ceil(ax.get_xbound()[1]/round_to)
#     lowerbound = np.floor(ax.get_xbound()[0]/round_to)
#     dx = upperbound - lowerbound
#     fit = np.floor(dx/(ticks - 1)) + 1
#     dx_new = (ticks - 1)*fit
#     if center:
#         offset = np.floor((dx_new - dx)/2)
#         lowerbound = lowerbound - offset
#     values = np.linspace(lowerbound, lowerbound + dx_new, ticks)
#     values = values*round_to
#     return values*round_to

def calculate_xticks(ax, ticks, round_to=0.1, center=False):
    """ 
    I made a function based on the solution above that makes sure that the labels don't end up to be something with a lot of decimals:
    From: https://stackoverflow.com/questions/20243683/matplotlib-align-twinx-tick-marks
    """
    upperbound = np.ceil(ax.get_xbound()[1]/round_to)
    lowerbound = np.floor(ax.get_xbound()[0]/round_to)
    dx = upperbound - lowerbound
    fit = np.floor(dx/(ticks - 1)) + 1
    dx_new = (ticks - 1)*fit
    if center:
        offset = np.floor((dx_new - dx)/2)
        lowerbound = lowerbound - offset
    values = np.linspace(lowerbound, lowerbound + dx_new, ticks)
    values = values*round_to
    return values*round_to

def calculate_yticks(ax, ticks, round_to=0.1, center=False):
    """ 
    I made a function based on the solution above that makes sure that the labels don't end up to be something with a lot of decimals:
    From: https://stackoverflow.com/questions/20243683/matplotlib-align-twinx-tick-marks
    """
    upperbound = np.ceil(ax.get_ybound()[1]/round_to)
    lowerbound = np.floor(ax.get_ybound()[0]/round_to)
    dy = upperbound - lowerbound
    fit = np.floor(dy/(ticks - 1)) + 1
    dy_new = (ticks - 1)*fit
    if center:
        offset = np.floor((dy_new - dy)/2)
        lowerbound = lowerbound - offset
    values = np.linspace(lowerbound, lowerbound + dy_new, ticks)
    return values*round_to