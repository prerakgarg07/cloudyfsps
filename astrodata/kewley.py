import numpy as np
import matplotlib.pyplot as plt
import cubehelix
import matplotlib.colors as mpl_colors
import matplotlib.cm as cmx

def NII_OIII_sf_lim(ax=None, color='k', **kwargs):
    '''
    Kewley (2006) division for SF/composite
    [NII]6584/Ha vs [OIII]5007/Hb
    '''
    if ax is None:
        ax = plt.gca()
    NII = np.linspace(-1.3, 0.0)
    OIII = 0.61/(NII - 0.05) + 1.3
    ax.plot(NII, OIII, '--', color=color, lw=2)
    return
def NII_OIII_agn_lim(ax=None, color='k', **kwargs):
    '''
    Kewley (2006) division for composite/agn
    [NII]6584/Ha vs [OIII]5007/Hb
    '''
    if ax is None:
        ax = plt.gca()
    NII = np.linspace(-2.0, 0.35)
    OIII = 0.61/(NII - 0.47) + 1.19
    ax.plot(NII, OIII, color=color, lw=2)
    return
def SII_OIII_agn_lim(ax=None, color='k', **kwargs):
    '''
    Kewley (2006)
    [SII]6717,31/Ha vs [OIII]5007/Hb
    '''
    if ax is None:
        ax = plt.gca()
    SII = np.linspace(-2.0, 0.5)
    OIII = 0.72/(SII - 0.32) + 1.3
    ax.plot(SII, OIII, color=color, lw=2)
    return
def OI_OIII_agn(ax=None, color='k', **kwargs):
    '''
    Kewley (2006)
    [OI]6300/Ha vs [OIII]5007/Hb
    '''
    if ax is None:
        ax = plt.gca()
    OI = np.linspace(-2.0, 0.0)
    OIII = 0.73/(OI + 0.59) + 1.33
    ax.plot(OI, OIII, color=color, lw=2)
    return
