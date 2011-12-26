''' Wind Manipulation Routines '''
import math
from sharppy.sharptab import interp
from sharppy.sharptab.constants import *

__all__ = ['mean_wind']


def mean_wind(pbot, ptop, profile, psteps=20):
    '''
    Calculates a pressure-weighted mean wint through a layer. The default
    layer is 850 to 200 hPa.

    Inputs
    ------
        pbot    (float)             Pressure of the bottom level (hPa)
        ptop    (float)             Pressure of the top level (hPa)
        profile (profile object)    Profile Object
        psteps  (int)               Number of steps to loop through (int)

    Returns
    -------
        mnu      (float)            U-component
        mnv      (float)            V-component
    '''
    if pbot == -1: lower = 850.
    if ptop == -1: upper = 200.
    pinc = int((pbot - ptop) / psteps)
    if pinc < 1:
        u1 = interp.interp_from_pres(pbot, profile, 4) * pbot
        u2 = interp.interp_from_pres(pbot, profile, 4) * ptop
        v1 = interp.interp_from_pres(pbot, profile, 5) * pbot
        v2 = interp.interp_from_pres(pbot, profile, 5) * ptop
        usum = u1 + u2
        vsum = v1 + v2
        wgt = pbot + ptop
    else:
        wgt = 0
        usum = 0
        vsum = 0
        for p in range(int(pbot), int(ptop), -pinc):
            usum += interp.interp_from_pres(p, profile, 4) * p
            vsum += interp.interp_from_pres(p, profile, 5) * p
            wgt += p

    return float(usum / wgt), float(vsum / wgt)
