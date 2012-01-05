''' Library for Composite Indices '''
import math
from sharppy.sharptab import *
from sharppy.sharptab.constants import *


__all__ = ['scp']


def scp(prof, **kwargs):
    '''
    Computes the Supercell Composite Parameter; uses the effective layer

    Inputs
    ------
        prof    (profile object)        Profile Object

    Returns
    -------
        supercell composite parameter   (float)
    '''

    # If MUPCL provided, use it, otherwise create MUPCL
    if 'mupcl' in kwargs:
        mupcl = kwargs.get('mupcl')
    else:
        mulplvals = params.DefineParcel(3, prof, pres=400)
        mupcl = params.parcelx(-1, -1, mulplvals.pres, mulplvals.temp,
            mulplvals.dwpt, prof, mulplvals.flag)

    # If EPCL provided, use it, otherwise create EPCL
    if 'epcl' in kwargs:
        epcl = kwargs.get('epcl')
    else:
        elplvals = params.DefineParcel(6, prof)
        epcl = params.parcelx(-1, -1, elplvals.pres, elplvals.temp,
            elplvals.dwpt, prof, lplvals=elplvals)

    cape = mupcl.bplus
    elhght = mupcl.elhght
    if not QC(cape): return RMISSD
    rstu, rstv = params.bunkers_storm_motion(prof, mupcl=mupcl)[:2]
    if cape < 100. or elhght < 0.:
        base = 0.
        ztop = 6000.
        pbot = interp.pres(interp.msl(base, prof), prof)
        ptop = interp.pres(interp.msl(ztop, prof), prof)
        shru, shrv = winds.wind_shear(pbot, ptop, prof)
        shrmag = vector.mag(shru, shrv)
        base = 0
        ztop = 6000.
    else:
        base = interp.agl(interp.hght(epcl.lplvals.pbot, prof), prof)
        depth = elhght - base
        top = base + (depth / 2.)
        ptmp = interp.pres(top, prof)
        shru, shrv = winds.wind_shear(epcl.lplvals.pbot, ptmp, prof)
        shrmag = vector.mag(shru, shrv)
        ztop = interp.agl(interp.hght(epcl.lplvals.ptop, prof), prof)

    esrh = winds.helicity(base, ztop, prof, rstu, rstv)[0]

    if not QC(esrh): return RMISSD
    if shrmag < 20: eshear = 0.
    elif shrmag > 40: eshear = 1.
    else: eshear = shrmag / 40.

    scp = eshear * (esrh / 50.) * (cape / 1000.)
    if scp <= 0: scp = 0
    return scp

