''' Interpolation Routines '''
import math
from sharppy.sharptab import vector
from sharppy.sharptab import thermo
from sharppy.sharptab.constants import *



__all__ = ['i_pres', 'i_hght', 'i_temp', 'i_dwpt', 'i_vec', 'i_vtmp',
           'interp_from_pres', 'interp_from_hght', 'agl', 'msl']


def pres(h, prof):
    '''
    Interpolates the given data to calculate a pressure at a given height

    Inputs
    ------
        h           (float)                 Height (m) of a level
        prof        (profile object)        Profile object

    Returns
    -------
        Pressure (hPa) at the given height
    '''
    return interp_from_hght(h, prof, prof.pind)


def hght(p, prof):
    '''
    Interpolates the given data to calculate a height at a given pressure

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        prof        (profile object)        Profile object

    Returns
    -------
        Height (m) at the given pressure
    '''
    return interp_from_pres(p, prof, prof.zind)


def temp(p, prof):
    '''
    Interpolates the given data to calculate a temperature at a given pressure

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        prof        (profile object)        Profile object

    Returns
    -------
        Temperature (C) at the given pressure
    '''
    return interp_from_pres(p, prof, prof.tind)


def dwpt(p, prof):
    '''
    Interpolates the given data to calculate a dew point at a given pressure

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        prof        (profile object)        Profile object

    Returns
    -------
        Dew point (C) at the given pressure
    '''
    return interp_from_pres(p, prof, prof.tdind)


def vtmp(p, prof):
    '''
    Interpolates the given data to calculate a virtual temperature at a
    given pressure

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        prof        (profile object)        Profile object

    Returns
    -------
        Virtual temperature (C) at the given pressure
    '''
    t = temp(p, prof)
    td = dwpt(p, prof)
    return thermo.virtemp(p, t, td)


def components(p, prof):
    '''
    Interpolates the given data to calculate the U and V components at a
    given pressure

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        prof        (profile object)        Profile object

    Returns
    -------
        U/V components at the given pressure
    '''
    U = interp_from_pres(p, prof, prof.uind)
    V = interp_from_pres(p, prof, prof.vind)
    return U, V


def vec(p, prof):
    '''
    Interpolates the given component data to a given pressure level
    and returns the interpolated direction and speed

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        prof        (profile object)        Profile object

    Returns
    -------
        dir         (float)                 Direction
        mag         (float)                 Magnitude
    '''
    u = interp_from_pres(p, prof, prof.uind)
    v = interp_from_pres(p, prof, prof.vind)
    return vector.comp2vec(u, v)



def interp_from_hght(h, prof, ind):
    '''
    General interpolation routine for height coordinates

    Inputs
    ------
        h           (float)                 Height (m) of a level
        prof        (profile object)        Profile object
        ind         (integer)               Index of variable to interpolate

    Returns
    -------
        Interpolated variable
    '''
    if not QC(h): return RMISSD
    tptr = -1
    bptr = 0
    for i in range(prof.gNumLevels):
        if not QC(prof.gSndg[i][ind]): continue
        if prof.gSndg[i][prof.zind] < h:
            bptr = i
        elif prof.gSndg[i][prof.zind] == h:
            return prof.gSndg[i][ind]
        else:
            tptr = i
            nm1 = h - prof.gSndg[bptr][prof.zind]
            nm2 = prof.gSndg[tptr][prof.zind] - prof.gSndg[bptr][prof.zind]
            val = prof.gSndg[bptr][ind] + ((nm1 / nm2) * \
                (prof.gSndg[tptr][ind] - prof.gSndg[bptr][ind]))
            return val


def interp_from_pres(p, prof, ind):
    '''
    General interpolation routine for pressure coordinates

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        prof        (profile object)        Profile object
        ind         (integer)               Index of variable to interpolate

    Returns
    -------
        Interpolated variable
    '''
    if not QC(p): return RMISSD
    tptr = -1
    bptr = 0
    for i in range(prof.gNumLevels):
        if not QC(prof.gSndg[i][ind]): continue
        if prof.gSndg[i][prof.pind] > p:
            bptr = i
        elif prof.gSndg[i][prof.pind] == p:
            return prof.gSndg[i][ind]
        else:
            tptr = i
            nm1 = prof.gSndg[tptr][ind] - prof.gSndg[bptr][ind]
            nm2 = math.log(prof.gSndg[bptr][prof.pind] / \
                           prof.gSndg[tptr][prof.pind])
            nm3 = math.log(prof.gSndg[bptr][prof.pind] / p)
            return prof.gSndg[bptr][ind] + ((nm3 / nm2) * nm1)


def agl(h, prof):
    '''
    Convert a height from mean sea-level (MSL) to above ground-level (AGL)

    Inputs
    ------
        h           (float)                 Height of a level
        prof        (profile object)        Profile object

    Returns
    -------
        Converted height
    '''
    if not QC(h): return RMISSD
    return h - prof.gSndg[prof.sfc][prof.zind]


def msl(h, prof):
    '''
    Convert a height from above ground-level (AGL) to mean sea-level (MSL)

    Inputs
    ------
        h           (float)                 Height of a level
        prof        (profile object)        Profile object

    Returns
    -------
        Converted height
    '''
    if not QC(h): return RMISSD
    return h + prof.gSndg[prof.sfc][prof.zind]


