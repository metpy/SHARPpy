''' Interpolation Routines '''
import math
from sharppy.sharptab import vector
from sharppy.sharptab import thermo
from sharppy.sharptab.constants import *



__all__ = ['i_pres', 'i_hght', 'i_temp', 'i_dwpt', 'i_vec', 'i_vtmp',
           'interp_from_pres', 'interp_from_hght', 'agl', 'msl']


def pres(h, profile):
    '''
    Interpolates the given data to calculate a pressure at a given height

    Inputs
    ------
        h           (float)                 Height (m) of a level
        profile     (profile object)        Profile object

    Returns
    -------
        Pressure (hPa) at the given height
    '''
    return interp_from_hght(h, profile, 0)


def hght(p, profile):
    '''
    Interpolates the given data to calculate a height at a given pressure

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        profile     (profile object)        Profile object

    Returns
    -------
        Height (m) at the given pressure
    '''
    return interp_from_pres(p, profile, 1)


def temp(p, profile):
    '''
    Interpolates the given data to calculate a temperature at a given pressure

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        profile     (profile object)        Profile object

    Returns
    -------
        Temperature (C) at the given pressure
    '''
    return interp_from_pres(p, profile, 2)


def dwpt(p, profile):
    '''
    Interpolates the given data to calculate a dew point at a given pressure

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        profile     (profile object)        Profile object

    Returns
    -------
        Dew point (C) at the given pressure
    '''
    return interp_from_pres(p, profile, 3)


def vtmp(p, profile):
    '''
    Interpolates the given data to calculate a virtual temperature at a
    given pressure

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        profile     (profile object)        Profile object

    Returns
    -------
        Virtual temperature (C) at the given pressure
    '''
    t = temp(p, profile)
    td = dwpt(p, profile)
    return virtemp(p, t, td)


def components(p, profile):
    '''
    Interpolates the given data to calculate the U and V components at a
    given pressure

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        profile     (profile object)        Profile object

    Returns
    -------
        U/V components at the given pressure
    '''
    U = interp_from_pres(p, profile, 4)
    V = interp_from_pres(p, profile, 5)
    return U, V


def vec(p, U, V):
    '''
    Interpolates the given component data to a given pressure level
    and returns the interpolated direction and speed

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        U           (list of floats)        U-components
        V           (list of floats)        V-components

    Returns
    -------
        dir         (float)                 Direction
        mag         (float)                 Magnitude
    '''
    u = interp_from_pres(p, U)
    v = interp_from_pres(p, V)
    return vector.comp2vec(u, v)



def interp_from_hght(h, profile, ind):
    '''
    General interpolation routine for height coordinates

    Inputs
    ------
        h           (float)                 Height (m) of a level
        profile     (profile object)        Profile object
        ind         (integer)               Index of variable to interpolate

    Returns
    -------
        Interpolated variable
    '''
    tptr = -1
    bptr = 0
    for i in range(profile.gNumLevels):
        if not QC(profile.gSndg[i][ind]): continue
        if profile.gSndg[i][1] < h:
            bptr = i
        elif profile.gSndg[i][1] == h:
            return profile.gSndg[i][ind]
        else:
            tptr = i
            nm1 = h - profile.gSndg[bptr][1]
            nm2 = profile.gSndg[tptr][1] - profile.gSndg[bptr][1]
            val = profile.gSndg[bptr][ind] + ((nm1 / nm2) * \
                (profile.gSndg[tptr][ind] - profile.gSndg[bptr][ind]))
            return val


def interp_from_pres(p, profile, ind):
    '''
    General interpolation routine for pressure coordinates

    Inputs
    ------
        p           (float)                 Pressure (hPa) of a level
        profile     (profile object)        Profile object
        ind         (integer)               Index of variable to interpolate

    Returns
    -------
        Interpolated variable
    '''
    tptr = -1
    bptr = 0
    for i in range(profile.gNumLevels):
        if not QC(profile.gSndg[i][ind]): continue
        if profile.gSndg[i][0] > p:
            bptr = i
        elif profile.gSndg[i][0] == p:
            return profile.gSndg[i][ind]
        else:
            tptr = i
            nm1 = profile.gSndg[tptr][ind] - profile.gSndg[bptr][ind]
            nm2 = math.log(profile.gSndg[bptr][0] / profile.gSndg[tptr][0])
            nm3 = math.log(profile.gSndg[bptr][0] / p)
            return profile.gSndg[bptr][ind] + ((nm3 / nm2) * nm1)


def agl(h, profile):
    '''
    Convert a height from mean sea-level (MSL) to above ground-level (AGL)

    Inputs
    ------
        h           (float)                 Height of a level
        profile     (profile object)        Profile object

    Returns
    -------
        Converted height
    '''
    if not QC(h): return RMISSD
    return h - profile.gSndg[profile.sfc][1]


def msl(h, profile):
    '''
    Convert a height from above ground-level (AGL) to mean sea-level (MSL)

    Inputs
    ------
        h           (float)                 Height of a level
        profile     (profile object)        Profile object

    Returns
    -------
        Converted height
    '''
    if not QC(h): return RMISSD
    return h + profile.gSndg[profile.sfc][1]


