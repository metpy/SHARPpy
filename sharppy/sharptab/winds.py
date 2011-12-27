''' Wind Manipulation Routines '''
import math
from sharppy.sharptab import interp, vector
from sharppy.sharptab.constants import *

__all__ = ['mean_wind', 'mean_wind_npw', 'sr_wind', 'sr_wind_npw',
           'wind_shear', 'helicity', 'max_wind', 'corfidi_mcs_motion',
           'bunkers_right_mover', 'bunkers_left_mover', 'mbe_vectors']


def mean_wind(pbot, ptop, profile, psteps=20, stu=0, stv=0):
    '''
    Calculates a pressure-weighted mean wind through a layer. The default
    layer is 850 to 200 hPa.

    Inputs
    ------
        pbot    (float)             Pressure of the bottom level (hPa)
        ptop    (float)             Pressure of the top level (hPa)
        profile (profile object)    Profile Object
        psteps  (int; optional)     Number of steps to loop through (int)
        stu     (float; optional)   U-component of storm-motion vector
        stv     (float; optional)   V-component of storm-motion vector

    Returns
    -------
        mnu      (float)            U-component
        mnv      (float)            V-component
    '''
    if pbot == -1: lower = 850.
    if ptop == -1: upper = 200.
    pinc = int((pbot - ptop) / psteps)
    if psteps < 1:
        u1 = (interp.interp_from_pres(pbot, profile, 4) - stu) * pbot
        u2 = (interp.interp_from_pres(ptop, profile, 4)-stu) * ptop
        v1 = (interp.interp_from_pres(pbot, profile, 5)-stv) * pbot
        v2 = (interp.interp_from_pres(ptop, profile, 5)-stv) * ptop
        usum = u1 + u2
        vsum = v1 + v2
        wgt = pbot + ptop
    else:
        wgt = 0
        usum = 0
        vsum = 0
        for p in range(int(pbot), int(ptop), -pinc):
            usum += (interp.interp_from_pres(p, profile, 4) - stu) * p
            vsum += (interp.interp_from_pres(p, profile, 5) - stv) * p
            wgt += p

    return float(usum / wgt), float(vsum / wgt)


def mean_wind_npw(pbot, ptop, profile, psteps=20, stu=0, stv=0):
    '''
    Calculates a pressure-weighted mean wind through a layer. The default
    layer is 850 to 200 hPa.

    Inputs
    ------
        pbot    (float)             Pressure of the bottom level (hPa)
        ptop    (float)             Pressure of the top level (hPa)
        profile (profile object)    Profile Object
        psteps  (int; optional)     Number of steps to loop through (int)
        stu     (float; optional)   U-component of storm-motion vector
        stv     (float; optional)   V-component of storm-motion vector

    Returns
    -------
        mnu      (float)            U-component
        mnv      (float)            V-component
    '''
    if pbot == -1: lower = 850.
    if ptop == -1: upper = 200.
    pinc = int((pbot - ptop) / psteps)
    if pinc < 1:
        u1 = (interp.interp_from_pres(pbot, profile, 4) - stu) * pbot
        u2 = (interp.interp_from_pres(pbot, profile, 4) - stu) * ptop
        v1 = (interp.interp_from_pres(pbot, profile, 5) - stv) * pbot
        v2 = (interp.interp_from_pres(pbot, profile, 5) - stv) * ptop
        usum = u1 + u2
        vsum = v1 + v2
        wgt = 2
    else:
        wgt = 0
        usum = 0
        vsum = 0
        for p in range(int(pbot), int(ptop), -pinc):
            usum += (interp.interp_from_pres(p, profile, 4) - stu)
            vsum += (interp.interp_from_pres(p, profile, 5) - stv)
            wgt += 1

    return float(usum / wgt), float(vsum / wgt)


def sr_wind(pbot, ptop, stu, stv, profile, psteps=20):
    '''
    Calculates a pressure-weighted mean storm-relative wind through a layer.
    The default layer is 850 to 200 hPa. This is a thin wrapper around
    mean_wind().

    Inputs
    ------
        pbot    (float)             Pressure of the bottom level (hPa)
        ptop    (float)             Pressure of the top level (hPa)
        stu     (float)             U-component of storm-motion vector
        stv     (float)             V-component of storm-motion vector
        profile (profile object)    Profile Object
        psteps  (int; optional)     Number of steps to loop through (int)

    Returns
    -------
        mnu      (float)            U-component
        mnv      (float)            V-component
    '''
    return mean_wind(pbot, ptop, profile, psteps, stu, stv)


def sr_wind_npw(pbot, ptop, stu, stv, profile, psteps=20):
    '''
    Calculates a non-pressure-weighted storm-relative mean wind through a
    layer. The default layer is 850 to 200 hPa. This is a thin wrapper
    around mean_wind_npw().

    Inputs
    ------
        pbot    (float)             Pressure of the bottom level (hPa)
        ptop    (float)             Pressure of the top level (hPa)
        stu     (float)             U-component of storm-motion vector
        stv     (float)             V-component of storm-motion vector
        profile (profile object)    Profile Object
        psteps  (int; optional)     Number of steps to loop through (int)

    Returns
    -------
        mnu      (float)            U-component
        mnv      (float)            V-component
    '''
    return mean_wind_npw(pbot, ptop, profile, psteps, stu, stv)


def wind_shear(pbot, ptop, profile):
    '''
    Calculates the shear between the wind at (pbot) and (ptop).

    Inputs
    ------
        pbot    (float)             Pressure of the bottom level (hPa)
        ptop    (float)             Pressure of the top level (hPa)
        profile (profile object)    Profile Object

    Returns
    -------
        shu      (float)            U-component
        shv      (float)            V-component
    '''
    ubot = interp.interp_from_pres(pbot, profile, 4)
    utop = interp.interp_from_pres(ptop, profile, 4)
    shu = utop - ubot
    vbot = interp.interp_from_pres(pbot, profile, 5)
    vtop = interp.interp_from_pres(ptop, profile, 5)
    shv = vtop - vbot
    return shu, shv


def helicity(lower, upper, profile, stu=0, stv=0):
    '''
    Calculates the relative helicity (m2/s2) of a layer from lower to upper.
    If storm-motion vector is supplied, storm-relative helicity, both
    positve and negative, is returned.

    Inputs
    ------
        lower       (float)             Bottom level of layer (m, AGL)
        upper       (float)             Top level of layer (m, AGL)
        profile     (profile object)    Profile Object
        stu         (float; optional)   U-component of storm-motion
        stv         (float; optional)   V-component of storm-motion

    Returns
    -------
        phel+nhel   (float)             Combined Helicity (m2/s2)
        phel        (float)             Positive Helicity (m2/s2)
        nhel        (float)             Negative Helicity (m2/s2)
    '''
    lower = interp.msl(lower, profile)
    upper = interp.msl(upper, profile)
    plower = interp.interp_from_hght(lower, profile, 0)
    pupper = interp.interp_from_hght(upper, profile, 0)
    phel = 0
    nhel = 0

    # Find lower and upper ind bounds for looping
    i = 0
    while interp.msl(profile.gSndg[i][1], profile) < lower: i+=1
    lptr = i
    while interp.msl(profile.gSndg[i][1], profile) < upper: i+=1
    uptr = i

    # Integrate from interpolated bottom level to iptr level
    sru1 = KTS2MS(interp.interp_from_pres(plower, profile, 4) - stu)
    srv1 = KTS2MS(interp.interp_from_pres(plower, profile, 5) - stv)

    # Loop through levels
    for i in range(lptr, uptr+1):
        lyrh = 0
        if QC(profile.gSndg[i][4]) and QC(profile.gSndg[i][5]):
            sru2 = KTS2MS(profile.gSndg[i][4] - stu)
            srv2 = KTS2MS(profile.gSndg[i][5] - stv)

            lyrh = (sru2 * srv1) - (sru1 * srv2)
            if lyrh > 0: phel += lyrh
            else: nhel += lyrh
            sru1 = sru2
            srv1 = srv2

    # Integrate from tptr level to interpolated top level
    sru2 = KTS2MS(interp.interp_from_pres(pupper, profile, 4) - stu)
    srv2 = KTS2MS(interp.interp_from_pres(pupper, profile, 5) - stv)

    lyrh = (sru2 * srv1) - (sru1 * srv2)
    if lyrh > 0: phel += lyrh
    else: nhel += lyrh

    return phel+nhel, phel, nhel


def max_wind(lower, upper, profile):
    '''
    Finds the maximum wind speed of the layer given by lower and upper levels.
    In the event of the maximum wind speed occurring at multiple levels, the
    lowest level it occurs is returned.

    Inputs
    ------
        lower       (float)             Bottom level of layer (m, AGL)
        upper       (float)             Top level of layer (m, AGL)
        profile     (profile object)    Profile Object

    Returns
    -------
        p           (float)             Pressure level (hPa) of max wind speed
        maxu        (float)             Maximum U-component
        maxv        (float)             Maximum V-component
    '''
    if lower == -1: lower = profile.gSndg[profile.sfc][0]
    if upper == -1: upper = profile.gSndg[profile.gNumLevels-1][0]

    # Find lower and upper ind bounds for looping
    i = 0
    while profile.gSndg[i][0] > lower: i+=1
    lptr = i
    while profile.gSndg[i][0] > upper: i+=1
    uptr = i

    # Start with interpolated bottom level
    maxu = interp.interp_from_pres(lower, profile, 4)
    maxv = interp.interp_from_pres(lower, profile, 5)
    maxspd = vector.comp2vec(maxu, maxv)[1]
    p = lower

    # Loop through all levels in layer
    for i in range(lptr, uptr+1):
        if QC(profile.gSndg[i][0]) and QC(profile.gSndg[i][4]) and \
           QC(profile.gSndg[i][5]):
            spd = vector.comp2vec(profile.gSndg[i][4], profile.gSndg[i][5])[1]
            if spd > maxspd:
                maxspd = spd
                maxu = profile.gSndg[i][4]
                maxv = profile.gSndg[i][5]
                p = profile.gSndg[i][0]

    # Finish with interpolated top level
    tmpu = interp.interp_from_pres(upper, profile, 4)
    tmpv = interp.interp_from_pres(upper, profile, 5)
    tmpspd = vector.comp2vec(tmpu, tmpv)[1]
    if tmpspd > maxspd:
        maxu = tmpu
        maxv = tmpv
        maxspd = tmpspd
        p = upper

    return p, maxu, maxv


def corfidi_mcs_motion(profile):
    '''
    Calculated the Meso-beta Elements (Corfidi) Vectors

    Inputs
    ------
        profile     (profile object)    Profile Object

    Returns
    -------
        upu         (float)             U-component of the upshear vector
        upv         (float)             V-component of the upshear vector
        dnu         (float)             U-component of the downshear vector
        dnv         (float)             V-component of the downshear vector
    '''
    # Compute the tropospheric (850hPa-300hPa) mean wind
    mnu1, mnv1 = mean_wind_npw(850., 300., profile)

    # Compute the low-level (SFC-1500m) mean wind
    p_1p5km = interp.interp_from_hght(interp.msl(1500., profile), profile, 0)
    mnu2, mnv2 = mean_wind_npw(profile.gSndg[profile.sfc][0], p_1p5km, profile)

    # Compute the upshear vector
    upu = mnu1 - mnu2
    upv = mnv1 - mnv2

    # Compute the downshear vector
    dnu = mnu1 + upu
    dnv = mnv1 + upv

    return upu, upv, dnu, dnv


def bunkers_right_mover(profile):
    '''
    Compute the Bunkers Storm Motion for a Right Moving Supercell

    Inputs
    ------
        profile     (profile object)        Profile Object

    Returns
    -------
        stu         (float)                 Storm Motion U-component
        stv         (float)                 Storm Motion V-component
    '''
    d = 7.5     # Deviation value emperically derived as 7.5 m/s
    msl6km = interp.msl(6000., profile)
    p6km = interp.interp_from_hght(msl6km, profile, 0)

    # SFC-6km Mean Wind
    mnu6, mnv6 = mean_wind_npw(profile.gSndg[profile.sfc][0],
        p6km, profile, 25)

    # SFC-6km Shear Vector
    shru6, shrv6 = wind_shear(profile.gSndg[profile.sfc][0],
        p6km, profile)

    # Bunkers Right Motion
    tmp = d / vector.comp2vec(shru6, shrv6)[1]
    stu = mnu6 + (tmp * shrv6)
    stv = mnv6 - (tmp * shru6)

    return stu, stv


def bunkers_left_mover(profile):
    '''
    Compute the Bunkers Storm Motion for a Left Moving Supercell

    Inputs
    ------
        profile     (profile object)        Profile Object

    Returns
    -------
        stu         (float)                 Storm Motion U-component
        stv         (float)                 Storm Motion V-component
    '''
    d = MS2KTS(7.5)     # Deviation value emperically derived as 7.5 m/s
    msl6km = interp.msl(6000., profile)
    p6km = interp.interp_from_hght(msl6km, profile, 0)

    # SFC-6km Mean Wind
    mnu6, mnv6 = mean_wind_npw(profile.gSndg[profile.sfc][0],
        p6km, profile, 25)

    # SFC-6km Shear Vector
    shru6, shrv6 = wind_shear(profile.gSndg[profile.sfc][0],
        p6km, profile)

    # Bunkers Left Motion
    tmp = d / vector.comp2vec(shru6, shrv6)[1]
    stu = mnu6 - (tmp * shrv6)
    stv = mnv6 + (tmp * shru6)

    return stu, stv


def mbe_vectors(profile):
    '''
    Thin wrapper around corfidi_mcs_motion()

    Inputs
    ------
        profile     (profile object)    Profile Object

    Returns
    -------
        upu         (float)             U-component of the upshear vector
        upv         (float)             V-component of the upshear vector
        dnu         (float)             U-component of the downshear vector
        dnv         (float)             V-component of the downshear vector
    '''
    return corfidi_mcs_motion(profile)