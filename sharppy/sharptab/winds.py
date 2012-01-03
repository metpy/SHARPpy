''' Wind Manipulation Routines '''
import math
from sharppy.sharptab import interp, vector
from sharppy.sharptab.constants import *

__all__ = ['mean_wind', 'mean_wind_npw', 'sr_wind', 'sr_wind_npw',
           'wind_shear', 'helicity', 'max_wind', 'corfidi_mcs_motion',
           'non_parcel_bunkers_motion', 'mbe_vectors']


def mean_wind(pbot, ptop, prof, psteps=20, stu=0, stv=0):
    '''
    Calculates a pressure-weighted mean wind through a layer. The default
    layer is 850 to 200 hPa.

    Inputs
    ------
        pbot    (float)             Pressure of the bottom level (hPa)
        ptop    (float)             Pressure of the top level (hPa)
        prof    (profile object)    Profile Object
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
        u1, v1 = interp.components(pbot, prof)
        u2, v2 = interp.components(ptop, prof)
        u1 = (u1 - stu) * pbot
        v1 = (v1 - stv) * pbot
        u2 = (u2 - stu) * ptop
        v2 = (v2 - stv) * ptop
        usum = u1 + u2
        vsum = v1 + v2
        wgt = pbot + ptop
    else:
        wgt = 0
        usum = 0
        vsum = 0
        for p in range(int(pbot), int(ptop)+1, -pinc):
            utmp, vtmp = interp.components(p, prof)
            usum += (utmp - stu) * p
            vsum += (vtmp - stv) * p
            wgt += p

    return float(usum / wgt), float(vsum / wgt)


def mean_wind_npw(pbot, ptop, prof, psteps=20, stu=0, stv=0):
    '''
    Calculates a pressure-weighted mean wind through a layer. The default
    layer is 850 to 200 hPa.

    Inputs
    ------
        pbot    (float)             Pressure of the bottom level (hPa)
        ptop    (float)             Pressure of the top level (hPa)
        prof    (profile object)    Profile Object
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
        u1, v1 = interp.components(pbot, prof)
        u2, v2 = interp.components(ptop, prof)
        u1 = (u1 - stu) * pbot
        v1 = (v1 - stv) * pbot
        u2 = (u2 - stu) * ptop
        v2 = (v2 - stv) * ptop
        usum = u1 + u2
        vsum = v1 + v2
        wgt = 2
    else:
        wgt = 0
        usum = 0
        vsum = 0
        for p in range(int(pbot), int(ptop), -pinc):
            utmp, vtmp = interp.components(p, prof)
            usum += (utmp - stu)
            vsum += (vtmp - stv)
            wgt += 1

    return float(usum / wgt), float(vsum / wgt)


def sr_wind(pbot, ptop, stu, stv, prof, psteps=20):
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
        prof    (profile object)    Profile Object
        psteps  (int; optional)     Number of steps to loop through (int)

    Returns
    -------
        mnu      (float)            U-component
        mnv      (float)            V-component
    '''
    return mean_wind(pbot, ptop, prof, psteps, stu, stv)


def sr_wind_npw(pbot, ptop, stu, stv, prof, psteps=20):
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
        prof    (profile object)    Profile Object
        psteps  (int; optional)     Number of steps to loop through (int)

    Returns
    -------
        mnu      (float)            U-component
        mnv      (float)            V-component
    '''
    return mean_wind_npw(pbot, ptop, prof, psteps, stu, stv)


def wind_shear(pbot, ptop, prof):
    '''
    Calculates the shear between the wind at (pbot) and (ptop).

    Inputs
    ------
        pbot    (float)             Pressure of the bottom level (hPa)
        ptop    (float)             Pressure of the top level (hPa)
        prof    (profile object)    Profile Object

    Returns
    -------
        shu      (float)            U-component
        shv      (float)            V-component
    '''
    ubot, vbot = interp.components(pbot, prof)
    utop, vtop = interp.components(ptop, prof)
    shu = utop - ubot
    shv = vtop - vbot
    return shu, shv


def helicity(lower, upper, prof, stu=0, stv=0):
    '''
    Calculates the relative helicity (m2/s2) of a layer from lower to upper.
    If storm-motion vector is supplied, storm-relative helicity, both
    positve and negative, is returned.

    Inputs
    ------
        lower       (float)             Bottom level of layer (m, AGL)
        upper       (float)             Top level of layer (m, AGL)
        prof        (profile object)    Profile Object
        stu         (float; optional)   U-component of storm-motion
        stv         (float; optional)   V-component of storm-motion

    Returns
    -------
        phel+nhel   (float)             Combined Helicity (m2/s2)
        phel        (float)             Positive Helicity (m2/s2)
        nhel        (float)             Negative Helicity (m2/s2)
    '''
    lower = interp.msl(lower, prof)
    upper = interp.msl(upper, prof)
    plower = interp.pres(lower, prof)
    pupper = interp.pres(upper, prof)
    phel = 0
    nhel = 0

    # Find lower and upper ind bounds for looping
    i = 0
    while interp.msl(prof.gSndg[i][prof.zind], prof) < lower: i+=1
    lptr = i
    if interp.msl(prof.gSndg[i][prof.zind], prof) == lower: lptr+=1
    while interp.msl(prof.gSndg[i][prof.zind], prof) < upper: i+=1
    uptr = i


    # Integrate from interpolated bottom level to iptr level
    sru1, srv1 = interp.components(plower, prof)
    sru1 = KTS2MS(sru1 - stu)
    srv1 = KTS2MS(srv1 - stv)

    # Loop through levels
    for i in range(lptr, uptr+1):
        lyrh = 0
        if QC(prof.gSndg[i][prof.uind]) and QC(prof.gSndg[i][prof.vind]):
            sru2, srv2 = interp.components(prof.gSndg[i][prof.pind], prof)
            sru2 = KTS2MS(sru2 - stu)
            srv2 = KTS2MS(srv2 - stv)

            lyrh = (sru2 * srv1) - (sru1 * srv2)
            if lyrh > 0: phel += lyrh
            else: nhel += lyrh
            sru1 = sru2
            srv1 = srv2

    # Integrate from tptr level to interpolated top level
    sru2, srv2 = interp.components(pupper, prof)
    sru2 = KTS2MS(sru2 - stu)
    srv2 = KTS2MS(srv2 - stv)

    lyrh = (sru2 * srv1) - (sru1 * srv2)
    if lyrh > 0: phel += lyrh
    else: nhel += lyrh

    return phel+nhel, phel, nhel


def max_wind(lower, upper, prof):
    '''
    Finds the maximum wind speed of the layer given by lower and upper levels.
    In the event of the maximum wind speed occurring at multiple levels, the
    lowest level it occurs is returned.

    Inputs
    ------
        lower       (float)             Bottom level of layer (m, AGL)
        upper       (float)             Top level of layer (m, AGL)
        prof        (profile object)    Profile Object

    Returns
    -------
        p           (float)             Pressure level (hPa) of max wind speed
        maxu        (float)             Maximum U-component
        maxv        (float)             Maximum V-component
    '''
    if lower == -1: lower = prof.gSndg[prof.sfc][prof.pind]
    if upper == -1: upper = prof.gSndg[prof.gNumLevels-1][prof.pind]

    # Find lower and upper ind bounds for looping
    i = 0
    while prof.gSndg[i][prof.pind] > lower: i+=1
    lptr = i
    while prof.gSndg[i][prof.pind] > upper: i+=1
    uptr = i

    # Start with interpolated bottom level
    maxu, maxv = interp.components(lower, prof)
    maxspd = vector.comp2vec(maxu, maxv)[1]
    p = lower

    # Loop through all levels in layer
    for i in range(lptr, uptr+1):
        if QC(prof.gSndg[i][prof.pind]) and QC(prof.gSndg[i][prof.uind]) and \
           QC(prof.gSndg[i][prof.vind]):
            spd = vector.comp2vec(prof.gSndg[i][prof.uind],
                prof.gSndg[i][prof.vind])[1]
            if spd > maxspd:
                maxspd = spd
                maxu = prof.gSndg[i][prof.uind]
                maxv = prof.gSndg[i][prof.vind]
                p = prof.gSndg[i][prof.pind]

    # Finish with interpolated top level
    tmpu, tmpv = interp.components(upper, prof)
    tmpspd = vector.comp2vec(tmpu, tmpv)[1]
    if tmpspd > maxspd:
        maxu = tmpu
        maxv = tmpv
        maxspd = tmpspd
        p = upper

    return p, maxu, maxv


def corfidi_mcs_motion(prof):
    '''
    Calculated the Meso-beta Elements (Corfidi) Vectors

    Inputs
    ------
        prof        (profile object)    Profile Object

    Returns
    -------
        upu         (float)             U-component of the upshear vector
        upv         (float)             V-component of the upshear vector
        dnu         (float)             U-component of the downshear vector
        dnv         (float)             V-component of the downshear vector
    '''
    # Compute the tropospheric (850hPa-300hPa) mean wind
    mnu1, mnv1 = mean_wind_npw(850., 300., prof)

    # Compute the low-level (SFC-1500m) mean wind
    p_1p5km = interp.pres(interp.msl(1500., prof), prof)
    mnu2, mnv2 = mean_wind_npw(prof.gSndg[prof.sfc][prof.pind], p_1p5km, prof)

    # Compute the upshear vector
    upu = mnu1 - mnu2
    upv = mnv1 - mnv2

    # Compute the downshear vector
    dnu = mnu1 + upu
    dnv = mnv1 + upv

    return upu, upv, dnu, dnv


def non_parcel_bunkers_motion(prof):
    '''
    Compute the Bunkers Storm Motion for a Right Moving Supercell

    Inputs
    ------
        prof         (profile object)    Profile Object

    Returns
    -------
        rstu         (float)            Right Storm Motion U-component
        rstv         (float)            Right Storm Motion V-component
        lstu         (float)            Left Storm Motion U-component
        lstv         (float)            Left Storm Motion V-component
    '''
    d = MS2KTS(7.5)         # Deviation value emperically derived as 7.5 m/s
    msl6km = interp.msl(6000., prof)
    p6km = interp.pres(msl6km, prof)

    # SFC-6km Mean Wind
    mnu6, mnv6 = mean_wind_npw(prof.gSndg[prof.sfc][prof.pind],
        p6km, prof, 20)

    # SFC-6km Shear Vector
    shru6, shrv6 = wind_shear(prof.gSndg[prof.sfc][prof.pind],
        p6km, prof)

    # Bunkers Right Motion
    tmp = d / vector.comp2vec(shru6, shrv6)[1]
    rstu = mnu6 + (tmp * shrv6)
    rstv = mnv6 - (tmp * shru6)
    lstu = mnu6 - (tmp * shrv6)
    lstv = mnv6 + (tmp * shru6)

    return rstu, rstv, lstu, lstv


def mbe_vectors(prof):
    '''
    Thin wrapper around corfidi_mcs_motion()

    Inputs
    ------
        prof        (profile object)    Profile Object

    Returns
    -------
        upu         (float)             U-component of the upshear vector
        upv         (float)             V-component of the upshear vector
        dnu         (float)             U-component of the downshear vector
        dnv         (float)             V-component of the downshear vector
    '''
    return corfidi_mcs_motion(prof)