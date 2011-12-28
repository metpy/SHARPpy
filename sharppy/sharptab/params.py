''' Thermodynamic Parameter Routines '''
import math
from sharppy.sharptab import interp, vector, thermo, winds
from sharppy.sharptab import parcel as ppcl
from sharppy.sharptab.constants import *

__all__ = ['k_index']


def k_index(profile):
    '''
    Calculates the K-Index from a profile object

    Inputs
    ------
        profile (profile object)    Profile Object

    Returns
    -------
        kind    (float)             K-Index
    '''
    t8 = interp.temp(850., profile)
    t7 = interp.temp(700., profile)
    t5 = interp.temp(500., profile)
    td8 = interp.dwpt(850., profile)
    td7 = interp.dwpt(700., profile)
    if not QC(t8) or not QC(t7) or not QC(t5) or not QC(td8) or not QC(td7):
        return RMISSD
    else:
        return t8 - t5 + td8 - (t7 - td7)


def t_totals(profile):
    '''
    Calculates the Total Totals Index from data in profile object. This is
    done by summing the Cross Totals and Vertical Totals.

    Inputs
    ------
        profile (profile object)    Profile Object

    Returns
    -------
        t_totals    (float)         Total Totals
    '''
    ct = c_totals(profile)
    vt = v_totals(profile)
    return ct + vt


def c_totals(profile):
    '''
    Calculates the Cross Totals Index from data in profile object.

    Inputs
    ------
        profile (profile object)    Profile Object

    Returns
    -------
        c_totals    (float)         Cross Totals
    '''
    t5 = interp.temp(500., profile)
    td8 = interp.dwpt(850., profile)
    return td8 - t5


def v_totals(profile):
    '''
    Calculates the Vertical Totals Index from data in profile object.

    Inputs
    ------
        profile (profile object)    Profile Object

    Returns
    -------
        v_totals    (float)         Vertical Totals
    '''
    t8 = interp.temp(850., profile)
    t5 = interp.temp(500., profile)
    return t8 - t5


def precip_water(lower, upper, profile):
    '''
    Calculates the precipitable water from a profile object within the
    specified layer. The default layer (lower=-1 & upper=-1) is defined to
    be surface to 400 hPa.

    Inputs
    ------
        lower       (float)             Lower pressure level
        upper       (float)             Upper pressure level
        profile     (profile object)    Profile Object

    Returns
    -------
        pwat        (float)             Precipitable Water (in)
    '''
    if lower == -1: lower = profile.gSndg[profile.sfc][0]
    if upper == -1: upper = 400.

    # Find lower and upper ind bounds for looping
    i = 0
    while profile.gSndg[i][0] > lower: i+=1
    lptr = i
    while profile.gSndg[i][0] > upper: i+=1
    uptr = i

    # Start with interpolated bottom level
    d1 = interp.dwpt(lower, profile)
    p1 = lower
    w1 = thermo.mixratio(p1, d1)

    # Loop through every level that has a dew point
    pwat = 0
    for i in range(lptr, uptr+1):
        if QC(profile.gSndg[i][3]):
            p2 = profile.gSndg[i][0]
            w2 = thermo.mixratio(p2, profile.gSndg[i][3])
            pwat += ((w1 + w2) / 2.) * (p1 - p2)
            p1 = p2
            w1 = w2

    # Finish with interpolated top level
    d2 = interp.dwpt(upper, profile)
    p2 = upper
    w2 = thermo.mixratio(p2, d2)
    pwat += ((w1 + w2) / 2.) * (p1 - p2)

    return pwat * 0.00040173


def parcel(lower, upper, lplvals, profile):
    '''
    Lifts the specified parcel, calculated various levels and parameters from
    the profile object. B+/B- are calculated based on the specified layer.

    !! All calculations use the virtual temperature correction unless noted. !!

    Inputs
    ------
        lower       (float)                 Lower-bound lifting level (hPa)
        upper       (float)                 Upper-bound lifting level
        lplvals     (defined parcel object) Defined Parcel Object
        profile     (profile object)        Profile Object
    Returns
    -------
        pcl         (parcel object)         Parcel Object
    '''
    pcl = ppcl.Parcel(-1, -1, lplvals)
    pres = pcl.pres
    temp = pcl.temp
    dwpt = pcl.dwpt
    if profile.gNumLevels < 1: return pcl

    lyre = -1
    cap_strength = RMISSD
    cap_strengthpres = RMISSD
    li_max = RMISSD
    li_maxpres = RMISSD
    totp = 0.
    totn = 0.
    tote = 0.

    # See if default layer is specified
    if lower == -1:
        lower = profile.gSndg[profile.sfc][0]
        pcl.blayer = lower
    if upper == -1:
        upper = profile.gSndg[profile.gNumLevels-1][0]
        pcl.tlayer = upper

    # Make sure that this is a valid layer
    if lower > pres:
        lower = pres
        pcl.blayer = lower
    if not QC(interp.vtmp(lower, profile)) or \
       not QC(interp.vtmp(upper, profile)):
        return RMISSD

    # Begin with the Mixing Layer
    te1 = interp.vtmp(pres, profile)
    pe1 = lower
    h1 = interp.hght(pe1, profile)
    tp1 = thermo.virtemp(pres, temp, dwpt)
    te1 = tp1

    # Lift parcel and return LCL pres (hPa) and LCL temp (c)
    pe2, tp2 = thermo.drylift(pres, temp, dwpt)
    blupper = pe2       # Define top of layer as LCL pres
    h2 = interp.hght(pe2, profile)
    te2 = interp.vtmp(pe2, profile)
    pcl.lclpres = pe2
    pcl.lclhght = interp.agl(h2, profile)

    # Calculate lifted parcel theta for use in iterative CINH loop below
    # RECALL: lifted parcel theta is CONSTANT from LPL to LCL
    theta_parcel = thermo.theta(pe2, tp2, 1000.)

    # Environmental theta and mixing ratio at LPL
    bltheta = thermo.theta(pres, interp.temp(pres, profile), 1000.)
    blmr = thermo.mixratio(pres, dwpt)

    # ACCUMULATED CINH IN MIXING LAYER BELOW THE LCL
    # This will be done in 10mb increments, and will use the virtual
    # temperature correction where possible
    pinc = -10
    for pp in range(int(lower), int(blupper), int(pinc)):
        pp1 = pp
        pp2 = pp + pinc
        if pp2 < blupper: pp2 = blupper
        dz = interp.hght(pp2, profile) - interp.hght(pp1, profile)

        # Calculate difference between Tv_parcel and Tv_environment at top
        # and bottom of 10mb layers. Make use of constant lifted parcel
        # theta and mixing ratio from LPL to LCL
        tv_env_bot = thermo.virtemp(pp1, thermo.theta(pp1,
            interp.temp(pp1, profile), 1000.), interp.dwpt(pp1, profile))
        tdef1 = (thermo.virtemp(pp1, theta_parcel,
            thermo.temp_at_mixrat(blmr, pp1))) / (thermo.ctok(tv_env_bot))
        tv_env_top = thermo.virtemp(pp2, thermo.theta(pp2,
            interp.temp(pp2, profile), 1000.), interp.dwpt(pp2, profile))
        tdef2 = (thermo.virtemp(pp2, theta_parcel,
            thermo.temp_at_mixrat(blmr, pp2))) / (thermo.ctok(tv_env_bot))
        lyre = G * (tdef1 + tdef2) / 2. * dz
        if lyre < 0:
            totn += lyre

    # Move the bottom layer to the top of the boundary layer
    if lower > pe2:
        pcl.blayer = pe2

    # Calculate height of various temperature levels
    p10c = temp_lvl(-10., profile)
    p30c = temp_lvl(-30., profile)
    hgt10c = interp.hght(p10c, profile)
    hgt30c = interp.hght(p30c, profile)

    # Find lowest observation in layer
    i = 0
    while profile.gSndg[i][0] > lower:
        if i == profile.gNumLevels-1: break
        i += 1
    while not QC(profile.gSndg[i][3]):
        if i == profile.gNumLevels-1: break
        i += 1
    lptr = i
    if profile.gSndg[i][0] == lower:
        if i != profile.gNumLevels-1: lptr += 1

    # Find highest observation in layer
    i = profile.gNumLevels-1
    while profile.gSndg[i][0] < upper:
        if i < lptr: break
        i -= 1
    uptr = i
    if profile.gSndg[i][0] == upper:
        if i > lptr: uptr -= 1

    # START WITH INTERPOLATED BOTTOM LAYER
    # Begin moist ascent from lifted parcel LCL (pe2, tp2)
    tp1 = thermo.wetlift(pe2, tp2, pe1)
    lyre = 0
    lyrlast = 0
    for i in range(lptr, profile.gNumLevels):
        if not QC(profile.gSndg[i][2]): continue
        pe2 = profile.gSndg[i][0]
        h2 = profile.gSndg[i][1]
        te2 = interp.vtmp(pe2, profile)
        tp2 = thermo.wetlift(pe1, tp1, pe2)
        tdef1 = (thermo.virtemp(pe1, tp1, tp1) - te1) / thermo.ctok(te1)
        tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / thermo.ctok(te2)
        lyrlast = lyre
        lyre = G * (tdef1 + tdef2) / 2. * (h2 - h1)

        # Add layer energy to total positive if lyre > 0
        if lyre > 0:
            totp += lyre
        # Add layer energy to total negative if lyre < 0, only up to EL
        else:
            if pe2 > 500.: totn += lyre

        # Check for Max LI
        mli = thermo.virtemp(pe2, tp2, tp2) - te2
        if  mli > li_max:
            li_max = mli
            li_maxpres = pe2

        # Check for Max Cap Strength
        mcap = te2 - mli
        if mcap > cap_strength:
            cap_strength = mcap
            cap_strengthpres = pe2

        tote += lyre
        pelast = pe1
        pe1 = pe2
        h1 = h2
        te1 = te2
        tp1 = tp2

        # Is this the top of the specified layer
        if i >= uptr and not QC(pcl.bplus):
            pe3 = pe1
            h3 = h1
            te3 = te1
            tp3 = tp1
            lyrf = lyre
            if lyrf > 0:
                pcl.bplus = totp - lyrf
                pcl.bminus = totn
            else:
                pcl.bplus = totp
                if pe2 > 500.: pcl.bminus = totn + lyrf
                else: pcl.bminus = totn

            pe2 = upper
            h2 = interp.hght(pe2, profile)
            te2 = interp.vtmp(pe2, profile)
            tp2 = thermo.wetlift(pe3, tp3, pe2)
            tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / thermo.ctok(te3)
            tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / thermo.ctok(te2)
            lyrf = G * (tdef3 + tdef2) / 2. * (h2 - h3)
            if lyrf > 0:
                pcl.bplus += lyrf
            else:
                if pe2 > 500.: pcl.bminus += lyrf

            if pcl.bplus == 0: pcl.bminus = 0.

        # Is this the freezing level
        if te2 < 0. and not QC(pcl.bfzl):
            pe3 = pelast
            h3 = interp.hght(pe3, profile)
            te3 = interp.vtmp(pe3, profile)
            tp3 = thermo.wetlift(pe1, tp1, pe3)
            lyrf = lyre
            if lyrf > 0.: pcl.bfzl = totp - lyrf
            else: pcl.bfzl = totp
            pe2 = temp_lvl(0., profile)
            if QC(pe2):
                h2 = interp.hght(pe2, profile)
                te2 = interp.vtmp(pe2, profile)
                tp2 = thermo.wetlift(pe3, tp3, pe2)
                tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / \
                    thermo.ctok(te3)
                tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / \
                    thermo.ctok(te2)
                lyrf = G * (tdef3 + tdef2) / 2. * (h2 - h3)
                if lyrf > 0: pcl.bfzl += lyrf

        # Is this the -10C level
        if te2 < -10. and not QC(pcl.wm10c):
            pe3 = pelast
            h3 = interp.hght(pe3, profile)
            te3 = interp.vtmp(pe3, profile)
            tp3 = thermo.wetlift(pe1, tp1, pe3)
            lyrf = lyre
            if lyrf > 0.: pcl.wm10c = totp - lyrf
            else: pcl.wm10c = totp
            pe2 = temp_lvl(-10., profile)
            if QC(pe2):
                h2 = interp.hght(pe2, profile)
                te2 = interp.vtmp(pe2, profile)
                tp2 = thermo.wetlift(pe3, tp3, pe2)
                tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / \
                    thermo.ctok(te3)
                tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / \
                    thermo.ctok(te2)
                lyrf = G * (tdef3 + tdef2) / 2. * (h2 - h3)
                if lyrf > 0: pcl.wm10c += lyrf

        # Is this the -30C level
        if te2 < -30. and not QC(pcl.wm30c):
            pe3 = pelast
            h3 = interp.hght(pe3, profile)
            te3 = interp.vtmp(pe3, profile)
            tp3 = thermo.wetlift(pe1, tp1, pe3)
            lyrf = lyre
            if lyrf > 0.: pcl.wm30c = totp - lyrf
            else: pcl.wm30c = totp
            pe2 = temp_lvl(-30., profile)
            if QC(pe2):
                h2 = interp.hght(pe2, profile)
                te2 = interp.vtmp(pe2, profile)
                tp2 = thermo.wetlift(pe3, tp3, pe2)
                tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / \
                    thermo.ctok(te3)
                tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / \
                    thermo.ctok(te2)
                lyrf = G * (tdef3 + tdef2) / 2. * (h2 - h3)
                if lyrf > 0: pcl.wm30c += lyrf

        # Is this the 3km level
        if pcl.lclhght < 3000.:
            h = interp.agl(interp.hght(pe2, profile), profile)
            if h >= 3000. and not QC(pcl.b3km):
                pe3 = pelast
                h3 = interp.hght(pe3, profile)
                te3 = interp.vtmp(pe3, profile)
                tp3 = thermo.wetlift(pe1, tp1, pe3)
                lyrf = lyre
                if lyrf > 0: pcl.b3km = totp - lyrf
                else: pcl.b3km = totp
                h2 = interp.msl(3000., profile)
                pe2 = interp.pres(h2, profile)
                if QC(pe2):
                    te2 = interp.vtmp(pe2, profile)
                    tp2 = thermo.wetlift(pe3, tp3, pe2)
                    tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / \
                        thermo.ctok(te3)
                    tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / \
                        thermo.ctok(te2)
                    lyrf = G * (tdef3 + tdef2) / 2. * (h2 - h3)
                    if lyrf > 0: pcl.b3km += lyrf
        else: pcl.b3km = 0.

        # Is this the 6km level
        if pcl.lclhght < 6000.:
            h = interp.agl(interp.hght(pe2, profile), profile)
            if h >= 6000. and not QC(pcl.b6km):
                pe3 = pelast
                h3 = interp.hght(pe3, profile)
                te3 = interp.vtmp(pe3, profile)
                tp3 = thermo.wetlift(pe1, tp1, pe3)
                lyrf = lyre
                if lyrf > 0: pcl.b6km = totp - lyrf
                else: pcl.b6km = totp
                h2 = interp.msl(6000., profile)
                pe2 = interp.pres(h2, profile)
                if QC(pe2):
                    te2 = interp.vtmp(pe2, profile)
                    tp2 = thermo.wetlift(pe3, tp3, pe2)
                    tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / \
                        thermo.ctok(te3)
                    tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / \
                        thermo.ctok(te2)
                    lyrf = G * (tdef3 + tdef2) / 2. * (h2 - h3)
                    if lyrf > 0: pcl.b6km += lyrf
        else: pcl.b6km = 0.

        # LFC Possibility
        if lyre >= 0. and lyrlast <= 0.:
            tp3 = tp1
            te3 = te1
            pe2 = pe1
            pe3 = pelast
            while interp.vtmp(pe3, profile) > thermo.virtemp(pe3,
                thermo.wetlift(pe2, tp3, pe3), thermo.wetlift(pe2, tp3, pe3)):
                    pe3 -= 5
            pcl.lfcpres = pe3
            pcl.lfchght = interp.agl(interp.hght(pe3, profile), profile)
            cinh_old = totn
            tote = 0.
            pcl.elpres = RMISSD

            if cap_strength < 0.: cap_strength = 0.
            pcl.cap = cap_strength
            pcl.cappres = cap_strengthpres

        # EL Possibility
        if lyre <=0. and lyrlast >= 0.:
            tp3 = tp1
            te3 = te1
            pe2 = pe1
            pe3 = pelast
            while interp.vtmp(pe3, profile) < thermo.virtemp(pe3,
                thermo.wetlift(pe2, tp3, pe3), thermo.wetlift(pe2, tp3, pe3)):
                    pe3 -= 5
            pcl.elpres = pe3
            pcl.elhght = interp.agl(interp.hght(pe3, profile), profile)
            pcl.mplpres = RMISSD
            pcl.limax = -li_max
            pcl.limaxpress = li_maxpres

        # MPL Possibility
        if tote < 0. and not QC(pcl.mplpres) and QC(pcl.elpres):
            pe3 = pelast
            h3 = interp.hght(pe3, profile)
            te3 = interp.vtmp(pe3, profile)
            tp3 = thermo.wetlift(pe1, tp1, pe3)
            totx = tote - lyre
            pe2 = pelast
            while totx > 0:
                pe2 -= 1
                te2 = interp.vtmp(pe2, profile)
                tp2 = thermo.wetlift(pe3, tp3, pe2)
                h2 = interp.hght(pe2, profile)
                tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / \
                    thermo.ctok(te3)
                tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / \
                    thermo.ctok(te2)
                lyrf = G * (tdef3 + tdef2) / 2. * (h2 - h3)
                totx += lyrf
                tp3 = tp2
                te3 = te2
                pe3 = pe2
            pcl.mplpres = pe2
            pcl.mplhght = interp.agl(interp.hght(pe2, profile), profile)

        # 500 hPa Lifted Index
        if profile.gSndg[i][0] <= 500. and pcl.li5 == RMISSD:
            a = interp.vtmp(500., profile)
            b = thermo.wetlift(pe1, tp1, 500.)
            pcl.li5 = a - thermo.virtemp(500, b, b)

        # 300 hPa Lifted Index
        if profile.gSndg[i][0] <= 300. and pcl.li3 == RMISSD:
            a = interp.vtmp(300., profile)
            b = thermo.wetlift(pe1, tp1, 300.)
            pcl.li3 = a - thermo.virtemp(300, b, b)

    # Calculate BRN if available
    pcl = bulk_rich(pcl, profile)

    pcl.bminus = cinh_old
    if pcl.bplus == 0: pcl.bminus = 0.
    return pcl


def temp_lvl(temp, profile):
    '''
    Calculates the level (hPa) of the first occurrence of the specified
    temperature.

    Inputs
    ------
        temp        (float)             Temperature being searched (C)
        profile     (profile object)    Profile Object

    Returns
    -------
        Level of the temperature (hPa)
    '''
    for i in range(profile.gNumLevels):
        if QC(profile.gSndg[i][2]) and profile.gSndg[i][2] <= temp:
            if i == 0: return RMISSD
            if profile.gSndg[i][2] == temp: return profile.gSndg[i][0]
            p0 = profile.gSndg[i-1][0]
            t0 = profile.gSndg[i-1][2]
            nm1 = temp - t0
            nm2 = profile.gSndg[i][2] - t0
            nm3 = math.log(profile.gSndg[i][0] / p0)
            return p0 * math.exp((nm1 / nm2) * nm3)
    return RMISSD


def bulk_rich(pcl, profile):
    '''
    Calculates the Bulk Richardson Number for a given parcel.

    Inputs
    ------
        pcl         (parcel object)         Parcel Object
        profile     (profile object)        Profile Object

    Returns
    -------
        Bulk Richardson Number
    '''
    # Make sure parcel is initialized
    if pcl.lplvals.flag == RMISSD: return RMISSD
    elif pcl.lplvals.flag > 0 and pcl.lplvals.flag < 4:
        ptop = interp.pres(interp.msl(6000., profile), profile)
        pbot = profile.gSndg[profile.sfc][0]
    else:
        h0 = interp.hght(pcl.lplvals.pres, profile)
        pbot = interp.pres(h0 - 500., profile)
        if not QC(pbot): pbot = profile.gSndg[profile.sfc][0]
        h1 = interp.pres(h1, profile)
        ptop = interp.pres(h1 + 6000., profile)

    # Calculate lowest 500m mean wind
    p = interp.pres(interp.hght(pbot, profile) + 500., profile)
    mnlu, mnlv = winds.mean_wind(pbot, p, profile)

    # Calculate the 6000m mean wind
    mnuu, mnuv = winds.mean_wind(pbot, ptop, profile)

    # Make sure CAPE and Shear are available
    if not QC(pcl.bplus) or not QC(mnlu) or not QC(mnuu): return RMISSD

    # Calculate shear between levels
    dx = mnuu - mnlu
    dy = mnuv - mnlv

    pcl.brnshear = KTS2MS(vector.comp2vec(dx, dy)[1])
    pcl.brnshear = pcl.brnshear**2 / 2.
    pcl.brn = pcl.bplus / pcl.brnshear
    return pcl




