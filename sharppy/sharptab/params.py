''' Thermodynamic Parameter Routines '''
import math
from sharppy.sharptab import interp, vector, thermo, winds
from sharppy.sharptab.constants import *

__all__ = ['DefineParcel', 'Parcel', 'k_index', 't_totals', 'c_totals',
           'v_totals', 'precip_water', 'parcel', 'temp_lvl', 'bulk_rich',
           'max_temp', 'mean_mixratio', 'mean_theta', 'unstable_level',
           'effective_inflow_layer', 'bunkers_storm_motion', 'convective_temp',
           'esfc', 'lapse_rate']


class DefineParcel(object):
    '''
    Create a parcel from a supplied profile object.

    flag : int
        Parcel Selection
           -1: Use Previous Selection      # Not implemented as of 4/9/11
            1: Observed Surface Parcel
            2: Forecast Surface Parcel
            3: Most Unstable Parcel
            4: Mean Mixed Layer Parcel
            5: User Defined Parcel
            6: Mean Effective Layer Parcel

    pres : float
        Variable Pressure level (mb)
    '''
    def __init__(self, prof, flag, **kwargs):
        self.flag = flag
        if flag == 1:
            self.presval = kwargs.get('pres', prof.gSndg[prof.sfc])
            self.__sfc(prof, **kwargs)
        elif flag == 2:
            self.presval = kwargs.get('pres', 100)
            self.__fcst(prof, **kwargs)
        elif flag == 3:
            self.presval = kwargs.get('pres', 400)
            self.__mu(prof, **kwargs)
        elif flag == 4:
            self.presval = kwargs.get('pres', 100)
            self.__ml(prof, **kwargs)
        elif flag == 5:
            self.presval = kwargs.get('pres', 100)
            self.__user(prof, **kwargs)
        elif flag == 6:
            self.presval = kwargs.get('pres', 100)
            self.__effective(prof, **kwargs)
        else:
            print 'Defaulting to Surface Parcel'
            self.presval = kwargs.get('pres', prof.gSndg[prof.sfc])
            self.__sfc(prof, **kwargs)
        return

    def __sfc(self, prof, **kwargs):
        ''' Create a parcel using surface conditions '''
        self.desc = 'Surface Parcel'
        self.temp = prof.gSndg[prof.sfc][prof.tind]
        self.dwpt = prof.gSndg[prof.sfc][prof.tdind]
        self.pres = prof.gSndg[prof.sfc][prof.pind]
        return

    def __fcst(self, prof, **kwargs):
        ''' Create a parcel using forecast conditions '''
        self.desc = 'Forecast Surface Parcel'
        self.temp = max_temp(prof, -1)
        mmr = mean_mixratio(prof, -1, -1)
        self.dwpt = thermo.temp_at_mixrat(mmr, prof.gSndg[prof.sfc][prof.pind])
        self.pres = prof.gSndg[prof.sfc][prof.pind]
        return

    def __mu(self, prof, **kwargs):
        ''' Create the most unstable parcel within defined level '''
        self.desc = 'Most Unstable Parcel in Lowest %.2f hPa' % self.presval
        diff = prof.gSndg[prof.sfc][prof.pind] - self.presval
        self.pres = unstable_level(prof, -1, diff)
        self.temp = interp.temp(self.pres, prof)
        self.dwpt = interp.dwpt(self.pres, prof)
        return

    def __ml(self, prof, **kwargs):
        ''' Create the mixed-layer parcel; mixing over defined pressure '''
        self.desc = '%.2f hPa Mixed Layer Parcel' % self.presval
        self.pres = prof.gSndg[prof.sfc][prof.pind]
        diff = prof.gSndg[prof.sfc][prof.pind] - self.presval
        mtha = mean_theta(prof, -1, diff)
        mmr = mean_mixratio(prof, -1, diff)
        self.temp = thermo.theta(1000., mtha, self.pres)
        self.dwpt = thermo.temp_at_mixrat(mmr, self.pres)
        return

    def __user(self, prof, **kwargs):
        ''' Create a user-defined parcel '''
        self.desc = '%.2f hPa Parcel' % self.presval
        self.pres = self.presval
        if 'temp' in kwargs: self.temp = kwargs.get('temp')
        else: self.temp = interp.temp(self.pres, prof)
        if 'dwpt' in kwargs: self.dwpt = kwargs.get('dwpt')
        else: self.dwpt = interp.dwpt(self.pres, prof)
        return

    def __effective(self, prof):
        ''' Create the mean-effective layer parcel '''
        pbot, ptop = effective_inflow_layer(100, -250, prof)
        if pbot > 0:
            self.desc = '%.2f hPa Mean Effective Layer Centered at %.2f' % (
                pbot-ptop, (pbot+ptop)/2.)
            mtha = mean_theta(prof, pbot, ptop)
            mmr = mean_mixratio(prof, pbot, ptop)
            self.pres = (ptop + pbot) / 2.
            self.temp = thermo.theta(1000., mtha, self.pres)
            self.dwpt = thermo.temp_at_mixrat(mmr, self.pres)
        else:
            self.desc = 'Defaulting to Surface Layer'
            self.pres = prof.gSndg[prof.sfc][prof.pind]
            self.temp = prof.gSndg[prof.sfc][prof.tind]
            self.dwpt = prof.gSndg[prof.sfc][prof.tdind]
        if QC(pbot): self.pbot = pbot
        else: self.pbot = RMISSD
        if QC(ptop): self.ptop = ptop
        else: self.ptop = RMISSD
        return


class Parcel(object):
    ''' Initialize variables for a parcel '''
    def __init__(self, lower, upper, pres, temp, dwpt, missing=RMISSD,
        **kwargs):
        self.pres = pres
        self.temp = temp
        self.dwpt = dwpt
        self.blayer = lower
        self.tlayer = upper
        self.entrain = 0.
        self.lclpres = RMISSD
        self.lclhght = RMISSD
        self.lfcpres = RMISSD
        self.lfchght = RMISSD
        self.elpres = RMISSD
        self.elhght = RMISSD
        self.mplpres = RMISSD
        self.mplhght = RMISSD
        self.bplus = RMISSD
        self.bminus = RMISSD
        self.bfzl = RMISSD
        self.b3km = RMISSD
        self.b6km = RMISSD
        self.p0c = RMISSD
        self.pm10c = RMISSD
        self.pm20c = RMISSD
        self.pm30c = RMISSD
        self.hght0c = RMISSD
        self.hghtm10c = RMISSD
        self.hghtm20c = RMISSD
        self.hghtm30c = RMISSD
        self.wm10c = RMISSD
        self.wm20c = RMISSD
        self.wm30c = RMISSD
        self.li5 = RMISSD
        self.li3 = RMISSD
        self.brnshear = RMISSD
        self.brn = RMISSD
        self.limax = RMISSD
        self.limaxpres = RMISSD
        self.cap = RMISSD
        self.cappres = RMISSD
        for kw in kwargs: setattr(self, kw, kwargs.get(kw))
        return


def k_index(prof):
    '''
    Calculates the K-Index from a profile object

    Inputs
    ------
        prof    (profile object)    Profile Object

    Returns
    -------
        kind    (float)             K-Index
    '''
    t8 = interp.temp(850., prof)
    t7 = interp.temp(700., prof)
    t5 = interp.temp(500., prof)
    td7 = interp.dwpt(700., prof)
    td8 = interp.dwpt(850., prof)
    if not QC(t8) or not QC(t7) or not QC(t5) or not QC(td8) or not QC(td7):
        return RMISSD
    else:
        return t8 - t5 + td8 - (t7 - td7)


def t_totals(prof):
    '''
    Calculates the Total Totals Index from data in profile object. This is
    done by summing the Cross Totals and Vertical Totals.

    Inputs
    ------
        prof        (profile object)    Profile Object

    Returns
    -------
        t_totals    (float)             Total Totals
    '''
    ct = c_totals(prof)
    vt = v_totals(prof)
    return ct + vt


def c_totals(prof):
    '''
    Calculates the Cross Totals Index from data in profile object.

    Inputs
    ------
        prof        (profile object)    Profile Object

    Returns
    -------
        c_totals    (float)             Cross Totals
    '''
    t5 = interp.temp(500., prof)
    td8 = interp.dwpt(850., prof)
    return td8 - t5


def v_totals(prof):
    '''
    Calculates the Vertical Totals Index from data in profile object.

    Inputs
    ------
        prof        (profile object)    Profile Object

    Returns
    -------
        v_totals    (float)             Vertical Totals
    '''
    t8 = interp.temp(850., prof)
    t5 = interp.temp(500., prof)
    return t8 - t5


def precip_water(lower, upper, prof):
    '''
    Calculates the precipitable water from a profile object within the
    specified layer. The default layer (lower=-1 & upper=-1) is defined to
    be surface to 400 hPa.

    Inputs
    ------
        lower       (float)             Lower pressure level
        upper       (float)             Upper pressure level
        prof        (profile object)    Profile Object

    Returns
    -------
        pwat        (float)             Precipitable Water (in)
    '''
    if lower == -1: lower = prof.gSndg[prof.sfc][prof.pind]
    if upper == -1: upper = 400.

    # Find lower and upper ind bounds for looping
    i = 0
    while prof.gSndg[i][prof.pind] > lower: i+=1
    lptr = i
    while prof.gSndg[i][prof.pind] > upper: i+=1
    uptr = i

    # Start with interpolated bottom level
    d1 = interp.dwpt(lower, prof)
    p1 = lower
    w1 = thermo.mixratio(p1, d1)

    # Loop through every level that has a dew point
    pwat = 0
    for i in range(lptr, uptr+1):
        if QC(prof.gSndg[i][prof.tdind]):
            p2 = prof.gSndg[i][prof.pind]
            w2 = thermo.mixratio(p2, prof.gSndg[i][prof.tdind])
            pwat += ((w1 + w2) / 2.) * (p1 - p2)
            p1 = p2
            w1 = w2

    # Finish with interpolated top level
    d2 = interp.dwpt(upper, prof)
    p2 = upper
    w2 = thermo.mixratio(p2, d2)
    pwat += ((w1 + w2) / 2.) * (p1 - p2)

    return pwat * 0.00040173


def parcelx(lower, upper, pres, temp, dwpt, prof, **kwargs):
    '''
    Lifts the specified parcel, calculated various levels and parameters from
    the profile object. B+/B- are calculated based on the specified layer.

    !! All calculations use the virtual temperature correction unless noted. !!

    Inputs
    ------
        lower       (float)                 Lower-bound lifting level (hPa)
        upper       (float)                 Upper-bound lifting level
        pres        (float)                 Pressure of parcel to lift (hPa)
        temp        (float)                 Temperature of parcel to lift (C)
        dwpt        (float)                 Dew Point of parcel to lift (C)
        prof        (profile object)        Profile Object

    Returns
    -------
        pcl         (parcel object)         Parcel Object
    '''
    pcl = Parcel(-1, -1, pres, temp, dwpt)
    if 'lplvals' in kwargs: pcl.lplvals = kwargs.get('lplvals')
    else:
        lplvals = DefineParcel(prof, 5, pres=pres, temp=temp, dwpt=dwpt)
        pcl.lplvals = lplvals

    if prof.gNumLevels < 1: return pcl

    lyre = -1
    cap_strength = RMISSD
    cap_strengthpres = RMISSD
    li_max = RMISSD
    li_maxpres = RMISSD
    totp = 0.
    totn = 0.
    tote = 0.
    cinh_old = 0.

    # See if default layer is specified
    if lower == -1:
        lower = prof.gSndg[prof.sfc][prof.pind]
        pcl.blayer = lower
    if upper == -1:
        upper = prof.gSndg[prof.gNumLevels-1][prof.pind]
        pcl.tlayer = upper

    # Make sure that this is a valid layer
    if lower > pres:
        lower = pres
        pcl.blayer = lower
    if not QC(interp.vtmp(lower, prof)) or \
       not QC(interp.vtmp(upper, prof)):
        return RMISSD

    # Begin with the Mixing Layer
    te1 = interp.vtmp(pres, prof)
    pe1 = lower
    h1 = interp.hght(pe1, prof)
    tp1 = thermo.virtemp(pres, temp, dwpt)
    # te1 = tp1

    # Lift parcel and return LCL pres (hPa) and LCL temp (c)
    pe2, tp2 = thermo.drylift(pres, temp, dwpt)
    blupper = pe2       # Define top of layer as LCL pres
    h2 = interp.hght(pe2, prof)
    te2 = interp.vtmp(pe2, prof)
    pcl.lclpres = pe2
    pcl.lclhght = interp.agl(h2, prof)

    # Calculate lifted parcel theta for use in iterative CINH loop below
    # RECALL: lifted parcel theta is CONSTANT from LPL to LCL
    theta_parcel = thermo.theta(pe2, tp2, 1000.)

    # Environmental theta and mixing ratio at LPL
    bltheta = thermo.theta(pres, interp.temp(pres, prof), 1000.)
    blmr = thermo.mixratio(pres, dwpt)

    # ACCUMULATED CINH IN MIXING LAYER BELOW THE LCL
    # This will be done in 10mb increments, and will use the virtual
    # temperature correction where possible
    pinc = -10
    a = int(lower)
    b = int(blupper)
    for pp in range(a, b, int(pinc)):
        pp1 = pp
        pp2 = pp + pinc
        if pp2 < blupper: pp2 = blupper
        dz = interp.hght(pp2, prof) - interp.hght(pp1, prof)

        # Calculate difference between Tv_parcel and Tv_environment at top
        # and bottom of 10mb layers. Make use of constant lifted parcel
        # theta and mixing ratio from LPL to LCL
        tv_env_bot = thermo.virtemp(pp1, thermo.theta(pp1,
            interp.temp(pp1, prof), 1000.), interp.dwpt(pp1, prof))
        tdef1 = (thermo.virtemp(pp1, theta_parcel,
            thermo.temp_at_mixrat(blmr, pp1)) - tv_env_bot) / \
            (thermo.ctok(tv_env_bot))

        tv_env_top = thermo.virtemp(pp2, thermo.theta(pp2,
            interp.temp(pp2, prof), 1000.), interp.dwpt(pp2, prof))
        tdef2 = (thermo.virtemp(pp2, theta_parcel,
            thermo.temp_at_mixrat(blmr, pp2)) - tv_env_top) / \
            (thermo.ctok(tv_env_bot))

        lyre = G * (tdef1 + tdef2) / 2. * dz
        if lyre < 0: totn += lyre

    # Move the bottom layer to the top of the boundary layer
    if lower > pe2:
        lower = pe2
        pcl.blayer = lower

    # Calculate height of various temperature levels
    p0c = temp_lvl(0., prof)
    pm10c = temp_lvl(-10., prof)
    pm20c = temp_lvl(-20., prof)
    pm30c = temp_lvl(-30., prof)
    hgt0c = interp.hght(p0c, prof)
    hgtm10c = interp.hght(pm10c, prof)
    hgtm20c = interp.hght(pm20c, prof)
    hgtm30c = interp.hght(pm30c, prof)
    pcl.p0c = p0c
    pcl.pm10c = pm10c
    pcl.pm20c = pm20c
    pcl.pm30c = pm30c
    pcl.hght0c = hgt0c
    pcl.hghtm10c = hgtm10c
    pcl.hghtm20c = hgtm20c
    pcl.hghtm30c = hgtm30c

    # Find lowest observation in layer
    i = 0
    while prof.gSndg[i][prof.pind] > lower:
        if i == prof.gNumLevels-1: break
        i += 1
    while not QC(prof.gSndg[i][prof.tdind]):
        if i == prof.gNumLevels-1: break
        i += 1
    lptr = i
    if prof.gSndg[i][prof.pind] == lower:
        if i != prof.gNumLevels-1: lptr += 1

    # Find highest observation in layer
    i = prof.gNumLevels-1
    while prof.gSndg[i][prof.pind] < upper:
        if i < lptr: break
        i -= 1
    uptr = i
    if prof.gSndg[i][prof.pind] == upper:
        if i > lptr: uptr -= 1

    # START WITH INTERPOLATED BOTTOM LAYER
    # Begin moist ascent from lifted parcel LCL (pe2, tp2)
    pe1 = lower
    h1 = interp.hght(pe1, prof)
    te1 = interp.vtmp(pe1, prof)
    tp1 = thermo.wetlift(pe2, tp2, pe1)
    lyre = 0
    lyrlast = 0
    for i in range(lptr, prof.gNumLevels):
        if not QC(prof.gSndg[i][prof.tind]): continue
        pe2 = prof.gSndg[i][prof.pind]
        h2 = prof.gSndg[i][prof.zind]
        te2 = interp.vtmp(pe2, prof)
        tp2 = thermo.wetlift(pe1, tp1, pe2)
        tdef1 = (thermo.virtemp(pe1, tp1, tp1) - te1) / thermo.ctok(te1)
        tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / thermo.ctok(te2)
        lyrlast = lyre
        lyre = G * (tdef1 + tdef2) / 2. * (h2 - h1)

        # Add layer energy to total positive if lyre > 0
        if lyre > 0: totp += lyre
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
            h2 = interp.hght(pe2, prof)
            te2 = interp.vtmp(pe2, prof)
            tp2 = thermo.wetlift(pe3, tp3, pe2)
            tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / thermo.ctok(te3)
            tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / thermo.ctok(te2)
            lyrf = G * (tdef3 + tdef2) / 2. * (h2 - h3)
            if lyrf > 0: pcl.bplus += lyrf
            else:
                if pe2 > 500.: pcl.bminus += lyrf
            if pcl.bplus == 0: pcl.bminus = 0.

        # Is this the freezing level
        if te2 < 0. and not QC(pcl.bfzl):
            pe3 = pelast
            h3 = interp.hght(pe3, prof)
            te3 = interp.vtmp(pe3, prof)
            tp3 = thermo.wetlift(pe1, tp1, pe3)
            lyrf = lyre
            if lyrf > 0.: pcl.bfzl = totp - lyrf
            else: pcl.bfzl = totp
            if not QC(p0c) or p0c > pe3:
                pcl.bfzl = 0
            elif QC(pe2):
                te2 = interp.vtmp(pe2, prof)
                tp2 = thermo.wetlift(pe3, tp3, pe2)
                tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / \
                    thermo.ctok(te3)
                tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / \
                    thermo.ctok(te2)
                lyrf = G * (tdef3 + tdef2) / 2. * (hgt0c - h3)
                if lyrf > 0: pcl.bfzl += lyrf

        # Is this the -10C level
        if te2 < -10. and not QC(pcl.wm10c):
            pe3 = pelast
            h3 = interp.hght(pe3, prof)
            te3 = interp.vtmp(pe3, prof)
            tp3 = thermo.wetlift(pe1, tp1, pe3)
            lyrf = lyre
            if lyrf > 0.: pcl.wm10c = totp - lyrf
            else: pcl.wm10c = totp
            if not QC(pm10c) or pm10c > pcl.lclpres:
                pcl.wm10c = 0
            elif QC(pe2):
                te2 = interp.vtmp(pe2, prof)
                tp2 = thermo.wetlift(pe3, tp3, pe2)
                tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / \
                    thermo.ctok(te3)
                tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / \
                    thermo.ctok(te2)
                lyrf = G * (tdef3 + tdef2) / 2. * (hgtm10c - h3)
                if lyrf > 0: pcl.wm10c += lyrf

        # Is this the -20C level
        if te2 < -20. and not QC(pcl.wm20c):
            pe3 = pelast
            h3 = interp.hght(pe3, prof)
            te3 = interp.vtmp(pe3, prof)
            tp3 = thermo.wetlift(pe1, tp1, pe3)
            lyrf = lyre
            if lyrf > 0.: pcl.wm20c = totp - lyrf
            else: pcl.wm20c = totp
            if not QC(pm20c) or pm20c > pcl.lclpres:
                pcl.wm20c = 0
            elif QC(pe2):
                te2 = interp.vtmp(pe2, prof)
                tp2 = thermo.wetlift(pe3, tp3, pe2)
                tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / \
                    thermo.ctok(te3)
                tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / \
                    thermo.ctok(te2)
                lyrf = G * (tdef3 + tdef2) / 2. * (hgtm20c - h3)
                if lyrf > 0: pcl.wm20c += lyrf

        # Is this the -30C level
        if te2 < -30. and not QC(pcl.wm30c):
            pe3 = pelast
            h3 = interp.hght(pe3, prof)
            te3 = interp.vtmp(pe3, prof)
            tp3 = thermo.wetlift(pe1, tp1, pe3)
            lyrf = lyre
            if lyrf > 0.: pcl.wm30c = totp - lyrf
            else: pcl.wm30c = totp
            if not QC(pm30c) or pm30c > pcl.lclpres:
                pcl.wm30c = 0
            elif QC(pe2):
                te2 = interp.vtmp(pe2, prof)
                tp2 = thermo.wetlift(pe3, tp3, pe2)
                tdef3 = (thermo.virtemp(pe3, tp3, tp3) - te3) / \
                    thermo.ctok(te3)
                tdef2 = (thermo.virtemp(pe2, tp2, tp2) - te2) / \
                    thermo.ctok(te2)
                lyrf = G * (tdef3 + tdef2) / 2. * (hgtm30c - h3)
                if lyrf > 0: pcl.wm30c += lyrf

        # Is this the 3km level
        if pcl.lclhght < 3000.:
            h = interp.agl(interp.hght(pe2, prof), prof)
            if h >= 3000. and not QC(pcl.b3km):
                pe3 = pelast
                h3 = interp.hght(pe3, prof)
                te3 = interp.vtmp(pe3, prof)
                tp3 = thermo.wetlift(pe1, tp1, pe3)
                lyrf = lyre
                if lyrf > 0: pcl.b3km = totp - lyrf
                else: pcl.b3km = totp
                h2 = interp.msl(3000., prof)
                pe2 = interp.pres(h2, prof)
                if QC(pe2):
                    te2 = interp.vtmp(pe2, prof)
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
            h = interp.agl(interp.hght(pe2, prof), prof)
            if h >= 6000. and not QC(pcl.b6km):
                pe3 = pelast
                h3 = interp.hght(pe3, prof)
                te3 = interp.vtmp(pe3, prof)
                tp3 = thermo.wetlift(pe1, tp1, pe3)
                lyrf = lyre
                if lyrf > 0: pcl.b6km = totp - lyrf
                else: pcl.b6km = totp
                h2 = interp.msl(6000., prof)
                pe2 = interp.pres(h2, prof)
                if QC(pe2):
                    te2 = interp.vtmp(pe2, prof)
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
            while interp.vtmp(pe3, prof) > thermo.virtemp(pe3,
                thermo.wetlift(pe2, tp3, pe3), thermo.wetlift(pe2, tp3, pe3)):
                    pe3 -= 5
            pcl.lfcpres = pe3
            pcl.lfchght = interp.agl(interp.hght(pe3, prof), prof)
            cinh_old = totn
            tote = 0.
            pcl.elpres = RMISSD
            li_max = RMISSD

            if cap_strength < 0.: cap_strength = 0.
            pcl.cap = cap_strength
            pcl.cappres = cap_strengthpres
            # Hack to force LFC to be at least at the LCL
            if pcl.lfcpres > pcl.lclpres:
                pcl.lfcpres = pcl.lclpres
                pcl.lfchght = pcl.lclhght

        # EL Possibility
        if lyre <= 0. and lyrlast >= 0.:
            tp3 = tp1
            te3 = te1
            pe2 = pe1
            pe3 = pelast
            while interp.vtmp(pe3, prof) < thermo.virtemp(pe3,
                thermo.wetlift(pe2, tp3, pe3), thermo.wetlift(pe2, tp3, pe3)):
                    pe3 -= 5
            pcl.elpres = pe3
            pcl.elhght = interp.agl(interp.hght(pe3, prof), prof)
            pcl.mplpres = RMISSD
            pcl.limax = -li_max
            pcl.limaxpress = li_maxpres

        # MPL Possibility
        if tote < 0. and not QC(pcl.mplpres) and QC(pcl.elpres):
            pe3 = pelast
            h3 = interp.hght(pe3, prof)
            te3 = interp.vtmp(pe3, prof)
            tp3 = thermo.wetlift(pe1, tp1, pe3)
            totx = tote - lyre
            pe2 = pelast
            while totx > 0:
                pe2 -= 1
                te2 = interp.vtmp(pe2, prof)
                tp2 = thermo.wetlift(pe3, tp3, pe2)
                h2 = interp.hght(pe2, prof)
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
            pcl.mplhght = interp.agl(interp.hght(pe2, prof), prof)

        # 500 hPa Lifted Index
        if prof.gSndg[i][prof.pind] <= 500. and pcl.li5 == RMISSD:
            a = interp.vtmp(500., prof)
            b = thermo.wetlift(pe1, tp1, 500.)
            pcl.li5 = a - thermo.virtemp(500, b, b)

        # 300 hPa Lifted Index
        if prof.gSndg[i][prof.pind] <= 300. and pcl.li3 == RMISSD:
            a = interp.vtmp(300., prof)
            b = thermo.wetlift(pe1, tp1, 300.)
            pcl.li3 = a - thermo.virtemp(300, b, b)

    # Calculate BRN if available
    pcl = bulk_rich(pcl, prof)

    pcl.bminus = cinh_old
    if pcl.bplus == 0: pcl.bminus = 0.
    return pcl


def temp_lvl(temp, prof):
    '''
    Calculates the level (hPa) of the first occurrence of the specified
    temperature.

    Inputs
    ------
        temp        (float)             Temperature being searched (C)
        prof        (profile object)    Profile Object

    Returns
    -------
        Level of the temperature (hPa)
    '''
    for i in range(prof.gNumLevels):
        if QC(prof.gSndg[i][prof.tind]) and prof.gSndg[i][prof.tind] <= temp:
            if i == 0: return RMISSD
            if prof.gSndg[i][prof.tind] == temp:
                return prof.gSndg[i][prof.pind]
            p0 = prof.gSndg[i-1][prof.pind]
            t0 = prof.gSndg[i-1][prof.tind]
            nm1 = temp - t0
            nm2 = prof.gSndg[i][prof.tind] - t0
            nm3 = math.log(prof.gSndg[i][prof.pind] / p0)
            return p0 * math.exp((nm1 / nm2) * nm3)
    return RMISSD


def bulk_rich(pcl, prof):
    '''
    Calculates the Bulk Richardson Number for a given parcel.

    Inputs
    ------
        pcl         (parcel object)         Parcel Object
        prof        (profile object)        Profile Object

    Returns
    -------
        Bulk Richardson Number
    '''
    # Make sure parcel is initialized
    if pcl.lplvals.flag == RMISSD:
        pbot = RMISSD
    elif pcl.lplvals.flag > 0 and pcl.lplvals.flag < 4:
        ptop = interp.pres(interp.msl(6000., prof), prof)
        pbot = prof.gSndg[prof.sfc][prof.pind]
    else:
        h0 = interp.hght(pcl.pres, prof)
        try:
            pbot = interp.pres(h0 - 500., prof)
        except:
            pbot = RMISSD
        if not QC(pbot): pbot = prof.gSndg[prof.sfc][prof.pind]
        h1 = interp.hght(pbot, prof)
        ptop = interp.pres(h1 + 6000., prof)

    if not QC(pbot) or not QC(ptop):
        pcl.brnshear = RMISSD
        pcl.brn = RMISSD
        return pcl

    # Calculate lowest 500m mean wind
    p = interp.pres(interp.hght(pbot, prof) + 500., prof)
    mnlu, mnlv = winds.mean_wind(pbot, p, prof)

    # Calculate the 6000m mean wind
    mnuu, mnuv = winds.mean_wind(pbot, ptop, prof)

    # Make sure CAPE and Shear are available
    if not QC(pcl.bplus) or not QC(mnlu) or not QC(mnuu):
        pcl.brnshear = RMISSD
        pcl.brn = RMISSD
        return pcl

    # Calculate shear between levels
    dx = mnuu - mnlu
    dy = mnuv - mnlv

    pcl.brnshear = KTS2MS(vector.comp2vec(dx, dy)[1])
    pcl.brnshear = pcl.brnshear**2 / 2.
    pcl.brn = pcl.bplus / pcl.brnshear
    return pcl


def max_temp(prof, mixlyr=-1):
    '''
    Calculates a maximum temperature forecast based on the depth of the mixing
    layer and low-level temperatures

    Inputs
    ------
        prof        (profile object)    Profile Object
        mixlyr      (float)             Top of layer over which to "mix" (hPa)

    Returns
    -------
        mtemp       (float)             Forecast Maximum Temperature
    '''
    sfcpres = prof.gSndg[prof.sfc][prof.pind]
    if mixlyr == -1: mixlyr = sfcpres - 100.
    temp = thermo.ctok(interp.temp(mixlyr, prof)) + 2.
    return thermo.ktoc(temp * (sfcpres / mixlyr)**ROCP)


def mean_mixratio(prof, lower=-1, upper=-1):
    '''
    Calculates the mean mixing ratio from a profile object within the
    specified layer.

    Inputs
    ------
        prof        (profile object)    Profile Object
        lower       (float)             Bottom level (hPa) [-1=SFC]
        upper       (float)             Top level (hPa) [-1=SFC-100hPa]

    Returns
    -------
        Mean Mixing Ratio   (float)
    '''
    if lower == -1: lower = prof.gSndg[prof.sfc][prof.pind]
    if upper == -1: upper = prof.gSndg[prof.sfc][prof.pind] - 100.

    if not QC(interp.temp(upper, prof)): mmw = RMISSD
    if not QC(interp.temp(lower, prof)): prof.gSndg[prof.sfc][prof.pind]

    # Find lowest observations in the layer
    i = 0
    while prof.gSndg[i][prof.pind] > lower: i+=1
    while not QC(prof.gSndg[i][prof.tdind]): i+=1
    lptr = i
    if prof.gSndg[i][prof.pind] == lower: lptr+=1

    # Find highest observations in the layer
    i = prof.gNumLevels - 1
    while prof.gSndg[i][prof.pind] < upper: i-=1
    uptr = i
    if prof.gSndg[i][prof.pind] == upper: uptr-=1

    totd = 0
    totp = 0

    # Start with interpolated bottom layer
    p1 = lower
    dp1 = interp.dwpt(p1, prof)
    num = 1

    # Calculate every level that reports a dew point
    for i in range(lptr, uptr+1):
        if QC(prof.gSndg[i][prof.tdind]):
            dp2 = prof.gSndg[i][prof.tdind]
            p2 = prof.gSndg[i][prof.pind]
            dpbar = (dp1 + dp2) / 2.
            pbar = (p1 + p2) / 2.
            totd += dpbar
            totp += pbar
            dp1 = dp2
            p1 = p2
            num += 1

    # Finish with top layer
    dp2 = interp.dwpt(upper, prof)
    p2 = upper
    dbar = (dp1 + dp2) / 2.
    pbar = (p1 + p2) / 2.
    totd += dbar
    totp += pbar
    return thermo.mixratio(totp/num, totd/num)


def mean_theta(prof, lower=-1, upper=-1):
    '''
    Calculates the mean theta from a profile object within the
    specified layer.

    Inputs
    ------
        prof        (profile object)    Profile Object
        lower       (float)             Bottom level (hPa) [-1=SFC]
        upper       (float)             Top level (hPa) [-1=SFC-100hPa]

    Returns
    -------
        Mean Theta   (float)
    '''
    if lower == -1: lower = prof.gSndg[prof.sfc][prof.pind]
    if upper == -1: upper = prof.gSndg[prof.sfc][prof.pind] - 100.

    if not QC(interp.temp(upper, prof)): mmw = RMISSD
    if not QC(interp.temp(lower, prof)): prof.gSndg[prof.sfc][prof.pind]

    # Find lowest observations in the layer
    i = 0
    while prof.gSndg[i][prof.pind] > lower: i+=1
    while not QC(prof.gSndg[i][prof.tind]): i+=1
    lptr = i
    if prof.gSndg[i][prof.pind] == lower: lptr+=1

    # Find highest observations in the layer
    i = prof.gNumLevels - 1
    while prof.gSndg[i][prof.pind] < upper: i-=1
    uptr = i
    if prof.gSndg[i][prof.pind] == upper: uptr-=1

    tott = 0

    # Start with interpolated bottom layer
    t1 = thermo.theta(lower, interp.temp(lower, prof), 1000.)
    num = 1

    # Calculate every level that reports a dew point
    for i in range(lptr, uptr+1):
        if QC(prof.gSndg[i][prof.tind]):
            t2 = thermo.theta(prof.gSndg[i][prof.pind],
                prof.gSndg[i][prof.tind], 1000.)
            tbar = (t1 + t2) / 2.
            tott += tbar
            t1 = t2
            num += 1

    # Finish with top layer
    t2 = thermo.theta(upper, interp.temp(upper, prof), 1000.)
    tbar = (t1 + t2) / 2.
    tott += tbar
    return tott / num


def unstable_level(prof, lower, upper):
    '''
    Finds the most unstable level between the lower and upper levels.

    Inputs
    ------
        prof        (profile object)    Profile Object
        lower       (float)             Bottom level (hPa) [-1=SFC]
        upper       (float)             Top level (hPa) [-1=SFC-100hPa]

    Returns
    -------
        Pressure Level of most unstable level   (float [hPa])
    '''
    if lower == -1: lower = prof.gSndg[prof.sfc][prof.pind]
    if upper == -1: upper = prof.gSndg[prof.sfc][prof.pind] - 400.

    # Make sure this is a valid layer
    while not QC(interp.dwpt(upper, prof)): upper += 50.
    if not QC(interp.temp(lower, prof)): lower = prof.gSndg[prof.sfc][0]

    # Find lowest observations in the layer
    i = 0
    while prof.gSndg[i][prof.pind] > lower: i+=1
    while not QC(prof.gSndg[i][prof.tind]): i+=1
    lptr = i
    if prof.gSndg[i][prof.pind] == lower: lptr+=1

    # Find highest observations in the layer
    i = prof.gNumLevels - 1
    while prof.gSndg[i][prof.pind] < upper: i-=1
    uptr = i
    if prof.gSndg[i][prof.pind] == upper: uptr-=1

    # Start with interpolated bottom layer
    p1 = lower
    t1 = interp.temp(p1, prof)
    td1 = interp.dwpt(p1, prof)
    p2, t2 = thermo.drylift(p1, t1, td1)
    tmax = thermo.wetlift(p2, t2, 1000.)
    pmax = p1

    # Calculate every level that reports a dew point
    for i in range(lptr, uptr+1):
        if QC(prof.gSndg[i][prof.tdind]):
            p1 = prof.gSndg[i][prof.pind]
            t1 = prof.gSndg[i][prof.tind]
            td1 = prof.gSndg[i][prof.tdind]
            p2, t2 = thermo.drylift(p1, t1, td1)
            t1 = thermo.wetlift(p2, t2, 1000.)
            if t1 > tmax:
                tmax = t1
                pmax = p1

    # Finish with interpolated top layer
    p1 = upper
    t1 = interp.temp(p1, prof)
    td1 = interp.dwpt(p1, prof)
    p2, t2 = thermo.drylift(p1, t1, td1)
    t1 = thermo.wetlift(p2, t2, 1000.)
    if t1 > tmax:
        pmax = prof.gSndg[i][prof.pind]

    return pmax


def effective_inflow_layer(ecape, ecinh, prof, **kwargs):
    '''
    Calculates the top and bottom of the effective inflow layer based on
    research by Thompson et al. (2004)

    Inputs
    ------
        ecape       (float)                 CAPE Threshold
        ecinh       (float)                 CINH Threshold
        prof        (profile object)        Profile Object
    Returns
    -------
        pbot        (float)                 Pressure at bottom level (hPa)
        ptop        (float)                 Pressure at top level (hPa)
    '''
    if 'mupcl' in kwargs:
        mupcl = kwargs.get('mupcl')
    else:
        mulplvals = DefineParcel(prof, 3, pres=400)
        mupcl = parcelx(-1, -1, mulplvals.pres, mulplvals.temp,
            mulplvals.dwpt, prof, lplvals=mulplvals)
    mucape = mupcl.bplus
    mucinh = mupcl.bminus

    # Scenario where shallow buoyancy present for
    # parcel with lesser theta near the ground
    mu2lplvals = DefineParcel(prof, 3, pres=300)
    mu2pcl = parcelx(-1, -1, mulplvals.pres, mulplvals.temp, mulplvals.dwpt,
        prof, lplvals=mu2lplvals)
    if mu2pcl.bplus > mucape:
        mucape = mu2pcl.bplus
        mucinh = mu2pcl.bminus

    pbot = RMISSD
    ptop = RMISSD

    if mucape >= ecape and mucinh > ecinh:
        # Begin at surface and search upward for effective surface
        for i in range(prof.sfc, prof.gNumLevels-1):
            pcl = parcelx(-1, -1, prof.gSndg[i][prof.pind],
                prof.gSndg[i][prof.tind], prof.gSndg[i][prof.tdind], prof)
            if pcl.bplus >= ecape and pcl.bminus >= ecinh:
                pbot = prof.gSndg[i][prof.pind]
                break
        if pbot == RMISSD: return pbot, ptop
        bptr = i
        # Keep searching upward for the effective top
        for i in range(bptr+1, prof.gNumLevels-1):
            pcl = parcelx(-1, -1, prof.gSndg[i][prof.pind],
                prof.gSndg[i][prof.tind], prof.gSndg[i][prof.tdind], prof)
            if pcl.bplus <= ecape or pcl.bminus <= ecinh:
                j = 1
                while not QC(prof.gSndg[i-j][prof.tind]) and \
                      not QC(prof.gSndg[i-j][prof.tdind]): j+=1
                ptop = prof.gSndg[i-j][prof.pind]
                if ptop > pbot: ptop = pbot
                break

    return pbot, ptop


def bunkers_storm_motion(prof, pbot=None, **kwargs):
    '''
    Compute the Bunkers Storm Motion for a Right Moving Supercell using
    a parcel based approach.

    Inputs
    ------
        prof        (profile object)    Profile Object
        pbot        (float)             Base of effective-inflow layer (hPa)

    Returns
    -------
        rstu        (float)             Right Storm Motion U-component
        rstv        (float)             Right Storm Motion V-component
        lstu        (float)             Left Storm Motion U-component
        lstv        (float)             Left Storm Motion V-component
    '''
    d = MS2KTS(7.5)     # Deviation value emperically derived as 7.5 m/s

    # If MUPCL provided, use it, otherwise create MUPCL
    if 'mupcl' in kwargs:
        mupcl = kwargs.get('mupcl')
    else:
        mulplvals = params.DefineParcel(3, prof, pres=400)
        mupcl = params.parcelx(-1, -1, mulplvals.pres, mulplvals.temp,
            mulplvals.dwpt, prof, lplvals=mulplvals)

    mucape = mupcl.bplus
    mucinh = mupcl.bminus
    muel = mupcl.elhght
    if not pbot:
        pbot, ptop = effective_inflow_layer(100, -250, prof)
    base = interp.agl(interp.hght(pbot, prof), prof)
    if mucape > 100. and QC(muel) and base >= 750:
        depth = muel - base
        htop = base + depth / 2.
        ptop = interp.pres(interp.msl(base + htop, prof), prof)
        mnu, mnv = winds.mean_wind_npw(pbot, ptop, prof)
        sru, srv = winds.wind_shear(pbot, ptop, prof)
        srmag = vector.mag(sru, srv)
        uchg = d / srmag * srv
        vchg = d / srmag * sru
        rstu = mnu + uchg
        rstv = mnv - vchg
        lstu = mnu - uchg
        lstv = mnv + vchg
    else:
        rstu, rstv, lstu, lstv =  winds.non_parcel_bunkers_motion(prof)

    return rstu, rstv, lstu, lstv


def convective_temp(prof, mincinh=-1.):
    '''
    Computes the convective temperature, assuming no change in the moisture
    profile. Parcels are iteratively lifted until only mincinh is left as a
    cap. The first guess is the observed surface temperature.

    Inputs
    ------
        prof        (profile object)    Profile Object
        mincinh     (float)             Amount of CINH left at CI

    Returns
    -------
        Convective Temperature (float [C])
    '''
    mmr = mean_mixratio(prof, -1, -1)
    sfcpres = prof.gSndg[prof.sfc][prof.pind]
    sfctemp = prof.gSndg[prof.sfc][prof.tind]
    sfcdwpt = thermo.temp_at_mixrat(mmr, sfcpres)

    # Do a quick search to find wheather to continue. If
    # If you need to heat up more than 25C, don't compute.
    pcl = parcelx(-1, -1, sfcpres, sfctemp+25., sfcdwpt, prof)
    if pcl.bplus == 0. or pcl.bminus < mincinh: return RMISSD

    excess = sfcdwpt - sfctemp
    if excess > 0: sfctemp = sfctemp + excess + 4.
    pcl = parcelx(-1, -1, sfcpres, sfctemp, sfcdwpt, prof)
    if pcl.bplus == 0.: pcl.bminus = RMISSD
    while pcl.bminus < mincinh:
        if pcl.bminus < -100.: sfctemp += 2.
        else: sfctemp += 0.5
        pcl = parcelx(-1, -1, sfcpres, sfctemp, sfcdwpt, prof)
        if pcl.bplus == 0.: pcl.bminus = RMISSD

    return sfctemp


def esfc(prof, val=50.):
    '''
    Calculates the effective surface for elevated convection. Assumes that
    the lowest layer with CAPE >= val is the "surface".

    Inputs
    ------
        prof        (profile object)    Profile Object
        val         (float)             Amount of CAPE necessary for CI

    Returns
    -------
        Level of the effective surface (float [hPa])
    '''
    for i in range(prof.sfc, prof.gNumLevels):
        if prof.gSndg[prof.sfc][prof.pind] < 500.: break
        pcl = parcelx(-1, -1, prof.gSndg[i][prof.pind],
            prof.gSndg[i][prof.tind], prof.gSndg[i][prof.tdind], prof)
        if pcl.bplus >= val: return prof.gSndg[i][prof.pind]
    return RMISSD


def lapse_rate(lower, upper, prof, pres=1):
    '''
    Calculates the lapse rate (C/km) from a profile object

    Inputs
    ------
        lower       (float)             Lower Bound of lapse rate
        upper       (float)             Upper Bound of lapse rate
        prof        (profile object)    Profile Object
        pres        (int 0/1)           Flag to know to convert pres to height

    Returns
    -------
        lapse rate  (float [C/km])
    '''
    if pres:
        p1 = lower
        p2 = upper
        z1 = interp.hght(lower, prof)
        z2 = interp.hght(upper, prof)
    else:
        z1 = interp.msl(lower, prof)
        z2 = interp.msl(upper, prof)
        p1 = interp.pres(z1, prof)
        p2 = interp.pres(z2, prof)

    tv1 = interp.vtmp(p1, prof)
    tv2 = interp.vtmp(p2, prof)
    if not QC(tv1) or not QC(tv2): return RMISSD
    return (tv2 - tv1) / (z2 - z1) * -1000.






