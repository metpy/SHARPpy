''' Thermodynamic Library '''
import math
from sharppy.sharptab.qc import qc
from sharppy.sharptab.constants import *


__all__ = ['lifted', 'drylift', 'lcltemp', 'thalvl', 'theta', 'wetlift',
           'wobf', 'satlift', 'temp_at_mixrat', 'mixratio', 'vappres',
           'wetbulb', 'thetaw', 'thetae', 'virtemp', 'esfc', 'c2f', 'relh',
           'effective_inflow_layer']


def lifted(p, t, td, lev):
    '''
    Calculated a lifted index for a parcel at a level.

    Inputs
    ------
        p       (float)         Pressure of initial parcel in hPa
        t       (float)         Temperature of initial parcel in C
        td      (float)         Dew Point of initial parcel in C
        lev     (float)         Pressure to which parcel is lifted in hPa

    Returns
    -------
        Temperature (C [float]) of lifted parcel
    '''
    if not qc(p) or not qc(t) or not qc(td): return RMISSD
    p2, t2 = drylift(p, t, td)
    return wetlift(p2, t2, lev)


def drylift(p, t, td):
    '''
    Lifts a parcel to the LCL and returns its new level and temperature.

    Inputs
    ------
        p       (float)         Pressure of initial parcel in hPa
        t       (float)         Temperature of inital parcel in C
        td      (float)         Dew Point of initial parcel in C

    Returns
    -------
        p2      (float)         LCL pressure in hPa
        t2      (float)         LCL Temperature in C
    '''
    p2 = RMISSD
    t2 = RMISSD
    if qc(p) or qc(t) or qc(td):
        t2 = lcltemp(t, td)
        p2 = thalvl(theta(p, t, 1000.), t2)
    return p2, t2


def lcltemp(t, td):
    '''
    Returns the temperature (C) of a parcel when raised to its LCL.

    Inputs
    ------
        t       (float)         Temperature of the parcel (C)
        td      (float)         Dew Point of the parcel (C)

    Returns
    -------
        Temperature (C [float]) of the parcel at its LCL
    '''
    if not qc(t) or not qc(td): return RMISSD
    s = t - td
    dlt = s * (1.2185 + 0.001278 * t + s * (-0.00219 + 1.173e-5 * s -
        0.0000052 * t))
    return t - dlt


def thalvl(thta, t):
    '''
    Returns the level (hPa) of a parcel.

    Inputs
    ------
        thta    (float)         Potential temperature of the parcel (C)
        t       (float)         Temperature of the parcel (C)

    Returns
    -------
        Pressure Level (hPa [float]) of the parcel
    '''
    if not qc(t) or not qc(thta): return RMISSD
    t += ZEROCNK
    thta += ZEROCNK
    return 1000. / ((thta / t)**(1/ROCP))


def theta(p, t, p2):
    '''
    Returns the potential temperature (C) of the parcel

    Inputs
    ------
        p       (float)         Pressure of the parcel (hPa)
        t       (float)         Temperature of the parcel (C)
        p2      (float)         Reference Pressure Level (hPa; Usually 1000.)

    Returns
    -------
        Potential Temperature (C [float])
    '''
    if not qc(p) or not qc(t) or not qc(p2): return RMISSD
    t += ZEROCNK
    return (t * (p2 / p)**ROCP) - ZEROCNK


def wetlift(p, t, p2):
    '''
    Lifts a parcel moist adiabatically to its new level.

    Inputs
    ------
        p       (float)         Pressure of initial parcel (hPa)
        t       (float)         Temperature of initial parcel (C)
        p2      (float)         Pressure of final level (hPa)

    Returns
    -------
        Temperature (C [float])
    '''
    if not qc(p) or not qc(t) or not qc(p2): return RMISSD
    thta = theta(p, t, 1000.)
    thm = thta - wobf(thta) + wobf(t)
    return satlift(p2, thm)


def wobf(t):
    '''
    Implementation of the Wobus Function for computing the moist adiabats.

    Inputs
    ------
        t       (float)         Temperature (C)

    Returns
    -------
        Correction to thta for calculation of saturated potential
        temperature
    '''
    if not qc(t): return RMISSD
    x = t - 20.
    if x <= 0:
        pol = 1 + x * (-8.841660499999999e-3 + x * ( 1.4714143e-4
              + x * (-9.671989000000001e-7 + x * (-3.2607217e-8
              + x * (-3.8598073e-10)))))
        return 15.13 / (pol**4)
    else:
        pol = x * (4.9618922e-07 + x * (-6.1059365e-09 +
              x * (3.9401551e-11 + x * (-1.2588129e-13 +
              x * (1.6688280e-16)))))
        pol = 1 + x * (3.6182989e-03 + x * (-1.3603273e-05 + pol))
        return (29.93 / (pol**4)) + (0.96 * x) - 14.8


def satlift(p, thm):
    '''
    Returns the temperature (C) of a saturated parcel (thm) when lifted to a
    new pressure level (hPa)

    Inputs
    ------
        p       (float)         Pressure to which parcel is raised (hPa)
        thm     (float)         Saturated Potential Temperature of parcel (C)

    Returns
    -------
        Temperature (C [float]) of saturated parcel at new level
    '''
    if not qc(p) or not qc(thm): return RMISSD
    if math.fabs(p - 1000.) - 0.001 <= 0: return thm
    eor = 999
    while math.fabs(eor) - 0.1 > 0:
        if eor == 999:                  # First Pass
            pwrp = (p / 1000.)**ROCP
            t1 = (thm + ZEROCNK) * pwrp - ZEROCNK
            e1 = wobf(t1) - wobf(thm)
            rate = 1
        else:                           # Successive Passes
            rate = (t2 - t1) / (e2 - e1)
            t1 = t2
            e1 = e2
        t2 = t1 - (e1 * rate)
        e2 = (t2 + ZEROCNK) / pwrp - ZEROCNK
        e2 += wobf(t2) - wobf(e2) - thm
        eor = e2 * rate
    return t2 - eor


def temp_at_mixrat(w, p):
    '''
    Returns the temperature (C) of air at the given mixing ratio (g/kg) and
    pressure (hPa)

    Inputs
    ------
        w       (float)         Mixing Ratio (g/kg)
        p       (float)         Pressure (hPa)

    Returns
    -------
        Temperature (C [float]) of air at given mixing ratio and pressure
    '''
    if not qc(w) or not qc(p): return RMISSD
    c1 = 0.0498646455
    c2 = 2.4082965
    c3 = 7.07475
    c4 = 38.9114
    c5 = 0.0915
    c6 = 1.2035
    x = math.log10(w * p / (622. + w))
    return (10.**((c1 * x) + c2) - c3 + (c4 * (10**(c5 * x) - c6)**2)) - ZEROCNK


def mixratio(p, t):
    '''
    Returns the mixing ratio (g/kg) of a parcel

    Inputs
    ------
        p       (float)         Pressure of the parcel (hPa)
        t       (float)         Temperature of the parcel (hPa)

    Returns
    -------
        Mixing Ratio (g/kg) of the given parcel
    '''
    if not qc(p) or not qc(t): return RMISSD
    x = 0.02 * (t - 12.5 + (7500. / p))
    wfw = 1. + (0.0000045 * p) + (0.0014 * x * x)
    fwesw = wfw * vappres(t)
    return 621.97 * (fwesw / (p - fwesw))


def vappres(t):
    '''
    Returns the vapor pressure of dry air at given temperature

    Inputs
    ------
        t       (float)         Temperature of the parcel (C)

    Returns
    -------
        Vapor Pressure of dry air
    '''
    if not qc(t): RMISSD
    pol = t * (1.1112018e-17 + (t * -3.0994571e-20))
    pol = t * (2.1874425e-13 + (t * (-1.789232e-15 + pol)))
    pol = t * (4.3884180e-09 + (t * (-2.988388e-11 + pol)))
    pol = t * (7.8736169e-05 + (t * (-6.111796e-07 + pol)))
    pol = 0.99999683 + (t * (-9.082695e-03 + pol))
    return 6.1078 / pol**8


def wetbulb(p, t, td):
    '''
    Calculates the wetbulb temperature (C) for the given parcel

    Inputs
    ------
        p       (float)         Pressure of parcel (hPa)
        t       (float)         Temperature of parcel (C)
        td      (float)         Dew Point of parcel (C)

    Returns
    -------
        Wetbulb temperature (C [float])
    '''
    if not qc(p) or not qc(t) or not qc(td): return RMISSD
    p2, t2 = drylift(p, t, td)
    return wetlift(p2, t2, p)


def thetaw(p, t, td):
    '''
    Calculates the wetbulb potential temperature for a given parcel

    Inputs
    ------
        p       (float)         Pressure of parcel (hPa)
        t       (float)         Temperature of parcel (C)
        td      (float)         Dew Point of parcel (C)

    Returns
    -------
        Wetbulb potential temperature (C [float])
    '''
    if not qc(p) or not qc(t) or not qc(td): return RMISSD
    p2, t2 = drylift(p, t, td)
    return wetlift(p2, t2, 1000.)


def thetae(p, t, td):
    '''
    Calculates the equivalent potential temperature for a given parcel

    Inputs
    ------
        p       (float)         Pressure of parcel (hPa)
        t       (float)         Temperature of parcel (C)
        td      (float)         Dew Point of parcel (C)

    Returns
    -------
        Equivalent potential temperature (C [float])
    '''
    if not qc(p) or not qc(t) or not qc(td): RMISSD
    p2, t2 = drylift(p, t, td)
    return theta(100., wetlift(p2, t2, 100.), 1000.)


def virtemp(p, t, td):
    '''
    Calculates the virtual temperature for a given parcel

    Inputs
    ------
        p       (float)         Pressure of parcel (hPa)
        t       (float)         Temperature of parcel (C)
        td      (float)         Dew Point of parcel (C)

    Returns
    -------
        Virtual temperature (C [float])
    '''
    if not qc(td): return t
    if not qc(p) or not qc(t): RMISSD
    eps = 0.62197
    tk = t + ZEROCNK
    w = 0.001 * mixratio(p, td)
    return (tk * (1. + w / eps) / (1. + w)) - ZEROCNK


def esfc():
    pass


def effective_inflow_layer():
    pass


def relh():
    pass


def ctof(t):
    '''
    Convert from Celsius to Fahrenheit

    Inputs
    ------
        t       (float)         Temperature (C)

    Returns
    -------
        Temperature (F [float])
    '''
    if (not qc(t)): return RMISSD
    return ((1.8 * t) + 32.0)







