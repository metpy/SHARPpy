''' Vector Manipulation Routines '''
import math
from sharppy.sharptab.constants import *

__all__ = ['vec2comp', 'comp2vec', 'mag']


def vec2comp(dir, spd):
    '''
    Convert direction and speed into U,V components

    Inputs
    ------
        dir    (float)         direction (degrees)
        spd    (float)         speed

    Returns
    -------
        u       (float)         U-component
        v       (float)         V-component
    '''
    if not QC(dir) or not QC(spd): return RMISSD, RMISSD
    u = spd * math.sin(math.radians(dir % 360.)) * -1
    v = spd * math.cos(math.radians(dir % 360.)) * -1
    return u, v


def comp2vec(u, v):
    '''
    Convert U,V components into direction and speed

    Inputs
    ------
        u       (float)         U-component
        v       (float)         V-component

    Returns
    -------
        dir    (float)         direction (degrees)
        spd    (float)         speed
    '''
    if not QC(u) or not QC(v): return RMISSD, RMISSD
    dir =  math.degrees(math.atan2(-u, -v))
    if dir < 0: dir += 360
    spd = mag(u, v)
    return dir, spd


def mag(u, v):
    '''
    Compute the magnitude of a vector

    Inputs
    ------
        u       (float)         U-component
        v       (float)         V-component

    Returns
    -------
        Returns the magnitude of a vector (float)
    '''
    if not QC(u) or not QC(v): return RMISSD
    return (u**2 + v**2)**0.5

