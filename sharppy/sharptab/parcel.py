''' Define Various Parcels '''
import math
from sharppy.sharptab import *
from sharppy.sharptab.constants import *

__all__ = ['DefineParcel', 'Parcel']


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

    def __init__(self, profile, flag, pres=None):
        if not pres:
            pres = profile.gSndg[profile.sfc][0]
        self.flag = flag
        self.presval = pres
        if flag == 1:
            self.__sfc(profile)
        elif flag == 2:
            self.__fcst(profile)
        elif flag == 3:
            self.__mu(profile)
        elif flag == 4:
            self.__ml(profile)
        elif flag == 5:
            self.__user(profile)
        elif flag == 6:
            self.__effective(profile)
        else:
            print 'Defaulting to Surface Parcel'
            self.__sfc()

        print self.desc


    def __sfc(self, profile):
        ''' Create a parcel using surface conditions '''
        self.desc = 'Surface Parcel'
        self.temp = profile.gSndg[profile.sfc][2]
        self.dwpt = profile.gSndg[profile.sfc][3]
        self.pres = profile.gSndg[profile.sfc][0]

    def __fcst(self, profile):
        ''' Create a parcel using forecast conditions '''
        self.desc = 'Forecast Surface Parcel'

    def __mu(self, profile):
        ''' Create the most unstable parcel within defined level '''
        self.desc = 'Most Unstable Parcel in Lowest %.2fmb' % self.presval

    def __ml(self, profile):
        ''' Create the mixed-layer parcel; mixing over defined pressure '''
        self.desc = '%.2fmb Mixed Layer Parcel' % self.presval

    def __user(self, profile):
        ''' Create a user-defined parcel '''
        self.desc = '%.2fmb Parcel' % self.presval

    def __effective(self):
        ''' Create the mean-effective layer parcel '''
        self.desc = 'Mean Effective Layer'


class Parcel(object):
    ''' Initialize variables for a parcel '''

    def __init__(self, lower, upper, lplvals, missing=RMISSD.):
        self.lplvals = lplvals
        self.pres = lplvals.pres
        self.temp = lplvals.temp
        self.dwpt = lplvals.dwpt
        self.blayer = lower
        self.tlayer = upper
        self.entrain = 0.
        self.lclpres = RMISSD
        self.lclhght = RMISSD
        self.lfcpres = RMISSD
        self.elpres = RMISSD
        self.mplpres = RMISSD
        self.bplus = RMISSD
        self.bminus = RMISSD
        self.bfzl = RMISSD
        self.cape3km = RMISSD
        self.wm10c = RMISSD
        self.li5 = RMISSD
        self.li3 = RMISSD
        self.brnshear = RMISSD
        self.brn = RMISSD
        self.limax = RMISSD
        self.limaxpres = RMISSD
        self.cap = RMISSD
        self.cappres = RMISSD













