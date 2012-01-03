''' Create a Profile Object '''
from sharppy.sharptab import interp, vector
from sharppy.sharptab.constants import *



class Profile(object):

    def __init__(self, **kwargs):
        if 'url' in kwargs:
            url = kwargs.get('url')
            print url
            import urllib

            snfile = urllib.urlopen(url).read().decode('utf-8').split("\n")
            for i in range(0, len(snfile)):
                if (snfile[i] == "%RAW%"): bgn = i+1
                if (snfile[i] == "%END%"): end = i-1
                if (snfile[i] == "%TITLE%"): ttl = i+1

            self.gSndg = []
            for i in range(bgn, end+1):
                vals = snfile[i].split(",")
                for j in range(0, len(vals)):
                    vals[j] = float(vals[j])
                self.gSndg.append(vals)

            self.gStation = snfile[ttl][0:6]
            self.gDate = snfile[ttl][7:]
            self.gNumLevels = end-bgn+1
        else:
            pres = kwargs.get('pres')
            hght = kwargs.get('hght')
            temp = kwargs.get('temp')
            dwpt = kwargs.get('dwpt')
            if 'wdir' in kwargs:
                wnd1 = kwargs.get('wdir')
                wnd2 = kwargs.get('wspd')
            else:
                wnd1 = kwargs.get('ucomp')
                wnd2 = kwargs.get('vcomp')

            self.gSndg = []
            for i in range(len(pres)):
                vals = [pres[i], hght[i], temp[i], dwpt[i], wnd1[i], wnd2[i]]
                for j in range(0, len(vals)):
                    vals[j] = float(vals[j])
                self.gSndg.append(vals)
            self.gStation = kwargs.get('stn', '????')
            self.gDate = kwargs.get('date', '-----')
            self.gNumLevels = len(pres)

        # Set Various Indices
        self.pind = kwargs.get('pind', 0)
        self.zind = kwargs.get('zind', 1)
        self.tind = kwargs.get('tind', 2)
        self.tdind = kwargs.get('tdind', 3)
        self.uind = kwargs.get('uind', 4)
        self.vind = kwargs.get('vind', 5)
        if 'uind' not in kwargs:
            self.wdirind = kwargs.get('wdirind', 4)
            self.wspdind = kwargs.get('wspdind', 5)
            self.dir2Comp()
        self.sfc = self.get_sfc()


    def get_sfc(self):
        if (self.gNumLevels < 3): return 0
        for i in range(0, self.gNumLevels):
            if (QC(self.gSndg[i][self.tind])): return i
        return 0


    def dir2Comp(self):
        for i in range(self.gNumLevels):
            self.gSndg[i][self.uind], self.gSndg[i][self.vind] = \
                vector.vec2comp(self.gSndg[i][self.wdirind],
                self.gSndg[i][self.wspdind])
