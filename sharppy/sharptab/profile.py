'''
Create a Profile Object
'''
from sharppy.sharptab import winds
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
            wdir = kwargs.get('wdir')
            wspd = kwargs.get('wspd')

            self.gSndg = []
            for i in range(len(pres)):
                vals = [pres[i], hght[i], temp[i], dwpt[i], wdir[i], wspd[i]]
                for j in range(0, len(vals)):
                    vals[j] = float(vals[j])
                self.gSndg.append(vals)
            self.gStation = kwargs.get('stn', '????')
            self.gDate = kwargs.get('date', '-----')
            self.gNumLevels = len(pres)

        self.sfc = self.get_sfc()


    def get_sfc(self):
        if (self.gNumLevels < 3): return 0
        for i in range(0, self.gNumLevels):
            if (QC(self.gSndg[i][2])): return i

        return 0


    def dir2Comp(self):
        pass
