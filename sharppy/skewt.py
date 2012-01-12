import math
import sharppy as sp
import sharppy.sharptab as tab
from sharppy.sharptab.constants import *


__all__ = ['SkewT']


class SkewT:

    def __init__(self, canvas, **kwargs):
        self.gCanvas = canvas

        # Canvas Widths
        self.rpad = 100
        self.bpad = 20
        self.wid = kwargs.get('width', 800) - self.rpad
        self.hgt = kwargs.get('height', 800) - self.bpad

        # Where on Canvas to start drawing the SkewT
        self.tlx = 30        # Top-Left X
        self.tly = 20        # Top-Left Y

        # Dimensions of the SkewT Frame
        self.brx = self.wid      # Bottom-Right X
        self.bry = self.hgt      # Bottom-Right Y

        # Maximum (bottom) & Minimum (top) Pressures
        self.pmin = 100.
        self.pmax = 1075.

        # Drawing Parameters
        self.bltemp = -55       # temperature at the bottom-left of the chart
        self.brtemp = 55        # temperature at the bottom-right of the chart

        # Rotation Angle in Degrees
        self.rot = 100/3.

        # SkewT Fonts
        self.font1 = ("Helvetica", 9)
        self.font2 = ("Helvetica", 11)
        self.font3 = ("Helvetica", 7)

        # Colors
        self.framefg = "#FFFFFF"
        self.framebg = "#000000"
        self.icolor = "#BBBBBB"         # Isobar
        self.tcolor = "#FF0000"         # Temperature Trace
        self.tdcolor = "#00FF00"        # Dewpoint  Trace
        self.twcolor = "#AAAAFF"        # Wetbulb Temperature Trace
        self.tpcolor = "#AAAA00"        # Parcel Trace
        self.tvpcolor = "#FFFF00"       # Virtual Parcel
        self.ithermcolor = "#333333"    # Isotherms
        self.ithermbold = "#CCCCCC"     # Bolded Isotherms
        self.adiabatcolor = "#333333"   # Dry Adiabat
        self.madiabatcolor = "#663333"  # Moist Adiabat
        self.mixratcolor = "#006600"    # Mixing Ratio
        self.stntextcolor = "#FF0000"   # Station ID Text
        self.tcolor = "#FF0000"         # Temperature Trace
        self.tdcolor = "#00FF00"        # Dew Point Trace
        self.twcolor = "#AAAAFF"        # Wetbulb Trace
        self.dgzcolor = "#00FFFF"       # DGZ Trace Color
        self.barbcolor = '#FFFFFF'      # Wind Barb Color

        # Lines to Plot
        self.dp = -25
        self.presrange = range(int(self.pmax), int(self.pmin-1), self.dp)
        self.isobars = [1000, 850, 700, 500, 300, 200, 100]
        self.isotherms = range(-160, 61, 10)
        self.thtas = range(-70, 350, 20)
        self.thtes = range(-160, 61, 10)
        self.mixws = [2] + range(4, 33, 4)
        self.wbot = self.pmax - 5           # Mixing Ratio Bottom Pressure
        self.wtop = 600                     # Mixing Ratio Top Pressure
        self.minTdgz = -18                  # Minimum temperature of DGZ
        self.maxTdgz = -12                  # Maximum temperature of DGZ
        self.tracewidth = 4                 # Tracewidth

        # Update All Keyword Arguments
        self.__dict__.update(kwargs)

        # Horizontal temperature spread (dT across the bottom of the chart)
        self.hspread = self.brtemp - self.bltemp

        # Vertical temperature spread (dT along the left edge of the chart)
        self.vspread = math.tan(math.radians(self.rot)) * self.hspread


    def drawSkewT(self):
        """ Draw the background SkewT """
        btm = int(self.pmax) / 50 * 50
        for p in range(btm, int(self.pmin), -50): self.drawIsobar(p, 1)
        for tw in self.thtes: self.drawMoistAdiabat(tw)
        for t in self.isotherms: self.drawIsotherm(t)
        for thta in self.thtas: self.drawDryAdiabat(thta)
        for w in self.mixws: self.drawMixRatioLine(w, self.font3)

        # Colorfill boxes around plotting area to mask lines
        self.gCanvas.create_rectangle((0, 0, self.tlx, self.bry),
            fill=self.framebg, outline=self.framebg)
        self.gCanvas.create_rectangle((0, self.pres2Pix(self.pmax), self.brx,
            self.bry), fill=self.framebg, outline=self.framebg)
        self.gCanvas.create_rectangle((self.brx, 0, self.wid+self.rpad,
            self.pres2Pix(self.pmax)), fill=self.framebg, outline=self.framebg)
        for isobar in self.isobars: self.drawIsobar(isobar, 0)

        # Plot frame around SkewT
        self.gCanvas.create_rectangle((self.tlx, self.tly, self.brx, self.bry),
            fill="", outline=self.framefg)


    def drawProfile(self, prof, **kwargs):
        ''' Draw the Sounding '''
        # Create the Sounding ID and Date/Time Header
        txt = prof.gStation + "  -  " + prof.gDate
        self.gCanvas.create_text(self.tlx, 2, fill=self.stntextcolor,
            text=txt, anchor='nw', font=self.font2)

        # Create the Model/Obs Header
        self.gCanvas.create_text(self.wid-150, 2, fill=self.stntextcolor,
            text=prof.gModel, anchor='nw', font=self.font2)

        # Add WetBulb to Profile
        prof = self.createWetBulb(prof)

        # Make the Drawings
        twwidth = kwargs.get('twwidth', 1)
        plottxt = kwargs.get('plottxt', True)
        self.__dict__.update(kwargs)
        self.drawTrace(prof, -1, color=self.twcolor, width=twwidth,
            plottxt=plottxt)
        self.drawTrace(prof, prof.tdind, self.tdcolor, width=self.tracewidth,
            plottxt=plottxt)
        self.drawTrace(prof, prof.tind, self.tcolor, width=self.tracewidth,
            plottxt=plottxt)
        self.drawDGZ(prof, self.dgzcolor, width=self.tracewidth)


    def drawBarbs(self, prof, color=None, **kwargs):
        ''' Draw the Wind Barbs '''
        if not color: color = self.barbcolor
        self.plevs = [prof.gSndg[prof.sfc][prof.pind]] + self.presrange
        self.__dict__.update(kwargs)

        if not self.plevs:
            self.plevs = [prof.gSndg[i][prof.sfc]
                for i in range(prof.gNumLevels)]

        for p in self.plevs:
            if p < self.pmin or p > self.pmax or \
                p > prof.gSndg[prof.sfc][prof.pind]: continue
            u, v = tab.interp.components(p, prof)
            y1 = self.pres2Pix(p)
            x1 = self.brx + self.rpad/2
            wdir, wspd = tab.vector.comp2vec(u, v)
            sp.Barb(self.gCanvas, x1, y1, wdir, wspd, color=color,
                **kwargs)


    def drawDGZ(self, prof, color=None, width=3):
        ''' Draw the Dendritic Snow Growth Zone '''
        if not color: color=self.dgzcolor
        if prof.gNumLevels < 3: return
        for i in range(prof.gNumLevels-1):
            if not QC(prof.gSndg[i][prof.tind]): continue
            if prof.gSndg[i][prof.tind] <= self.maxTdgz and \
               prof.gSndg[i][prof.tind] >= self.minTdgz and \
               prof.gSndg[i+1][prof.tind] <= self.maxTdgz and \
               prof.gSndg[i+1][prof.tind] >= self.minTdgz:
                rh = tab.thermo.relh(prof.gSndg[i][prof.pind],
                    prof.gSndg[i][prof.tind], prof.gSndg[i][prof.tdind])
                if rh >= 75:
                    rh2 = tab.thermo.relh(prof.gSndg[i+1][prof.pind],
                        prof.gSndg[i+1][prof.tind],
                        prof.gSndg[i+1][prof.tdind])
                    if rh2 >= 75:
                        x1 = self.temp2Pix(prof.gSndg[i][prof.tind],
                            prof.gSndg[i][prof.pind])
                        y1 = self.pres2Pix(prof.gSndg[i][prof.pind])
                        x2 = self.temp2Pix(prof.gSndg[i+1][prof.tind],
                            prof.gSndg[i+1][prof.pind])
                        y2 = self.pres2Pix(prof.gSndg[i+1][prof.pind])
                        self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                            width=width)


    def drawTrace(self, prof, ind, color, **kwargs):
        ''' Draw the Temperature Trace on the Sounding '''
        font = kwargs.get('font', self.font3)
        width = kwargs.get('width', 4)
        plottxt = kwargs.get('plottxt', True)
        if prof.gNumLevels < 3: return
        x1 = self.temp2Pix(prof.gSndg[prof.sfc][ind],
            prof.gSndg[prof.sfc][prof.pind])
        y1 = self.pres2Pix(prof.gSndg[prof.sfc][prof.pind])
        txt = "%.1f" % tab.thermo.ctof(prof.gSndg[prof.sfc][ind])
        xoff = int((float(len(txt)) / 2.) * font[1]) - 1
        yoff = font[1]
        x2 = 0; y2 = 0
        if plottxt:
            self.gCanvas.create_rectangle((x1-xoff, y1, x1+xoff, y1+2*yoff),
                fill=self.framebg)
            self.gCanvas.create_text(x1, y1+yoff, fill=color, text=txt,
                font=font)

        for i in range(prof.gNumLevels):
            if QC(prof.gSndg[i][ind]):
                x1 = x2
                y1 = y2
                if prof.gSndg[i][0] > self.pmin:
                    x2 = self.temp2Pix(prof.gSndg[i][ind],
                        prof.gSndg[i][prof.pind])
                    y2 = self.pres2Pix(prof.gSndg[i][prof.pind])
                    if x1 <= 0: continue
                    self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                        width=width)
                else:
                    v = tab.interp.interp_from_pres(self.pmin, prof, ind)
                    x2 = self.temp2Pix(v, self.pmin)
                    y2 = self.pres2Pix(self.pmin)
                    self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                        width=width)
                    break


    def drawParcelTrace(self, pcl, width=2, dash=(1,1), color=None):
        ''' Draw the trace of supplied parcel '''
        if not color: self.tpcolor
        p = pcl.pres
        t = pcl.temp
        td = pcl.dwpt
        x1 = self.temp2Pix(t, p)
        y1 = self.pres2Pix(p)
        p2, t2 = tab.thermo.drylift(p, t, td)
        x2 = self.temp2Pix(t2, p2)
        y2 = self.pres2Pix(p2)
        self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
            width=width, dash=dash)

        for i in range(int(p2 + self.dp), int(self.pmin-1), int(self.dp)):
            x1 = x2
            y1 = y2
            t3 = tab.thermo.wetlift(p2, t2, float(i))
            x2 = self.temp2Pix(t3, float(i))
            y2 = self.pres2Pix(float(i))
            if x2 < self.tlx: break
            self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                width=width, dash=dash)


    def drawVirtualParcelTrace(self, pcl, width=2, dash=(1,1), color=None):
        ''' Draw the trace of supplied parcel '''
        if not color: color = self.tvpcolor
        p = pcl.pres
        t = pcl.temp
        td = pcl.dwpt
        x1 = self.temp2Pix(tab.thermo.virtemp(p, t, td), p)
        y1 = self.pres2Pix(p)
        p2, t2 = tab.thermo.drylift(p, t, td)
        x2 = self.temp2Pix(tab.thermo.virtemp(p2, t2, t2), p2)
        y2 = self.pres2Pix(p2)
        self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
            width=width, dash=dash)

        for i in range(int(p2 + self.dp), int(self.pmin-1), int(self.dp)):
            x1 = x2
            y1 = y2
            t3 = tab.thermo.wetlift(p2, t2, float(i))
            x2 = self.temp2Pix(tab.thermo.virtemp(i, t3, t3), float(i))
            y2 = self.pres2Pix(float(i))
            if x2 < self.tlx: break
            self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                width=width, dash=dash)


    def drawDryAdiabat(self, thta):
        ''' Draw dry adiabats on background SkewT '''
        for p in self.presrange:
            t = ((thta + ZEROCNK) / ((1000. / p)**ROCP)) - ZEROCNK
            x = self.temp2Pix(t, p)
            y = self.pres2Pix(p)
            if p == self.pmax:
                x2 = x
                y2 = y
            else:
                x1 = x2
                y1 = y2
                x2 = x
                y2 = y
                self.gCanvas.create_line(x1, y1, x2, y2,
                    fill=self.adiabatcolor, width=1)


    def drawIsotherm(self, t):
        ''' Draw isotherms on background SkewT '''
        x1 = self.temp2Pix(t, self.pmax-5)
        x2 = self.temp2Pix(t, self.pmin)
        if t >= self.bltemp and t <= self.brtemp:
            self.gCanvas.create_text(x1-2, self.bry+2, fill=self.ithermbold,
                text=t, anchor="n", font=self.font1)
        self.gCanvas.create_line(x1, self.bry, x2, self.tly,
            fill=self.ithermcolor, dash=(4, 2), width=1)
        if t == 0 or t==-20:
            self.gCanvas.create_line(x1, self.bry, x2, self.tly,
                fill=self.ithermbold, dash=(4, 2), width=1)


    def drawMoistAdiabat(self, tw, width=1):
        ''' Draw moist adiabats on background SkewT '''
        for p in self.presrange:
            t = tab.thermo.wetlift(1000., tw, p)
            x = self.temp2Pix(t, p)
            y = self.pres2Pix(p)
            if p == self.pmax:
                x2 = x
                y2 = y
            else:
                x1 = x2
                y1 = y2
                x2 = x
                y2 = y
                self.gCanvas.create_line(x1, y1, x2, y2,
                    fill=self.madiabatcolor, width=width)


    def drawIsobar(self, p, pipflag, width=1):
        ''' Draw isobars on background SkewT '''
        y1 = self.pres2Pix(p)
        if pipflag == 0:
            self.gCanvas.create_line(self.tlx, y1, self.brx, y1,
                fill=self.icolor, width=width)
            self.gCanvas.create_text(self.tlx-2, y1,
                fill=self.framefg, text=p, anchor="e", font=self.font1)
        else:
            self.gCanvas.create_line(self.tlx, y1, self.tlx+5, y1,
                fill=self.icolor, width=width)
            self.gCanvas.create_line(self.brx, y1, self.brx-5, y1,
                fill=self.icolor, width=width)


    def drawMixRatioLine(self, w, font, width=1):
        ''' Function to draw mixing ratio lines '''
        t = tab.thermo.temp_at_mixrat(w, self.wbot)
        x1 = self.temp2Pix(t, self.wbot)
        y1 = self.pres2Pix(self.wbot)
        t = tab.thermo.temp_at_mixrat(w, self.wtop)
        x2 = self.temp2Pix(t, self.wtop)
        y2 = self.pres2Pix(self.wtop)
        self.gCanvas.create_line(x1, y1, x2, y2, fill=self.mixratcolor,
            width=width)

        self.gCanvas.create_rectangle((x2-font[1], y2-2*font[1],
            x2+font[1], y2), fill=self.framebg,
            outline=self.framebg)

        self.gCanvas.create_text(x2, y2-font[1], fill=self.mixratcolor,
            text=w, font=font)


    def createWetBulb(self, prof):
        ''' Create the Wetbulb Temperature Array '''
        for i in range(prof.gNumLevels):
            prof.gSndg[i].append(tab.thermo.wetbulb(prof.gSndg[i][prof.pind],
                prof.gSndg[i][prof.tind], prof.gSndg[i][prof.tdind]))
        return prof


    def temp2Pix(self, t, p):
        ''' Function to convert a temperature level to a pixel '''
        scl1 = self.brtemp - (((self.bry - self.pres2Pix(p)) /
                        (self.bry - self.tly)) * self.vspread)
        scl2 = self.brx - (((scl1 - t) / self.hspread) * (self.brx - self.tlx))
        return scl2


    def pres2Pix(self, p):
        ''' Function to convert a pressure level to a pixel level '''
        scl1 = math.log(self.pmax) - math.log(self.pmin)
        scl2 = math.log(self.pmax) - math.log(p)
        return (self.bry - (scl2/scl1) * (self.bry - self.tly))


    def pix2Pres(self, y):
        ''' Function to convert a pixel to a pressure level'''
        scl1 = math.log(self.pmax) - math.log(self.pmin)
        scl2 = self.bry - float(y)
        scl3 = self.bry - self.tly + 1
        return (self.pmax / math.exp((scl2/scl3) * scl1));
