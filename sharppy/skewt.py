from Tkinter import Canvas
import tkFont
import math
import sharppy.sharptab as tab
from sharppy.sharptab.constants import *


__all__ = ['SkewT']

class SkewT:

    def __init__(self, canvas, **kwargs):
        self.gCanvas = canvas

        # Canvas Widths
        self.rpad = 30
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
        self.bltemp = -50       # temperature at the bottom-left of the chart
        self.brtemp = 60        # temperature at the bottom-right of the chart

        # Horizontal temperature spread (dT across the bottom of the chart)
        self.hspread = self.brtemp - self.bltemp

        # Rotation Angle in Degrees
        self.rot = 100/3.

        # Vertical temperature spread (dT along the left edge of the chart)
        self.vspread = math.tan(math.radians(self.rot)) * self.hspread

        # (un-implemented option) Draw a skew-t or other thermodynamic diagram.
        self.type = 1

        # SkewT Fonts
        self.font1 = ("Helvetica", 9)
        self.font2 = ("Helvetica", 11)
        self.font3 = ("Helvetica", 7)

        # Colors
        self.framefg = "#FFFFFF"
        self.framebg = "#000000"
        self.icolor = "#FFFFFF"         # Isobar
        self.tcolor = "#FF0000"         # Temperature Trace
        self.tdcolor = "#00FF00"        # Dewpoint  Trace
        self.twcolor = "#AAAAFF"        # Wetbulb Temperature Trace
        self.tpcolor = "#AAAA00"        # Parcel Trace
        self.tvpcolor = "#FFFF00"       # Virtual Parcel
        self.ithermcolor = "#666666"    # Isotherms
        self.ithermbold = "#FFFFFF"     # Bolded Isotherms
        self.adiabatcolor = "#666666"   # Dry Adiabat
        self.madiabatcolor = "#553333"  # Moist Adiabat
        self.mixratcolor = "#009900"    # Mixing Ratio
        self.stntextcolor = "#FF0000"   # Station ID Text
        self.tcolor = "#FF0000"         # Temperature Trace
        self.tdcolor = "#00FF00"        # Dew Point Trace
        self.wcolor = "#AAAAFF"         # Wetbulb Trace
        self.dgzcolor = "#0000FF"       # DGZ Trace Color

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

        # Plot frame around SkewT
        self.gCanvas.create_rectangle((self.tlx, self.tly, self.brx, self.bry),
            fill="", outline=self.framefg)

        for isobar in self.isobars: self.drawIsobar(isobar, 0)


    def drawProfile(self, profile):
        ''' Draw the Sounding '''
        txt = profile.gStation + " - " + profile.gDate
        self.gCanvas.create_text(self.tlx, 2, fill=self.stntextcolor,
            text=txt, anchor='nw', font=self.font2)
        profile = self.createWetBulb(profile)
        self.drawTrace(profile, -1, self.wcolor, 1, font=self.font3)
        self.drawTrace(profile, 3, self.tdcolor, 3, font=self.font3)
        self.drawTrace(profile, 2, self.tcolor, 3, font=self.font3)
        self.drawDGZ(profile, 2, 3, self.dgzcolor, 2)


    def drawDGZ(self, profile, tind, tdind, color, width):
        ''' Draw the Dendritic Snow Growth Zone '''
        if profile.gNumLevels < 3: return
        for i in range(profile.gNumLevels-1):
            if not QC(profile.gSndg[i][tind]): continue
            if profile.gSndg[i][tind] <= self.maxTdgz and \
               profile.gSndg[i][tind] >= self.minTdgz and \
               profile.gSndg[i+1][tind] <= self.maxTdgz and \
               profile.gSndg[i+1][tind] >= self.minTdgz:
                rh = tab.thermo.relh(profile.gSndg[i][0],
                    profile.gSndg[i][tind], profile.gSndg[i][tdind])
                if rh >= 75:
                    rh2 = tab.thermo.relh(profile.gSndg[i+1][0],
                        profile.gSndg[i+1][tind], profile.gSndg[i+1][tdind])
                    if rh2 >= 75:
                        x1 = self.temp2Pix(profile.gSndg[i][tind],
                            profile.gSndg[i][0])
                        y1 = self.pres2Pix(profile.gSndg[i][0])
                        x2 = self.temp2Pix(profile.gSndg[i+1][tind],
                            profile.gSndg[i+1][0])
                        y2 = self.pres2Pix(profile.gSndg[i+1][0])
                        self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                            width=width)


    def drawTrace(self, profile, ind, color, width, font):
        ''' Draw the Temperature Trace on the Sounding '''
        if profile.gNumLevels < 3: return
        x1 = self.temp2Pix(profile.gSndg[profile.sfc][ind],
            profile.gSndg[profile.sfc][0])
        y1 = self.pres2Pix(profile.gSndg[profile.sfc][0])
        txt = "%.1f" % tab.thermo.ctof(profile.gSndg[profile.sfc][ind])
        xoff = int((float(len(txt)) / 2.) * font[1]) - 1
        yoff = font[1]
        x2 = 0
        y2 = 0
        self.gCanvas.create_rectangle((x1-xoff, y1, x1+xoff, y1+2*yoff),
            fill=self.framebg)
        self.gCanvas.create_text(x1, y1+yoff, fill=color, text=txt,
            font=font)

        for i in range(profile.gNumLevels):
            if QC(profile.gSndg[i][ind]):
                x1 = x2
                y1 = y2
                if profile.gSndg[i][0] > self.pmin:
                    x2 = self.temp2Pix(profile.gSndg[i][ind], profile.gSndg[i][0])
                    y2 = self.pres2Pix(profile.gSndg[i][0])
                    if x1 <= 0: continue
                    self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                        width=width)
                else:
                    v = tab.interp.interp_from_pres(self.pmin, profile, ind)
                    x2 = self.temp2Pix(v, self.pmin)
                    y2 = self.pres2Pix(self.pmin)
                    self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                        width=width)
                    break


    def drawParcelTrace(self, pcl):
        ''' Draw the trace of supplied parcel '''
        p = pcl.pres
        t = pcl.temp
        td = pcl.dwpt
        x1 = self.temp2Pix(t, p)
        y1 = self.pres2Pix(p)
        p2, t2 = tab.thermo.drylift(p, t, td)
        x2 = self.temp2Pix(t2, p2)
        y2 = self.pres2Pix(p2)
        self.gCanvas.create_line(x1, y1, x2, y2, fill=self.tpcolor,
            width=2, dash=(1,1))

        for i in range(int(p2 + self.dp), int(self.pmin-1), int(self.dp)):
            x1 = x2
            y1 = y2
            t3 = tab.thermo.wetlift(p2, t2, float(i))
            x2 = self.temp2Pix(t3, float(i))
            y2 = self.pres2Pix(float(i))
            if x2 < self.tlx: break
            self.gCanvas.create_line(x1, y1, x2, y2, fill=self.tpcolor,
                width=2, dash=(1,1))


    def drawVirtualParcelTrace(self, pcl, **kwargs):
        ''' Draw the trace of supplied parcel '''
        color = kwargs.get('color', self.tvpcolor)
        p = pcl.pres
        t = pcl.temp
        td = pcl.dwpt
        x1 = self.temp2Pix(tab.thermo.virtemp(p, t, td), p)
        y1 = self.pres2Pix(p)
        p2, t2 = tab.thermo.drylift(p, t, td)
        x2 = self.temp2Pix(tab.thermo.virtemp(p2, t2, t2), p2)
        y2 = self.pres2Pix(p2)
        self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
            width=2, dash=(1,1))

        for i in range(int(p2 + self.dp), int(self.pmin-1), int(self.dp)):
            x1 = x2
            y1 = y2
            t3 = tab.thermo.wetlift(p2, t2, float(i))
            x2 = self.temp2Pix(tab.thermo.virtemp(i, t3, t3), float(i))
            y2 = self.pres2Pix(float(i))
            if x2 < self.tlx: break
            self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                width=2, dash=(1,1))


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


    def drawMoistAdiabat(self, tw):
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
                    fill=self.madiabatcolor, width=1)


    def drawIsobar(self, p, pipflag):
        ''' Draw isobars on background SkewT '''
        y1 = self.pres2Pix(p)
        if pipflag == 0:
            self.gCanvas.create_line(self.tlx, y1, self.brx, y1,
                fill=self.icolor, width=1)
            self.gCanvas.create_text(self.tlx-2, y1,
                fill=self.icolor, text=p, anchor="e", font=self.font1)
        else:
            self.gCanvas.create_line(self.tlx, y1, self.tlx+5, y1,
                fill=self.icolor, width=1)
            self.gCanvas.create_line(self.brx, y1, self.brx-5, y1,
                fill=self.icolor, width=1)


    def drawMixRatioLine(self, w, font):
        ''' Function to draw mixing ratio lines '''
        t = tab.thermo.temp_at_mixrat(w, self.wbot)
        x1 = self.temp2Pix(t, self.wbot)
        y1 = self.pres2Pix(self.wbot)
        t = tab.thermo.temp_at_mixrat(w, self.wtop)
        x2 = self.temp2Pix(t, self.wtop)
        y2 = self.pres2Pix(self.wtop)
        self.gCanvas.create_line(x1, y1, x2, y2, fill=self.mixratcolor,
            width=1)

        self.gCanvas.create_rectangle((x2-font[1], y2-2*font[1],
            x2+font[1], y2), fill=self.framebg,
            outline=self.framebg)

        self.gCanvas.create_text(x2, y2-font[1], fill=self.mixratcolor,
            text=w, font=font)


    def createWetBulb(self, profile):
        ''' Create the Wetbulb Temperature Array '''
        for i in range(profile.gNumLevels):
            profile.gSndg[i].append(tab.thermo.wetbulb(profile.gSndg[i][0],
                profile.gSndg[i][2], profile.gSndg[i][3]))
        return profile


    def temp2Pix(self, t, p):
        ''' Function to convert a temperature level to a pixel '''
        if self.type == 1:
            scl1 = self.brtemp - (((self.bry - self.pres2Pix(p)) /
                (self.bry - self.tly)) * self.vspread)
        else:
            scl1 = self.brtemp
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
