''' Hodograph Class '''
import math
import sharppy.sharptab as tab
from sharppy.sharptab.constants import *


__all__ = ['Hodo']


class Hodo:

    def __init__(self, canvas, **kwargs):
        self.gCanvas = canvas

        # Canvas Widths
        self.rpad = 0
        self.bpad = 0
        self.wid = kwargs.get('width', 800) - self.rpad
        self.hgt = kwargs.get('height', 800) - self.bpad

        # Total speed difference along Y axis
        self.hodomag = 200

        # Highest level to plot (msl)
        self.hodotop = 15000

        # Where on Canvas to start drawing the Hodog
        self.tlx = self.rpad        # Top-Left X
        self.tly = self.bpad        # Top-Left Y

        # Dimensions of the Hodog Frame
        self.brx = self.wid      # Bottom-Right X
        self.bry = self.hgt      # Bottom-Right Y

        # Center (x,y) position in dir/spd (speed=0)
        self.centerx = self.wid/2 + self.rpad/2
        self.centery = self.hgt/2 + self.bpad/2
        if 'prof' in kwargs:
            prof = kwargs.get('prof')
            self.centerHodo(prof)
        else:
            self.scle = (self.brx - self.tlx) / self.hodomag

        # Fonts
        self.font1 = ("Helvetica", 9)
        self.font2 = ("Helvetica", 11)
        self.font3 = ("Helvetica", 7)

        # Colors
        self.framefg = "#FFFFFF"
        self.framebg = "#000000"
        self.acolor = "#FFFFFF"          # Axes
        self.rcolor = "#660000"          # Radial circles
        self.rcolor = "#444444"          # Radial circles
        self.hcolor = {}
        self.hcolor[1] = "#FF0000"       # 0-1km color
        self.hcolor[2] = "#CC0033"       # 1-2km color
        self.hcolor[3] = "#990066"       # 2-3km color
        self.hcolor[4] = "#660099"       # 3-4km color
        self.hcolor[5] = "#3300CC"       # 4-5km color
        self.hcolor[6] = "#0000FF"       # 5-6km color
        self.hcolor[7] = "#0033CC"       # 6-7km color
        self.hcolor[8] = "#006699"       # 7-8km color
        self.hcolor[9] = "#009966"       # 8-9km color
        self.hcolor[10] = "#00CC33"      # 9-10km color
        self.hcolor[11] = "#00FF00"      # 10-11km color
        self.hcolor[12] = "#33CC00"      # 11-12km color
        self.hcolor[13] = "#669900"      # 12-13km color
        self.hcolor[14] = "#996600"      # 13-14km color
        self.hcolor[15] = "#CC3300"      # 14-15km color
        self.hcolor[16] = "#333300"      # 15+ km color
        self.hcolor[17] = "#333300"      # 16+ km color
        self.hcolor[18] = "#333300"      # 17+ km color
        self.hcolor[19] = "#333300"      # 17+ km color
        self.hcolor[20] = "#333300"      # 17+ km color
        self.stntextcolor = "#FF0000"    # Station ID Text

        # Lines to Plot
        self.ringincrement = 10
        self.rings = range(self.ringincrement, 200, self.ringincrement)


    def drawHodo(self):
        """ Draw the background Hodograph """
        # Plot rings
        for s in self.rings: self.drawRing(s)
        # Plot axes
        self.drawAxes()

        # Colorfill boxes around plotting area to mask lines
        self.gCanvas.create_rectangle((0, 0, self.wid, self.tly),
            fill=self.framebg, outline=self.framebg)
        self.gCanvas.create_rectangle((0, 0, self.tlx, self.hgt),
            fill=self.framebg, outline=self.framebg)
        self.gCanvas.create_rectangle((0, self.bry, self.wid+self.rpad,
            self.hgt+self.bpad), fill=self.framebg, outline=self.framebg)
        self.gCanvas.create_rectangle((self.brx, 0, self.wid+self.rpad,
            self.hgt+self.bpad), fill=self.framebg, outline=self.framebg)

        # Plot frame around SkewT
        self.gCanvas.create_rectangle((self.tlx, self.tly, self.brx, self.bry),
            fill="", outline=self.framefg)


    def drawProfile(self, prof, width=3):
        ''' Draw the Shear Trace '''
        txt = prof.gStation + " - " + prof.gDate
        self.drawShearX(prof, width=width, font=self.font3)


    def drawRing(self, s):
        x2 = self.centerx
        y2 = self.centery
        uu, vv = tab.vector.vec2comp(0,s)
        uu = uu * self.scle
        vv = vv * self.scle
        x1, y1 = self.hodo2Pix(0, s)
        self.gCanvas.create_oval(x2-vv, y2-vv, x2+vv, y2+vv,
            outline=self.rcolor, dash=(4,4))
        self.gCanvas.create_text(x1-3, y1, fill=self.acolor, text=s,
            anchor='e', font=self.font1)
        x1, y1 = self.hodo2Pix(180, s)
        self.gCanvas.create_text(x1-3, y1, fill=self.acolor, text=s,
            anchor='e', font=self.font1)
        x1, y1 = self.hodo2Pix(90, s)
        self.gCanvas.create_text(x1, y1+3, fill=self.acolor, text=s,
            anchor='n', font=self.font1)
        x1, y1 = self.hodo2Pix(270, s)
        self.gCanvas.create_text(x1, y1+3, fill=self.acolor, text=s,
            anchor='n', font=self.font1)

    def drawAxes(self):
        x2 = self.centerx
        y2 = self.centery
        self.gCanvas.create_line(x2, self.tly, x2, self.bry,
            fill=self.acolor)
        self.gCanvas.create_line(self.tlx, y2, self.brx, y2,
            fill=self.acolor)


    def hodo2Pix(self, ang, spd):
        ''' Function to convert a wind dir/spd to x/y coords'''
        midx = self.centerx
        midy = self.centery
        uu, vv = tab.vector.vec2comp(ang,spd)
        xx = midx + (uu * self.scle)
        yy = midy + (vv * self.scle)
        return xx, yy


    def uv2Pix(self, u, v):
        ''' Function to convert a wind u/v to x/y coords'''
        midx = self.centerx
        midy = self.centery
        xx = midx + (u * self.scle)
        yy = midy + (v * self.scle)
        return xx, yy


    def drawShearX(self, prof, width=2, font=None):
        ''' Draw the wind shear vectors '''
        if not font: font=self.font3
        if prof.gNumLevels < 3: return
        dir, spd = tab.vector.comp2vec(prof.gSndg[prof.sfc][prof.uind],
            -prof.gSndg[prof.sfc][prof.vind])
        x2, y2 = self.hodo2Pix(dir, spd)
        lasth = prof.gSndg[prof.sfc][prof.zind]
        color = self.hcolor[1]
        xs = []
        ys = []
        ts = []
        for i in range(prof.gNumLevels):
            h = prof.gSndg[i][prof.zind]
            if lasth > self.hodotop: break
            if QC(prof.gSndg[i][prof.uind]):
                x1 = x2
                y1 = y2
                for z in range(20):
                    hcheck = (z+1) * 1000
                    if ((lasth <= hcheck) and (h >= hcheck)):
                        zz = z+1
                        p = tab.interp.pres(hcheck, prof)
                        # print p
                        u, v = tab.interp.components(p, prof)
                        dir, spd = tab.vector.comp2vec(u, -v)
                        x2, y2 = self.hodo2Pix(dir, spd)
                        self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                            width=width)
                        self.gCanvas.create_rectangle(x2-2,y2-2,x2+2,y2+2,
                            fill=self.acolor)
                        xs.append(x2+2)
                        ys.append(y2+2)
                        ts.append(zz)
                        color = self.hcolor[z+1]

                h = prof.gSndg[i][prof.zind]
                if (h < self.hodotop):
                    x1 = x2
                    y1 = y2
                    dir, spd = tab.vector.comp2vec(prof.gSndg[i][prof.uind],
                        -prof.gSndg[i][prof.vind])
                    x2, y2 = self.hodo2Pix(dir, spd)
                    self.gCanvas.create_line(x1, y1, x2, y2, fill=color,
                        width=width)
                lasth = h

        for x,y,z in zip(xs, ys, ts):
            self.gCanvas.create_text(x, y, fill=self.acolor,
                text=z, anchor='nw', font=font)


    def centerHodo(self, prof):
        minu, maxu = self.findMinMax(prof, prof.uind)
        minv, maxv = self.findMinMax(prof, prof.vind)
        self.hodomag = max(maxu-minu, maxv-minv) + 55
        self.scle = (self.brx - self.tlx) / self.hodomag
        pbot = tab.interp.pres(tab.interp.msl(0, prof), prof)
        ptop = tab.interp.pres(tab.interp.msl(self.hodotop, prof), prof)
        mnu, mnv = tab.winds.mean_wind(pbot, ptop, prof)
        xx, yy = self.uv2Pix(-mnu, mnv)
        self.centerx = xx
        self.centery = yy


    def findMinMax(self, prof, ind):
        maxval = 0
        minval = 0
        for i in range(prof.gNumLevels):
            if not QC(prof.gSndg[i][ind]): continue
            if prof.gSndg[i][prof.zind] > self.hodotop: break
            if maxval < prof.gSndg[i][ind]: maxval = prof.gSndg[i][ind]
            if minval > prof.gSndg[i][ind]: minval = prof.gSndg[i][ind]
        return minval, maxval
