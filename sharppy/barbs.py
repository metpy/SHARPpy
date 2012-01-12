''' Barbs Class '''
import sharppy.sharptab as tab


__all__ = ['Barb']


class Barb(object):

    def __init__(self, canvas, x, y, dir, spd, **kwargs):
        self.gCanvas = canvas                           # Canvas to draw upon
        self.centerx = x                                # Location to draw
        self.centery = y
        self.spd = spd                                  # Wind speed
        self.dir = dir                                  # Wind direction (deg)
        self.size = kwargs.get('size', 4)               # Scale factor to barb
        self.color = kwargs.get('color', "#FFFFFF")     # Color to draw
        self.width = kwargs.get('width', 1)             # Line width
        self.rotate = kwargs.get('rotate', 0)           # Rotation from north
        self.spacing = 9        # Num of elements that will fit on a backbone

        # unit vector
        self.iuu, self.ivv = tab.vector.vec2comp(self.dir + self.rotate, self.size)

        # Start with station circle
        self.gCanvas.create_oval(x-self.size/3, y-self.size/3, x+self.size/3, y+self.size/3,
            outline=self.color, width=self.width)

        # Backbone line
        x1 = self.centerx - (self.iuu)
        y1 = self.centery + (self.ivv)
        x2 = self.centerx - (self.iuu * 10)
        y2 = self.centery + (self.ivv * 10)
        self.gCanvas.create_line(x1, y1, x2, y2, fill=self.color, width=self.width)

        # Variable setup
        sped = spd
        dx = self.iuu * self.size
        dy = self.ivv * self.size
        spcx = -(x2-x1) / self.spacing
        spcy = (y2-y1) / self.spacing
        x1 = x2
        y1 = y2

        # Draw wind flags (increments of 50)
        flag=0
        fullbarb=0
        while (sped > 47):
            flag=1
            x2 = x1 - dy - 2
            y2 = y1 - dx - 2
            x3 = x1 + spcx
            y3 = y1 - spcy
            points = [x1, y1, x2, y2, x3, y3]
            self.gCanvas.create_polygon(points, fill=self.color)
            sped = sped - 50
            x1 = x3
            y1 = y3

        # Single barbs (increments of 10)
        while (sped > 7):
            fullbarb=1
            if flag == 0:
                x1 = x1 - spcx
                y1 = y1 + spcy
                flag = 1
            x2 = x1 - dy
            y2 = y1 - dx
            x1 = x1 + spcx
            y1 = y1 - spcy
            self.gCanvas.create_line(x1, y1, x2, y2, fill=self.color, width=self.width)
            sped = sped - 10

        # Half barb (if needed)
        if (sped > 2):
            if fullbarb == 1:
                x1 += spcx
                y1 -= spcy
            x2 = x1 - ((dy)/2) - (spcx/2)
            y2 = y1 - ((dx)/2) + (spcy/2)
            self.gCanvas.create_line(x1, y1, x2, y2, fill=self.color, width=self.width)
