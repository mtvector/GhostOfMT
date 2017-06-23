# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 14:29:31 2016

@author: mschmitz
"""

import matplotlib
import matplotlib.pyplot as plt
import math
import PIL
import numpy
import subprocess
from PIL import ImageDraw

class ImageFinder(object):
    """callback for matplotlib to display an annotation when points are
    clicked on.  The point which is closest to the click and within
    xtol and ytol is identified.

    Register this function like this:

    scatter(xdata, ydata)
    af = ImageFinder(xdata, ydata, annotes)
    connect('button_press_event', af)
    """

    def __init__(self, xdata, ydata, annotes,dataDF,path='~/', ax=None, xtol=None, ytol=None):   
        print("theinit")
        self.data = list(zip(xdata, ydata, annotes))
        if xtol is None:
            xtol = ((max(xdata) - min(xdata))/10.0)
        if ytol is None:
            ytol = ((max(ydata) - min(ydata))/10.0)
        self.xtol = xtol
        self.ytol = ytol
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax
        self.drawnAnnotations = {}
        self.dataDF=dataDF
        self.path=path	
        self.image = PIL.Image.open(path+"/Ch2Ch3.tif")
        
    def distance(self, x1, x2, y1, y2):
        """
        return the distance between two points
        """
        return(math.sqrt((x1 - x2)**2 + (y1 - y2)**2))

    def __call__(self, event):

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata
            if (self.ax is None) or (self.ax is event.inaxes):
                annotes = []
                print(event.xdata, event.ydata)
                for x, y, a in self.data:
                    if ((clickX-self.xtol < x < clickX+self.xtol) and
                            (clickY-self.ytol < y < clickY+self.ytol)):
                        annotes.append(
                            (self.distance(x, clickX, y, clickY), x, y, a))
                if annotes:
                    annotes.sort()
                    distance, x, y, annote = annotes[0]
                    self.drawAnnote(event.inaxes, x, y, annote)
                    self.drawImage(annote)

    def drawAnnote(self, ax, x, y, annote):
        """
        Draw the annotation on the plot
        """
        if (x, y) in self.drawnAnnotations:
            markers = self.drawnAnnotations[(x, y)]
            for m in markers:
                m.set_visible(not m.get_visible())
            self.ax.figure.canvas.draw_idle()
        else:
            t = ax.text(x, y, " - %s" % (annote),)
            m = ax.scatter([x], [y], marker='d', c='white', zorder=100)
            self.drawnAnnotations[(x, y)] = (t, m)
            self.ax.figure.canvas.draw_idle()

    def drawSpecificAnnote(self, annote):
        annotesToDraw = [(x, y, a) for x, y, a in self.data if a == annote]
        for x, y, a in annotesToDraw:
            self.drawAnnote(self.ax, x, y, a)        

    def drawImage(self, annote):
        print("pick Image")
        print(self.dataDF.loc[annote,'t'])
        self.image.seek(int(self.dataDF.loc[annote,'t']))
        im = self.image
        print(im)
        #p = subprocess.Popen(["display", self.image])
        print("draw image")
        draw = ImageDraw.Draw(im)
        print("drawn")
        ix=self.dataDF.loc[annote,'x']      
        iy=self.dataDF.loc[annote,'y']
        diam = self.dataDF.loc[annote,'estdiam']/2
        print("draw ellipse")
        draw.ellipse((ix-diam, iy-diam, ix+diam, iy+diam))
        im.show()
        print("subprocessing?")
        #raw_input("Press Enter to continue...")       
        # p.kill()

				
                                                            
#To use this functor you can simply do something like this:

#x = range(10)
#y = range(10)
#annotes = [0,1,2,3,4,5,6,7,8,9]
#annotes = range(len(outDF))
#fig, ax = plt.subplots()
#ax.scatter(x,y)
#af =  ImageFinder(x,y, annotes, dataDF=outDF,path=srcpath, ax=ax)
#fig.canvas.mpl_connect('button_press_event', af)
#plt.show()

