#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# Support module generated by PAGE version 4.17
# In conjunction with Tcl version 8.6
#    Oct 15, 2018 09:31:29 PM CEST  platform: Linux

import sys
import pyxtalmain_support
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_agg import FigureCanvasAgg
#https://matplotlib.org/gallery/user_interfaces/embedding_in_tk_canvas_sgskip.html
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg, NavigationToolbar2Tk)

try:
    from Tkinter import *
except ImportError:
    from tkinter import *

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True

def set_Tk_var():
    None
    #This function is a place holder, originally created by PAGE
    #It defined a zillion global Tk variables that were used for all
    #the widgets.  That functionality is now within the creation function.

    
def set_views_to_globals(viewer):
    viewer.whichImage.set(viewer.pmw.whichImage)
    viewer.invertImage.set(viewer.pmw.invertImage)
    viewer.showCircles.set(viewer.pmw.showCircles)
    viewer.showTriang.set(viewer.pmw.showTriang)
    viewer.showDefects.set(viewer.pmw.showDefects)
    viewer.showOrientation.set(viewer.pmw.showOrientation)
    

def changeVisibleAnnotations():
    print('pyxtalviewer_support.changeVisibleAnnotations')
    sys.stdout.flush()

def invertImageChange():
    print('pyxtalviewer_support.invertImageChange')
    sys.stdout.flush()

def showImageChange():
    print('pyxtalviewer_support.showImageChange')
    sys.stdout.flush()

def showStatsWin():
    print('pyxtalviewer_support.showStatsWin')
    sys.stdout.flush()

def xxx(p1):
    print('pyxtalviewer_support.xxx')
    print('p1 = {0}'.format(p1))
    sys.stdout.flush()

def load_images_and_locations(viewer):
    #Based on the input file type, this function reads the file.
    #If File is an image, it adds the location data.
    #If File is location data, it adds a fake "image" of spheres.
    #If File is assemblies, it calcultes both an image and location data.
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import trackpy as tp
    print("loading images")
    if viewer.pmw.inFileType.get() == "image":
        #use code from colloid group.
        print("finding particles")
        viewer.image = plt.imread(viewer.filename)
        
        #This gives dataframe with 8 columns. First two are y, x 
        full_locations = tp.locate(viewer.image, viewer.pmw.sphereSize[0])
        viewer.locations = np.array(full_locations)[:,0:2]
        print("found 'em")
        
    elif viewer.pmw.inFileType.get() == "particles":
        #read gsd file.
        None

def setup_images(viewer):

    viewer.fig, ax = plt.subplots()
    x = range(300)
    ax.plot(x, x, '--', linewidth=5, color='firebrick',zorder=1)

    viewer.rawimg = ax.imshow(viewer.image, extent=[0, 400, 0, 300],zorder=0)

    x2 = range(20,200,10)
    #self.dogimg = ax.imshow(thedog, extent=[0, 400, 0, 300],zorder=0.5)#, alpha=0.5) 
    blueplot = ax.scatter(x2, x2, color='blue',zorder=2)
    blueplot.set_visible(1)
    ax.axis('off')
    viewer.fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)

    #Here's where we actually put the plot on the tk canvas:
    viewer.imgCanvas = FigureCanvasTkAgg(viewer.fig, master=viewer.top)
    viewer.imgCanvas.draw()
#    viewer.imgCanvas.get_tk_widget().pack()
    viewer.imgCanvas.get_tk_widget().place(relx=0.012, rely=0.069, relheight=0.923
                , relwidth=0.973)
#    viewer.img_frame.pack()

    

def init(top, viewer, *args, **kwargs):
    viewer.top = top
    viewer.pmw = args[0]
    viewer.filename = args[1]
    viewer.idx = args[2]
    set_views_to_globals(viewer)
    viewer.top.title("Pyxtal Viewer: "
                     + viewer.filename 
                     + " [" + str(viewer.idx) + "]")
    viewer.top.protocol("WM_DELETE_WINDOW", lambda: destroy_viewer(viewer))

    viewer.top.update()    
    load_images_and_locations(viewer)
    setup_images(viewer)    
    viewer.top.update()    

def destroy_viewer(viewer):
    # Function which closes the window.
    viewer.pmw.viewers.remove(viewer) #remove from main list of viewers
    plt.close(viewer.fig) #keeps the plot from reappearing in the console.
    top = viewer.top
    top.destroy()

if __name__ == '__main__':
    #import pyxtalviewer
    #pyxtalviewer.vp_start_gui()
    print("This file is not runnable as main.  Run Pyxtalmain.py instead.")


