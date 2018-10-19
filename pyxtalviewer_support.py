#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# Support module generated by PAGE version 4.17
# In conjunction with Tcl version 8.6
#    Oct 15, 2018 09:31:29 PM CEST  platform: Linux

import sys
import pyxtalmain_support

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
    "hello"
#    global whichImage
#    whichImage = Stage
#    invertImage = BooleanVar()
#    global showCircles
#    showCircles = BooleanVar()
#    global showTriang
#    showTriang = BooleanVar()
#    global showDefects
#    showDefects = BooleanVar()
#    global showOrientation
#    showOrientation = BooleanVar()
#    global showTraject
#    showTraject = BooleanVar()
#    global showStats
#    showStats = BooleanVar()
    
def set_views_to_globals(gui, pmw):
    
    gui.whichImage.set(pmw.whichImage)
    gui.invertImage.set(pmw.invertImage)
    gui.showCircles.set(pmw.showCircles)
    gui.showTriang.set(pmw.showTriang)
    gui.showDefects.set(pmw.showDefects)
    gui.showOrientation.set(pmw.showOrientation)
    

def changeVisibleAnnotations():
    print('pyxtalviewer_support.changeVisibleAnnotations')
    sys.stdout.flush()

def defects_visiblechangeVisibleAnnotations():
    print('pyxtalviewer_support.defects_visiblechangeVisibleAnnotations')
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

def init(top, gui, *args, **kwargs):
#    global w, top_level, root
#    w = gui
#    top_level = top
#    root = top
    print("now doing init in pv support")
    pmw = args[0]
    filename = args[1]
    viewer_number = args[2]
    set_views_to_globals(gui, pmw)
    top.title("Pyxtal Viewer: " + filename + " [" + str(viewer_number) + "]")
    top.update()
    


def destroy_window():
    # Function which closes the window.
    global top_level
    top_level.destroy()
    top_level = None

if __name__ == '__main__':
    import pyxtalviewer
    pyxtalviewer.vp_start_gui()


