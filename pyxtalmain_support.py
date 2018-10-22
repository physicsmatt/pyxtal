#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# Support module generated by PAGE version 4.17
# In conjunction with Tcl version 8.6
#    Oct 15, 2018 11:55:06 AM CEST  platform: Linux
#    Oct 15, 2018 12:26:07 PM CEST  platform: Linux
#    Oct 16, 2018 11:21:03 AM CEST  platform: Linux

import sys
import pyxtalviewer
import pyxtalviewer_support

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

import tkinter.filedialog as fd
import os.path

def set_Tk_var():
    None
    #This function is a place holder, originally created by PAGE
    #It defined a zillion global Tk variables that were used for all
    #the widgets.  That functionality is now within the creation function.


def set_widget_state(to_state, widgets):
    #Input "to_state" is expected to be NORMAL or DISABLED.
    #Input "widget" is a single widget or list of widgets.
    #If widget is a list, or a frame, then 
    if isinstance(widgets, (list, tuple)):
        for item in widgets:
            set_widget_state(to_state, item)
    elif widgets.winfo_class() in ["frame", "Labelframe"]:
        for child in widgets.winfo_children():
            set_widget_state(to_state, child)
    else:
        widgets.configure(state=to_state)
    

def inFileTypeChange():
    global pmw
    inputtype = pmw.inFileType.get()
    if inputtype=='image':
       set_widget_state(NORMAL, [pmw.darkSpheresCheck, pmw.SphereSizeLabel, pmw.sphereEntry])
       set_widget_state(DISABLED, [pmw.framesframe, pmw.PartTypeLabel, pmw.partTypeEntry])
    if inputtype=='particles':
       set_widget_state(NORMAL, [pmw.framesframe])
       set_widget_state(DISABLED, [pmw.darkSpheresCheck, pmw.SphereSizeLabel, pmw.sphereEntry])
       set_widget_state(DISABLED, [pmw.PartTypeLabel, pmw.partTypeEntry])
    if inputtype=='assemblies':
       set_widget_state(NORMAL, [pmw.framesframe, pmw.SphereSizeLabel, pmw.sphereEntry])
       set_widget_state(NORMAL, [pmw.PartTypeLabel, pmw.partTypeEntry])
       set_widget_state(DISABLED, [pmw.darkSpheresCheck])


def validateInteger(p1, thestring, theinteger):
    #This function is used for Entry fields where an integer is required.
    #If user types a non-integer, like "abc" or 4.1, the value is rejected,
    #and the previous value is retained.
    #Some notes: I tried to implement this with the validatecommand apparatus
    #in tkinter, but that seemed more aimed at character-by-character checking.
    #For each entry (like sphereSize, for instance), there's both a string 
    #"sphereSizeStr" which is a tk StringVar, and "sphereSize" which is a list
    #that holds a single integer value.  The list is required so that it can be
    #effectively passed by reference, and thus changed inside this function.
    #STILL TO FIX sometime, maybe: figure out how to call this function when
    #user clicks on a different widget; that is, broaden the meaning of
    #"focusout" when this function is called.

    #print('pyxtalviewer_support.yyyyy', thestring, theinteger)
    #print('p1 = {0}'.format(p1))
    #sys.stdout.flush()
    try:
        theinteger[0]=int(thestring.get())
    #    print("after assignment, fromFrame:", theinteger[0])
    except ValueError:
    #    print("rejected")
        thestring.set(str(theinteger[0]))

def GoButtonCommand():
    sys.stdout.flush()
    pmw.numFiles=len(pmw.filelist)
    for i in range(0,pmw.numFiles):
        pmw.viewers.append(pyxtalviewer.create_Pyxtal_Viewer(
                root, pmw, pmw.filelist[i], i))

def killAllButtonCommand():
    for viewer in pmw.viewers.copy():
        pyxtalviewer_support.destroy_viewer(viewer)

def addButtonCommand():
    #Adds files to the filelist (the variable and the scroll box)...
    filenames = fd.askopenfilenames()
    for filename in filenames:
        pmw.filelist.append(os.path.basename(filename))
        pmw.fileListbox.insert(END, os.path.basename(filename))

    #...And updates the path to wherever the files came from.
    pmw.path = os.path.dirname(filenames[0])
    pmw.pathBox.delete(1.0, END)
    #I don't know why the correct index above is 1.0
    pmw.pathBox.insert(END, pmw.path)


def clearButtonCommand():
    pmw.filelist.clear()
    pmw.fileListbox.delete(0, END)
    #Note: I don't understand why the index for deleting is "0" here, but
    #1.0 for clearing the pathBox.


def saveButtonCommand():
    global pmw
    savefilename = pmw.path + "/pyxtalrc.py"  #(use os.path.join)
    print(savefilename)
    f = open(savefilename, 'w')
    
    f.write("#Global parameter file for Pyxtal\n")
    f.write("#Generated automatically by Pyxtal\n")
    #date and time stamp.   
    global inFileType, darkSpheres, partTypeStr, fromFrame, toFrame, byFrame
    global sphereSize, periodBound
    global outCircles, outTriang, outAll, imageSize, outMpeg, outLog
    global doOrientCorr, doTraject
    global retainWin, lockViews, lockZoom

    f.write('global inFileType, darkSpheres, partTypeStr, fromFrame, toFrame, byFrame\n')
    f.write('global sphereSize, periodBound\n')
    f.write('global outCircles, outTriang, outAll, imageSize, outMpeg, outLog\n')
    f.write('global doOrientCorr, doTraject\n')
    f.write('global retainWin, lockViews, lockZoom\n')


    f.write("byFrame = " + str(byFrame) + "\n")

    f.write('inFileType.set("' + inFileType.get() + '")\n')

    f.write("retainWin.set(" + str(retainWin.get()) + ")\n")
    f.close()
            
    

def loadButtonCommand():
    print("loading parameters.")
    global pmw
    loadfilename = pmw.path + "/pyxtalrc" #(use os.path.join)
    
    #global inFileType, darkSpheres, partTypeStr, fromFrame, toFrame, byFrame
    #global sphereSize, periodBound
    #global outCircles, outTriang, outAll, imageSize, outMpeg, outLog
    #global doOrientCorr, doTraject
    #global retainWin, lockViews, lockZoom
    
    import pyxtalrc
    #see: https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path


def init(top, gui, *args, **kwargs):
    global pmw, top_level, root
    pmw = gui
    top_level = top
    root = top
    pmw.top = top
    top.protocol("WM_DELETE_WINDOW", lambda: destroy_pyxtalmain(pmw))

    pmw.filelist = list()
    pmw.viewers = list()
    initialize_parameters(pmw)

    import os
    os.getcwd()
    pmw.path = os.getcwd()
    pmw.pathBox.insert(END, pmw.path)


def initialize_parameters(pmw):
    #set default views for viewers:
    pmw.whichImage = "raw"
    pmw.invertImage = False
    pmw.showCircles = True
    pmw.showTriang = False
    pmw.showDefects = True
    pmw.showOrientation = True
    pmw.showTraject = False
    pmw.showStats = False

    #initialize all of the Tk variables declared during creation: 
    pmw.inFileType.set("image") 
    pmw.darkSpheres.set(False)
    pmw.partTypeStr.set("")
    pmw.periodBound.set(False)
    pmw.outCircles.set(True)
    pmw.outTriang.set(False)
    pmw.outAll.set(True)
    pmw.outMpeg.set(False)
    pmw.outLog.set(True)
    pmw.doOrientCorr.set(False)
    pmw.doTraject.set(False)
    pmw.retainWin.set(True)
    pmw.lockViews.set(False)
    pmw.lockZoom.set(False)

    #These numeric values function as a way to save previous values
    #if the associated strings are changed to non-integer values.
    pmw.fromFrame = [0]
    pmw.toFrame = [-1]
    pmw.byFrame = [1]
    pmw.sphereSize = [7]
    pmw.imageSize = [-1]

    #Set the string values to the corresponding integers, above
    pmw.fromFrameStr.set (str(pmw.fromFrame[0]))
    pmw.toFrameStr.set (str(pmw.toFrame[0]))
    pmw.byFrameStr.set (str(pmw.byFrame[0]))
    pmw.sphereSizeStr.set (str(pmw.sphereSize[0]))
    pmw.imageSizeStr.set (str(pmw.imageSize[0]))

    inFileTypeChange()


def destroy_pyxtalmain(pmw):
    # Function which closes the window.
#    global top_level
    killAllButtonCommand()    
    pmw.top.destroy()
    top_level = None

if __name__ == '__main__':
    import pyxtalviewer
    import pyxtalmain
    pyxtalmain.vp_start_gui()

