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
#import gc
import tkinter.filedialog as fd
#import tkfilebrowser as fd  Maybe this is better?
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
        set_widget_state('normal', [pmw.darkSpheresCheck, pmw.sphereSizeLabel, pmw.sphereEntry])
        set_widget_state('disabled', [pmw.framesFrame, pmw.partTypeLabel, pmw.partTypeEntry])
    if inputtype=='particles':
        set_widget_state('normal', [pmw.framesFrame])
        set_widget_state('disabled', [pmw.darkSpheresCheck, pmw.sphereSizeLabel, pmw.sphereEntry])
        set_widget_state('disabled', [pmw.partTypeLabel, pmw.partTypeEntry])
    if inputtype=='assemblies':
        set_widget_state('normal', [pmw.framesFrame, pmw.sphereSizeLabel, pmw.sphereEntry])
        set_widget_state('normal', [pmw.partTypeLabel, pmw.partTypeEntry])
        set_widget_state('disabled', [pmw.darkSpheresCheck])


def batchmodeChange():
    global pmw
    if pmw.batchmode.get():
        pmw.retainWin.set(False)
        set_widget_state('disabled', [pmw.retainCheck])
    else:
        set_widget_state('normal', [pmw.retainCheck])


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

def create_logfiles(pmw, filename):
    import os
    splitpoint = filename.rfind('.')
    base = os.path.join(pmw.path, filename[0:splitpoint])
    if pmw.doOrientHist.get():
        pmw.orientHistfile = open(base + "_orientHist.txt", 'w')
        pmw.orientHistfile.write("# First column is time step\n")
        pmw.orientHistfile.write("# Next columns are area bins, 0 to 60 degrees, by 2 degrees\n")
    if pmw.doMeaningLife.get():
        pmw.meaningLifefile = open(base + "_meaningLife.txt", 'w')
        pmw.meaningLifefile.write("42\n")
        pmw.meaningLifefile.flush()
    if pmw.doZProfile.get():
        pmw.zProfilefile = open(base + "_zProfile.txt", 'w')

        
def GoButtonCommand():
    #This function actually creates the viewer windows for each file (if image)
    #or each image and snapshot, if input is a gsd file.
    #The info passed to create the viewer is the full filename (with path),
    #the index number of the viewer (0 to whatever) and the frame number.
    import gsd.hoomd
    numFiles = len(pmw.filelist)
    vieweridx = len(pmw.viewers)
    if pmw.inFileType.get() == "image":
        #For image input, only one log file created, for first image
        create_logfiles(pmw, pmw.filelist[0])
        for fileidx in range(0,numFiles):
            pmw.viewers.append(pyxtalviewer.create_Pyxtal_Viewer(root, pmw, 
                    pmw.filelist[fileidx], vieweridx, 0))
            vieweridx += 1
    else: #must be some kind of gsd file
        for fileidx in range(0,numFiles):
            pmw.fileMessage.configure(text="Files: "
                                      +str(fileidx+1) +"/" + str(numFiles))
            filename = pmw.filelist[fileidx]
            #For gsd input, each gsd file gets a separate logfile:
            create_logfiles(pmw, filename) 
            start = pmw.fromFrame[0]
            end = pmw.toFrame[0]
            if end == -1:
                full_filename = os.path.join(pmw.path, filename)
                s = gsd.hoomd.open(name=full_filename, mode='rb')
                end = len(s)
            by = pmw.byFrame[0]
            numframes = int((end - start ) / by)
            framecount = 0
            for frameidx in range(start, end, by):
                framecount += 1
                pmw.frameMessage.configure(text="Frames: "
                                      +str(framecount) +"/" + str(numframes))
                newv = pyxtalviewer.create_Pyxtal_Viewer(root, pmw, 
                        filename, vieweridx, frameidx)
                pmw.viewers.append(newv)
                if not pmw.retainWin.get():
                    pyxtalviewer_support.destroy_viewer(newv)
                vieweridx += 1
                

def killAllButtonCommand():
    for viewer in pmw.viewers.copy():
        pyxtalviewer_support.destroy_viewer(viewer)

def addButtonCommand():
    #Adds files to the filelist (the variable and the scroll box)...
    filenames = fd.askopenfilenames(initialdir=pmw.path)
    for filename in filenames:
        pmw.filelist.append(os.path.basename(filename))
        pmw.fileListbox.insert("end", os.path.basename(filename))

    #...And updates the path to wherever the files came from.
    if len(filenames) > 0:
        pmw.path = os.path.dirname(filenames[0])
        pmw.pathBox.delete(1.0, "end")
        #I don't know why the correct index above is 1.0
        pmw.pathBox.insert("end", pmw.path)


def clearButtonCommand():
    pmw.filelist.clear()
    pmw.fileListbox.delete(0, "end")
    #Note: I don't understand why the index for deleting is "0" here, but
    #1.0 for clearing the pathBox.


def saveButtonCommand():
    #This command will eventually save the various parameters to a file,
    #probably .pyxtalrc or something like that.  Not implemented yet;
    #does not work.
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
    #This command will eventually save the various parameters to a file,
    #probably .pyxtalrc or something like that.  Not implemented yet;
    #does not work.
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

def pmw_key_event(event, pmw):
    #print("pmw key event: ",event)
    if event.keysym == "Return":
        pmw.top.focus_force()


def hide_stats(pmw):
    pmw.stats.top.withdraw()
    pmw.stats.viewer.showStats.set(False)

def init(top, gui, *args, **kwargs):
    global pmw, top_level, root
    pmw = gui
    top_level = top
    root = top
    pmw.top = top #this "top" is really the tk root.
    top.protocol("WM_DELETE_WINDOW", lambda: destroy_pyxtalmain(pmw))

    import os
    os.getcwd()
    pmw.path = os.getcwd()
#    pmw.path = "/home/mtrawick/Documents/simulations/2d_diblock/half-loop"
#    pmw.path = "/home/mtrawick/Documents/simulations/half_loop_particles/take2"
#    pmw.path = "/home/mtrawick/Documents/simulations/2d_diblock/half-loop/take2"
    pmw.pathBox.insert("end", pmw.path)

    pmw.filelist = list()
    pmw.viewers = list()

    import pyxtalstats
    pmw.stats = pyxtalstats.create_pyxtal_stats_win(pmw.top)
    pmw.statsText = pmw.stats.Scrolledtext1
    for i in range(0,15):
        pmw.statsText.insert("end", "Number of sandwiches: " + str(i) + "\n")
    pmw.stats.top.withdraw()
    pmw.stats.top.protocol("WM_DELETE_WINDOW", lambda: hide_stats(pmw))

    initialize_parameters(pmw)

    #as this is a work in progress, I'm disabling controls that are
    #not implemented yet:
    set_widget_state('disabled', pmw.analysisFrame)
    set_widget_state('normal', pmw.orientHistCheck)
    set_widget_state('normal', pmw.meaningLifeCheck)
    set_widget_state('normal', pmw.zProfileCheck)
    set_widget_state('disabled', [pmw.outMpegCheck, 
                                  pmw.imageSizeLabel, pmw.imageSizeEntry])
    set_widget_state('disabled', [pmw.saveDefButton, pmw.loadDefButton])
    if pmw.batchmode.get(): set_widget_state('disabled', pmw.retainCheck)

    #For debugging and demonstration purposes, it's handy to have a default 
    #filename already loaded up.
#    filename = "double.tif"
#    filename = "hex1short.gsd"
    filename = "small_2d_test.gsd"
#    filename = "test_diblock1.gsd"
#    filename = "HLsph_uh6_uw10_se41_kT1.0_ts5000000.gsd"
    pmw.filelist.append(filename)
    pmw.fileListbox.insert("end", os.path.basename(filename))

    pmw.top.bind("<Key>", lambda e:pmw_key_event(e, pmw))


def initialize_parameters(pmw):
    #set default views for viewers:
    pmw.whichImage = "raw"
    pmw.invertImage = False
    pmw.showCircles = False
    pmw.showTriang = True
    pmw.showDefects = True
    pmw.showOrientation = True
    pmw.showTraject = False
    pmw.showStats = False
    pmw.global_zoom = 1.00
    pmw.global_corners = None
    pmw.global_corners_set = False

    #initialize all of the Tk variables declared during creation: 
#    pmw.inFileType.set("particles") 
    pmw.inFileType.set("assemblies") 
    pmw.darkSpheres.set(False)
    pmw.partTypeStr.set("B")
    pmw.periodBound.set(True)
    pmw.outCircles.set(False)
    pmw.outTriang.set(False)
    pmw.outAll.set(False)
    pmw.outMpeg.set(False)

    pmw.doZProfile.set(True)
    pmw.doSphereStats.set(False)
    pmw.doTraject.set(False)
    pmw.doDefectStats.set(False)
    pmw.doOrientHist.set(False)
    pmw.doOrientCorr.set(False)
    pmw.doMeaningLife.set(False)

    pmw.batchmode.set(False)
    pmw.retainWin.set(True)
    pmw.lockViews.set(True)
    pmw.lockZoom.set(True)

    #These numeric values function as a way to save previous values
    #if the associated strings are changed to non-integer values.
    pmw.fromFrame = [50]
    pmw.toFrame = [52]
    pmw.byFrame = [50]
    pmw.sphereSize = [7]
    pmw.imageSize = [-1]

    #Set the string values to the corresponding integers, above
    pmw.fromFrameStr.set (str(pmw.fromFrame[0]))
    pmw.toFrameStr.set (str(pmw.toFrame[0]))
    pmw.byFrameStr.set (str(pmw.byFrame[0]))
    pmw.sphereSizeStr.set (str(pmw.sphereSize[0]))
    pmw.imageSizeStr.set (str(pmw.imageSize[0]))

    #Run this procedure so that the active/disabled options reflect the defaults
    inFileTypeChange()


def destroy_pyxtalmain(pmw):
    # Function which closes the window.
#    global top_level
    #Not sure if I need to do this:
    try:
        pmw.orientHistfile.close()
    except:
        pass
    killAllButtonCommand()    
    pmw.top.destroy()
    top_level = None

if __name__ == '__main__':
    import pyxtalviewer
    import pyxtal
    pyxtal.vp_start_gui()

