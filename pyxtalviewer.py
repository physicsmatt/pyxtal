#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# GUI module generated by PAGE version 4.17
# In conjunction with Tcl version 8.6
#    Oct 19, 2018 11:13:56 AM CEST  platform: Linux

import sys

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

import pyxtalviewer_support

def vp_start_gui():
    '''Starting point when module is the main routine.'''
    global val, w, root
    root = Tk()
    pyxtalviewer_support.set_Tk_var()
    top = Pyxtal_Viewer (root)
    pyxtalviewer_support.init(root, top)
    root.mainloop()

w = None
def create_Pyxtal_Viewer(root, *args, **kwargs):
    '''Starting point when module is imported by another program.'''
    top = Toplevel(root)
    pyxtalviewer_support.set_Tk_var()
    w = Pyxtal_Viewer(top)
    pyxtalviewer_support.init(top, w, *args, **kwargs)
    return(w)

def destroy_Pyxtal_Viewer():
    global w
    w.destroy()
    print("destroying pixal viewer")
    w = None


class Pyxtal_Viewer:
    def __init__(self, top=None):
        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        _bgcolor = '#d9d9d9'  # X11 color: 'gray85'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9' # X11 color: 'gray85'
        _ana1color = '#d9d9d9' # X11 color: 'gray85' 
        _ana2color = '#d9d9d9' # X11 color: 'gray85' 

        top.geometry("822x867+1085+122")
        top.title("Pyxtal Viewer")
        top.configure(highlightcolor="black")

        #PAGE defined these as global variables, which didn't work for
        #multiple instances of the class.  If GUI is re-generated by PAGE,
        #This will have to be put back in manually
        self.whichImage = StringVar()
        self.invertImage = BooleanVar()
        self.showCircles = BooleanVar()
        self.showTriang = BooleanVar()
        self.showDefects = BooleanVar()
        self.showOrientation = BooleanVar()
        self.showTraject = BooleanVar()
        self.showStats = BooleanVar()


        self.imageframe = LabelFrame(top)
        self.imageframe.place(relx=0.012, rely=0.012, relheight=0.052
                , relwidth=0.341)
        self.imageframe.configure(relief=GROOVE)
        self.imageframe.configure(text='''Image''')
        self.imageframe.configure(width=280)

        self.rawButton = Radiobutton(self.imageframe)
        self.rawButton.place(relx=0.036, rely=0.444, relheight=0.378
                , relwidth=0.161, bordermode='ignore')
        self.rawButton.configure(activebackground="#d9d9d9")
        self.rawButton.configure(command=lambda: pyxtalviewer_support.changeVisibleImage(self))
        self.rawButton.configure(justify=LEFT)
        self.rawButton.configure(text='''Raw''')
        self.rawButton.configure(value="raw")
        self.rawButton.configure(variable=self.whichImage)

        self.filteredButton = Radiobutton(self.imageframe)
        self.filteredButton.place(relx=0.214, rely=0.444, relheight=0.378
                , relwidth=0.286, bordermode='ignore')
        self.filteredButton.configure(activebackground="#d9d9d9")
        self.filteredButton.configure(command=lambda: pyxtalviewer_support.changeVisibleImage(self))
        self.filteredButton.configure(justify=LEFT)
        self.filteredButton.configure(text='''Filtered''')
        self.filteredButton.configure(value="filtered")
        self.filteredButton.configure(variable=self.whichImage)

        self.noneButton = Radiobutton(self.imageframe)
        self.noneButton.place(relx=0.536, rely=0.444, relheight=0.378
                , relwidth=0.186, bordermode='ignore')
        self.noneButton.configure(activebackground="#d9d9d9")
        self.noneButton.configure(command=lambda: pyxtalviewer_support.changeVisibleImage(self))
        self.noneButton.configure(justify=LEFT)
        self.noneButton.configure(text='''None''')
        self.noneButton.configure(value="none")
        self.noneButton.configure(variable=self.whichImage)

        self.invertCheck = Checkbutton(self.imageframe)
        self.invertCheck.place(relx=0.75, rely=0.444, relheight=0.378
                , relwidth=0.225, bordermode='ignore')
        self.invertCheck.configure(activebackground="#d9d9d9")
        self.invertCheck.configure(command=lambda: pyxtalviewer_support.changeVisibleImage(self))
        self.invertCheck.configure(justify=LEFT)
        self.invertCheck.configure(text='''Invert''')
        self.invertCheck.configure(variable=self.invertImage)

        self.annotationsframe = LabelFrame(top)
        self.annotationsframe.place(relx=0.377, rely=0.012, relheight=0.052
                , relwidth=0.523)
        self.annotationsframe.configure(relief=GROOVE)
        self.annotationsframe.configure(text='''Annotations''')
        self.annotationsframe.configure(width=430)

        self.circlesCheck = Checkbutton(self.annotationsframe)
        self.circlesCheck.place(relx=0.023, rely=0.444, relheight=0.378
                , relwidth=0.163, bordermode='ignore')
        self.circlesCheck.configure(activebackground="#d9d9d9")
        self.circlesCheck.configure(anchor=W)
        self.circlesCheck.configure(command=lambda: pyxtalviewer_support.changeVisibleAnnotations(self))
        self.circlesCheck.configure(justify=LEFT)
        self.circlesCheck.configure(text='''circles''')
        self.circlesCheck.configure(variable=self.showCircles)

        self.triangulationCheck = Checkbutton(self.annotationsframe)
        self.triangulationCheck.place(relx=0.209, rely=0.444, relheight=0.378
                , relwidth=0.147, bordermode='ignore')
        self.triangulationCheck.configure(activebackground="#d9d9d9")
        self.triangulationCheck.configure(anchor=W)
        self.triangulationCheck.configure(command=lambda: pyxtalviewer_support.changeVisibleAnnotations(self))
        self.triangulationCheck.configure(justify=LEFT)
        self.triangulationCheck.configure(text='''triang''')
        self.triangulationCheck.configure(variable=self.showTriang)

        self.defectsCheck = Checkbutton(self.annotationsframe)
        self.defectsCheck.place(relx=0.372, rely=0.444, relheight=0.378
                , relwidth=0.163, bordermode='ignore')
        self.defectsCheck.configure(activebackground="#d9d9d9")
        self.defectsCheck.configure(anchor=W)
        self.defectsCheck.configure(command=lambda: pyxtalviewer_support.changeVisibleAnnotations(self))
        self.defectsCheck.configure(justify=LEFT)
        self.defectsCheck.configure(text='''defects''')
        self.defectsCheck.configure(variable=self.showDefects)

        self.orientationCheck = Checkbutton(self.annotationsframe)
        self.orientationCheck.place(relx=0.558, rely=0.444, relheight=0.378
                , relwidth=0.13, bordermode='ignore')
        self.orientationCheck.configure(activebackground="#d9d9d9")
        self.orientationCheck.configure(anchor=W)
        self.orientationCheck.configure(command=lambda: pyxtalviewer_support.changeVisibleAnnotations(self))
        self.orientationCheck.configure(justify=LEFT)
        self.orientationCheck.configure(text='''angle''')
        self.orientationCheck.configure(variable=self.showOrientation)

        self.trajectCheck = Checkbutton(self.annotationsframe)
        self.trajectCheck.place(relx=0.721, rely=0.444, relheight=0.378
                , relwidth=0.27, bordermode='ignore')
        self.trajectCheck.configure(activebackground="#d9d9d9")
        self.trajectCheck.configure(anchor=W)
        self.trajectCheck.configure(command=lambda: pyxtalviewer_support.changeVisibleAnnotations(self))
        self.trajectCheck.configure(justify=LEFT)
        self.trajectCheck.configure(text='''trajectories''')
        self.trajectCheck.configure(variable=self.showTraject)

        self.statsCheck = Checkbutton(top)
        self.statsCheck.place(relx=0.912, rely=0.023, relheight=0.032
                , relwidth=0.068)
        self.statsCheck.configure(activebackground="#d9d9d9")
        self.statsCheck.configure(command=pyxtalviewer_support.showStatsWin)
        self.statsCheck.configure(justify=LEFT)
        self.statsCheck.configure(text='''Show stats''')
        self.statsCheck.configure(variable=self.showStats)
        self.statsCheck.configure(wraplength="40")
        self.statsCheck.configure(command=lambda: pyxtalviewer_support.changeVisibleAnnotations(self))

        self.imgCanvas = Canvas(top)
        self.imgCanvas.place(relx=0.012, rely=0.069, relheight=0.923
                , relwidth=0.973)
        self.imgCanvas.configure(relief=RIDGE)
        self.imgCanvas.configure(selectbackground="#c4c4c4")
        self.imgCanvas.configure(width=841)
        self.imgCanvas.bind('<B1-Motion>',lambda e:pyxtalviewer_support.xxx(e))
        self.imgCanvas.bind('<Button-3>',lambda e:pyxtalviewer_support.xxx(e))
        self.imgCanvas.bind('<ButtonRelease-3>',lambda e:pyxtalviewer_support.xxx(e))
        self.imgCanvas.bind('<MouseWheel>',lambda e:pyxtalviewer_support.xxx(e))






if __name__ == '__main__':
    #vp_start_gui()
    #print("This file is not runnable as main.  Run Pyxtalmain.py instead.")
    import pyxtal
    pyxtal.vp_start_gui()


