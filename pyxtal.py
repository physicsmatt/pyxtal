#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# GUI module generated by PAGE version 4.18
#  in conjunction with Tcl version 8.6
#    Jan 16, 2019 02:52:59 PM CET  platform: Linux

import sys

try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True

import pyxtal_support

def vp_start_gui():
    '''Starting point when module is the main routine.'''
    global val, w, root
    root = tk.Tk()
    pyxtal_support.set_Tk_var()
    top = Pyxtal_Main_Controls (root)
    pyxtal_support.init(root, top)
    root.mainloop()

w = None
def create_Pyxtal_Main_Controls(root, *args, **kwargs):
    '''Starting point when module is imported by another program.'''
    global w, w_win, rt
    rt = root
    w = tk.Toplevel (root)
    pyxtal_support.set_Tk_var()
    top = Pyxtal_Main_Controls (w)
    pyxtal_support.init(w, top, *args, **kwargs)
    return (w, top)

def destroy_Pyxtal_Main_Controls():
    global w
    w.destroy()
    w = None

class Pyxtal_Main_Controls:
    def __init__(self, top=None):
        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        _bgcolor = '#d9d9d9'  # X11 color: 'gray85'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9' # X11 color: 'gray85'
        _ana1color = '#d9d9d9' # X11 color: 'gray85' 
        _ana2color = '#d9d9d9' # X11 color: 'gray85' 
        font10 = "-family gothic -size 9 -weight normal -slant roman "  \
            "-underline 0 -overstrike 0"
        font9 = "-family fixed -size 10 -weight bold -slant roman "  \
            "-underline 0 -overstrike 0"
        self.style = ttk.Style()
        if sys.platform == "win32":
            self.style.theme_use('winnative')
        self.style.configure('.',background=_bgcolor)
        self.style.configure('.',foreground=_fgcolor)
        self.style.map('.',background=
            [('selected', _compcolor), ('active',_ana2color)])

        top.geometry("601x728+472+131")
        top.title("Pyxtal Main Controls")
        top.configure(relief="sunken")
        top.configure(highlightcolor="black")

        #PAGE defined these as global variables, which didn't work for
        #multiple instances of the class.  If GUI is re-generated by PAGE,
        #This will have to be put back in manually
        self.inFileType = tk.StringVar()
        self.darkSpheres = tk.BooleanVar()
        self.partTypeStr = tk.StringVar()
        self.periodBound = tk.BooleanVar()
        self.outCircles = tk.BooleanVar()
        self.outTriang = tk.BooleanVar()
        self.outAll = tk.BooleanVar()
        self.outMpeg = tk.BooleanVar()
        self.doZProfile = tk.BooleanVar()
        self.doSphereStats = tk.BooleanVar()
        self.doTraject = tk.BooleanVar()
        self.doDefectStats = tk.BooleanVar()
        self.doOrientHist = tk.BooleanVar()
        self.doOrientCorr = tk.BooleanVar()
        self.doMeaningLife = tk.BooleanVar()
        self.retainWin = tk.BooleanVar()
        self.batchmode = tk.BooleanVar()
        self.lockViews = tk.BooleanVar()
        self.lockZoom = tk.BooleanVar()
        self.fromFrameStr = tk.StringVar()
        self.toFrameStr = tk.StringVar()
        self.byFrameStr = tk.StringVar()
        self.sphereSizeStr = tk.StringVar()
        self.imageSizeStr = tk.StringVar()


        self.inputFileFrame = tk.LabelFrame(top)
        self.inputFileFrame.place(relx=0.017, rely=0.11, relheight=0.46
                , relwidth=0.483)
        self.inputFileFrame.configure(relief='groove')
        self.inputFileFrame.configure(text='''Input Files''')
        self.inputFileFrame.configure(width=290)

        self.fileListbox = ScrolledListBox(self.inputFileFrame)
        self.fileListbox.place(relx=0.034, rely=0.478, relheight=0.507
                , relwidth=0.952, bordermode='ignore')
        self.fileListbox.configure(background="white")
        self.fileListbox.configure(font="TkFixedFont")
        self.fileListbox.configure(highlightcolor="#d9d9d9")
        self.fileListbox.configure(selectbackground="#c4c4c4")
        self.fileListbox.configure(width=10)

        self.addButton = tk.Button(self.inputFileFrame)
        self.addButton.place(relx=0.034, rely=0.09, height=33, width=98
                , bordermode='ignore')
        self.addButton.configure(activebackground="#d9d9d9")
        self.addButton.configure(command=pyxtal_support.addButtonCommand)
        self.addButton.configure(text='''Add Files''')

        self.clearButton = tk.Button(self.inputFileFrame)
        self.clearButton.place(relx=0.586, rely=0.09, height=33, width=108
                , bordermode='ignore')
        self.clearButton.configure(activebackground="#d9d9d9")
        self.clearButton.configure(command=pyxtal_support.clearButtonCommand)
        self.clearButton.configure(text='''Clear Files''')

        self.pathBox = tk.Text(self.inputFileFrame)
        self.pathBox.place(relx=0.034, rely=0.299, relheight=0.09, relwidth=0.917
                , bordermode='ignore')
        self.pathBox.configure(background="white")
        self.pathBox.configure(font="TkTextFont")
        self.pathBox.configure(selectbackground="#c4c4c4")
        self.pathBox.configure(width=266)
        self.pathBox.configure(wrap='word')

        self.filesLabel = tk.Label(self.inputFileFrame)
        self.filesLabel.place(relx=0.034, rely=0.418, height=17, width=46
                , bordermode='ignore')
        self.filesLabel.configure(activebackground="#f9f9f9")
        self.filesLabel.configure(font=font9)
        self.filesLabel.configure(text='''Files:''')

        self.pathLabel = tk.Label(self.inputFileFrame)
        self.pathLabel.place(relx=0.034, rely=0.239, height=17, width=46
                , bordermode='ignore')
        self.pathLabel.configure(activebackground="#f9f9f9")
        self.pathLabel.configure(font=font9)
        self.pathLabel.configure(text='''Path:''')

        self.saveDefButton = tk.Button(top)
        self.saveDefButton.place(relx=0.283, rely=0.041, height=33, width=108)
        self.saveDefButton.configure(activebackground="#d9d9d9")
        self.saveDefButton.configure(command=pyxtal_support.saveButtonCommand)
        self.saveDefButton.configure(text='''Save Defaults''')

        self.inputOptionsFrame = tk.LabelFrame(top)
        self.inputOptionsFrame.place(relx=0.017, rely=0.591, relheight=0.391
                , relwidth=0.483)
        self.inputOptionsFrame.configure(relief='groove')
        self.inputOptionsFrame.configure(text='''Input File Options''')
        self.inputOptionsFrame.configure(width=290)

        self.imageRadio = tk.Radiobutton(self.inputOptionsFrame)
        self.imageRadio.place(relx=0.034, rely=0.105, relheight=0.06
                , relwidth=0.497, bordermode='ignore')
        self.imageRadio.configure(activebackground="#d9d9d9")
        self.imageRadio.configure(anchor='w')
        self.imageRadio.configure(command=pyxtal_support.inFileTypeChange)
        self.imageRadio.configure(justify='left')
        self.imageRadio.configure(text='''Image file''')
        self.imageRadio.configure(value="image")
        self.imageRadio.configure(variable=self.inFileType)

        self.gsdAssemRadio = tk.Radiobutton(self.inputOptionsFrame)
        self.gsdAssemRadio.place(relx=0.034, rely=0.316, relheight=0.06
                , relwidth=0.548, bordermode='ignore')
        self.gsdAssemRadio.configure(activebackground="#d9d9d9")
        self.gsdAssemRadio.configure(anchor='w')
        self.gsdAssemRadio.configure(command=pyxtal_support.inFileTypeChange)
        self.gsdAssemRadio.configure(justify='left')
        self.gsdAssemRadio.configure(text='''gsd assemblies''')
        self.gsdAssemRadio.configure(value="assemblies")
        self.gsdAssemRadio.configure(variable=self.inFileType)

        self.gsdPartRadio = tk.Radiobutton(self.inputOptionsFrame)
        self.gsdPartRadio.place(relx=0.034, rely=0.211, relheight=0.06
                , relwidth=0.455, bordermode='ignore')
        self.gsdPartRadio.configure(activebackground="#d9d9d9")
        self.gsdPartRadio.configure(anchor='w')
        self.gsdPartRadio.configure(command=pyxtal_support.inFileTypeChange)
        self.gsdPartRadio.configure(justify='left')
        self.gsdPartRadio.configure(text='''gsd particles''')
        self.gsdPartRadio.configure(value="particles")
        self.gsdPartRadio.configure(variable=self.inFileType)

        self.sphereEntry = tk.Entry(self.inputOptionsFrame)
        self.sphereEntry.place(relx=0.414, rely=0.754, height=27, relwidth=0.228
                , bordermode='ignore')
        self.sphereEntry.configure(background="white")
        self.sphereEntry.configure(font="TkFixedFont")
        self.sphereEntry.configure(selectbackground="#c4c4c4")
        self.sphereEntry.configure(textvariable=self.sphereSizeStr)
        self.sphereEntry.bind('<FocusOut>',lambda e:pyxtal_support.validateInteger(e,
                                    self.sphereSizeStr,
                                    self.sphereSize))

        self.sphereSizeLabel = tk.Label(self.inputOptionsFrame)
        self.sphereSizeLabel.place(relx=0.052, rely=0.772, height=15, width=98
                , bordermode='ignore')
        self.sphereSizeLabel.configure(activebackground="#f9f9f9")
        self.sphereSizeLabel.configure(anchor='w')
        self.sphereSizeLabel.configure(text='''Sphere Size:''')

        self.partTypeLabel = tk.Label(self.inputOptionsFrame)
        self.partTypeLabel.place(relx=0.517, rely=0.316, height=17, width=68
                , bordermode='ignore')
        self.partTypeLabel.configure(activebackground="#f9f9f9")
        self.partTypeLabel.configure(anchor='w')
        self.partTypeLabel.configure(text='''of type:''')

        self.partTypeEntry = tk.Entry(self.inputOptionsFrame)
        self.partTypeEntry.place(relx=0.724, rely=0.298, height=27
                , relwidth=0.228, bordermode='ignore')
        self.partTypeEntry.configure(background="white")
        self.partTypeEntry.configure(font="TkFixedFont")
        self.partTypeEntry.configure(selectbackground="#c4c4c4")
        self.partTypeEntry.configure(textvariable=self.partTypeStr)

        self.periodicCheck = tk.Checkbutton(self.inputOptionsFrame)
        self.periodicCheck.place(relx=0.034, rely=0.877, relheight=0.06
                , relwidth=0.638, bordermode='ignore')
        self.periodicCheck.configure(activebackground="#d9d9d9")
        self.periodicCheck.configure(anchor='w')
        self.periodicCheck.configure(justify='left')
        self.periodicCheck.configure(offrelief="sunken")
        self.periodicCheck.configure(text='''Periodic Boundaries''')
        self.periodicCheck.configure(variable=self.periodBound)

        self.framesFrame = tk.LabelFrame(self.inputOptionsFrame)
        self.framesFrame.place(relx=0.034, rely=0.421, relheight=0.263
                , relwidth=0.931, bordermode='ignore')
        self.framesFrame.configure(relief='groove')
        self.framesFrame.configure(text='''frames''')
        self.framesFrame.configure(width=270)

        self.fromEntry = tk.Entry(self.framesFrame)
        self.fromEntry.place(relx=0.185, rely=0.467, height=27, relwidth=0.133
                , bordermode='ignore')
        self.fromEntry.configure(background="white")
        self.fromEntry.configure(font="TkFixedFont")
        self.fromEntry.configure(selectbackground="#c4c4c4")
        self.fromEntry.configure(textvariable=self.fromFrameStr)
        self.fromEntry.bind('<FocusOut>',lambda e:pyxtal_support.validateInteger(e,
                                    self.fromFrameStr,
                                    self.fromFrame))

        self.fromLabel = tk.Label(self.framesFrame)
        self.fromLabel.place(relx=0.037, rely=0.533, height=15, width=38
                , bordermode='ignore')
        self.fromLabel.configure(activebackground="#f9f9f9")
        self.fromLabel.configure(anchor='w')
        self.fromLabel.configure(text='''From''')

        self.toLabel = tk.Label(self.framesFrame)
        self.toLabel.place(relx=0.37, rely=0.533, height=15, width=38
                , bordermode='ignore')
        self.toLabel.configure(activebackground="#f9f9f9")
        self.toLabel.configure(anchor='w')
        self.toLabel.configure(text='''to''')

        self.toEntry = tk.Entry(self.framesFrame)
        self.toEntry.place(relx=0.481, rely=0.467, height=27, relwidth=0.133
                , bordermode='ignore')
        self.toEntry.configure(background="white")
        self.toEntry.configure(font="TkFixedFont")
        self.toEntry.configure(selectbackground="#c4c4c4")
        self.toEntry.configure(textvariable=self.toFrameStr)
        self.toEntry.bind('<FocusOut>',lambda e:pyxtal_support.validateInteger(e,
                                    self.toFrameStr,
                                    self.toFrame))

        self.byLabel = tk.Label(self.framesFrame)
        self.byLabel.place(relx=0.704, rely=0.533, height=15, width=38
                , bordermode='ignore')
        self.byLabel.configure(activebackground="#f9f9f9")
        self.byLabel.configure(anchor='w')
        self.byLabel.configure(text='''by''')

        self.byEntry = tk.Entry(self.framesFrame)
        self.byEntry.place(relx=0.815, rely=0.467, height=27, relwidth=0.133
                , bordermode='ignore')
        self.byEntry.configure(background="white")
        self.byEntry.configure(font="TkFixedFont")
        self.byEntry.configure(selectbackground="#c4c4c4")
        self.byEntry.configure(textvariable=self.byFrameStr)
        self.byEntry.bind('<FocusOut>',lambda e:pyxtal_support.validateInteger(e,
                                    self.byFrameStr,
                                    self.byFrame))

        self.darkSpheresCheck = tk.Checkbutton(self.inputOptionsFrame)
        self.darkSpheresCheck.place(relx=0.517, rely=0.105, relheight=0.06
                , relwidth=0.466, bordermode='ignore')
        self.darkSpheresCheck.configure(activebackground="#d9d9d9")
        self.darkSpheresCheck.configure(anchor='w')
        self.darkSpheresCheck.configure(justify='left')
        self.darkSpheresCheck.configure(text='''Dark Spheres''')
        self.darkSpheresCheck.configure(variable=self.darkSpheres)

        self.loadDefButton = tk.Button(top)
        self.loadDefButton.place(relx=0.05, rely=0.041, height=33, width=108)
        self.loadDefButton.configure(activebackground="#d9d9d9")
        self.loadDefButton.configure(command=pyxtal_support.loadButtonCommand)
        self.loadDefButton.configure(text='''Load Defaults''')

        self.outputFrame = tk.LabelFrame(top)
        self.outputFrame.place(relx=0.532, rely=0.027, relheight=0.254
                , relwidth=0.449)
        self.outputFrame.configure(relief='groove')
        self.outputFrame.configure(text='''Output File Options''')
        self.outputFrame.configure(width=270)

        self.imageSizeEntry = tk.Entry(self.outputFrame)
        self.imageSizeEntry.place(relx=0.593, rely=0.595, height=27
                , relwidth=0.281, bordermode='ignore')
        self.imageSizeEntry.configure(background="white")
        self.imageSizeEntry.configure(font="TkFixedFont")
        self.imageSizeEntry.configure(selectbackground="#c4c4c4")
        self.imageSizeEntry.configure(textvariable=self.imageSizeStr)
        self.imageSizeEntry.bind('<FocusOut>',lambda e:pyxtal_support.validateInteger(e,
                                    self.imageSizeStr,
                                    self.imageSize))

        self.imageSizeLabel = tk.Label(self.outputFrame)
        self.imageSizeLabel.place(relx=0.056, rely=0.649, height=15, width=128
                , bordermode='ignore')
        self.imageSizeLabel.configure(activebackground="#f9f9f9")
        self.imageSizeLabel.configure(anchor='w')
        self.imageSizeLabel.configure(text='''Output Image size:''')

        self.outTriangCheck = tk.Checkbutton(self.outputFrame)
        self.outTriangCheck.place(relx=0.037, rely=0.162, relheight=0.092
                , relwidth=0.648, bordermode='ignore')
        self.outTriangCheck.configure(activebackground="#d9d9d9")
        self.outTriangCheck.configure(anchor='w')
        self.outTriangCheck.configure(justify='left')
        self.outTriangCheck.configure(text='''Image and Circles''')
        self.outTriangCheck.configure(variable=self.outCircles)

        self.outTriangCheck = tk.Checkbutton(self.outputFrame)
        self.outTriangCheck.place(relx=0.037, rely=0.324, relheight=0.092
                , relwidth=0.648, bordermode='ignore')
        self.outTriangCheck.configure(activebackground="#d9d9d9")
        self.outTriangCheck.configure(anchor='w')
        self.outTriangCheck.configure(justify='left')
        self.outTriangCheck.configure(text='''Triangulation''')
        self.outTriangCheck.configure(variable=self.outTriang)

        self.outAllCheck = tk.Checkbutton(self.outputFrame)
        self.outAllCheck.place(relx=0.037, rely=0.486, relheight=0.092
                , relwidth=0.759, bordermode='ignore')
        self.outAllCheck.configure(activebackground="#d9d9d9")
        self.outAllCheck.configure(anchor='w')
        self.outAllCheck.configure(justify='left')
        self.outAllCheck.configure(text='''Image + angle + defects''')
        self.outAllCheck.configure(variable=self.outAll)

        self.outMpegCheck = tk.Checkbutton(self.outputFrame)
        self.outMpegCheck.place(relx=0.037, rely=0.811, relheight=0.092
                , relwidth=0.389, bordermode='ignore')
        self.outMpegCheck.configure(activebackground="#d9d9d9")
        self.outMpegCheck.configure(anchor='w')
        self.outMpegCheck.configure(justify='left')
        self.outMpegCheck.configure(text='''mpeg''')
        self.outMpegCheck.configure(variable=self.outMpeg)

        self.analysisFrame = tk.LabelFrame(top)
        self.analysisFrame.place(relx=0.532, rely=0.302, relheight=0.337
                , relwidth=0.449)
        self.analysisFrame.configure(relief='groove')
        self.analysisFrame.configure(text='''Analysis''')
        self.analysisFrame.configure(width=270)

        self.trajCheck = tk.Checkbutton(self.analysisFrame)
        self.trajCheck.place(relx=0.037, rely=0.367, relheight=0.069
                , relwidth=0.759, bordermode='ignore')
        self.trajCheck.configure(activebackground="#d9d9d9")
        self.trajCheck.configure(anchor='w')
        self.trajCheck.configure(justify='left')
        self.trajCheck.configure(text='''Sphere Trajectories''')
        self.trajCheck.configure(variable=self.doTraject)
        self.trajCheck.configure(width=205)

        self.orientCorrCheck = tk.Checkbutton(self.analysisFrame)
        self.orientCorrCheck.place(relx=0.037, rely=0.735, relheight=0.069
                , relwidth=0.907, bordermode='ignore')
        self.orientCorrCheck.configure(activebackground="#d9d9d9")
        self.orientCorrCheck.configure(anchor='w')
        self.orientCorrCheck.configure(justify='left')
        self.orientCorrCheck.configure(text='''Orientational Correlation func.''')
        self.orientCorrCheck.configure(variable=self.doOrientCorr)
        self.orientCorrCheck.configure(width=245)

        self.orientHistCheck = tk.Checkbutton(self.analysisFrame)
        self.orientHistCheck.place(relx=0.037, rely=0.612, relheight=0.069
                , relwidth=0.833, bordermode='ignore')
        self.orientHistCheck.configure(activebackground="#d9d9d9")
        self.orientHistCheck.configure(anchor='w')
        self.orientHistCheck.configure(justify='left')
        self.orientHistCheck.configure(text='''Area Orientation Histogram''')
        self.orientHistCheck.configure(variable=self.doOrientHist)
        self.orientHistCheck.configure(width=225)

        self.sphereStatsCheck = tk.Checkbutton(self.analysisFrame)
        self.sphereStatsCheck.place(relx=0.037, rely=0.245, relheight=0.069
                , relwidth=0.722, bordermode='ignore')
        self.sphereStatsCheck.configure(activebackground="#d9d9d9")
        self.sphereStatsCheck.configure(anchor='w')
        self.sphereStatsCheck.configure(justify='left')
        self.sphereStatsCheck.configure(text='''Sphere Statistics''')
        self.sphereStatsCheck.configure(variable=self.doSphereStats)
        self.sphereStatsCheck.configure(width=195)

        self.defectStatsCheck = tk.Checkbutton(self.analysisFrame)
        self.defectStatsCheck.place(relx=0.037, rely=0.49, relheight=0.069
                , relwidth=0.796, bordermode='ignore')
        self.defectStatsCheck.configure(activebackground="#d9d9d9")
        self.defectStatsCheck.configure(anchor='w')
        self.defectStatsCheck.configure(justify='left')
        self.defectStatsCheck.configure(text='''Chrystal Defect Statistics''')
        self.defectStatsCheck.configure(variable=self.doDefectStats)
        self.defectStatsCheck.configure(width=215)

        self.zProfileCheck = tk.Checkbutton(self.analysisFrame)
        self.zProfileCheck.place(relx=0.037, rely=0.122, relheight=0.069
                , relwidth=0.722, bordermode='ignore')
        self.zProfileCheck.configure(activebackground="#d9d9d9")
        self.zProfileCheck.configure(anchor='w')
        self.zProfileCheck.configure(justify='left')
        self.zProfileCheck.configure(text='''Depth Profile (z axis)''')
        self.zProfileCheck.configure(variable=self.doZProfile)

        self.meaningLifeCheck = tk.Checkbutton(self.analysisFrame)
        self.meaningLifeCheck.place(relx=0.037, rely=0.857, relheight=0.069
                , relwidth=0.907, bordermode='ignore')
        self.meaningLifeCheck.configure(activebackground="#d9d9d9")
        self.meaningLifeCheck.configure(anchor='w')
        self.meaningLifeCheck.configure(justify='left')
        self.meaningLifeCheck.configure(text='''Meaning of Life''')
        self.meaningLifeCheck.configure(variable=self.doMeaningLife)

        self.windowFrame = tk.LabelFrame(top)
        self.windowFrame.place(relx=0.532, rely=0.659, relheight=0.199
                , relwidth=0.449)
        self.windowFrame.configure(relief='groove')
        self.windowFrame.configure(text='''Window Control''')
        self.windowFrame.configure(width=270)

        self.retainCheck = tk.Checkbutton(self.windowFrame)
        self.retainCheck.place(relx=0.037, rely=0.414, relheight=0.117
                , relwidth=0.441, bordermode='ignore')
        self.retainCheck.configure(activebackground="#d9d9d9")
        self.retainCheck.configure(anchor='w')
        self.retainCheck.configure(justify='left')
        self.retainCheck.configure(text='''Retain Windows''')
        self.retainCheck.configure(variable=self.retainWin)

        self.goButton = tk.Button(self.windowFrame)
        self.goButton.place(relx=0.037, rely=0.621, height=43, width=108
                , bordermode='ignore')
        self.goButton.configure(activebackground="#d9d9d9")
        self.goButton.configure(command=pyxtal_support.GoButtonCommand)
        self.goButton.configure(font=font10)
        self.goButton.configure(text='''Go''')

        self.killAllButton = tk.Button(self.windowFrame)
        self.killAllButton.place(relx=0.556, rely=0.621, height=43, width=108
                , bordermode='ignore')
        self.killAllButton.configure(activebackground="#d9d9d9")
        self.killAllButton.configure(command=pyxtal_support.killAllButtonCommand)
        self.killAllButton.configure(text='''Kill All 
Image Windows''')

        self.lockZoomCheck = tk.Checkbutton(self.windowFrame)
        self.lockZoomCheck.place(relx=0.556, rely=0.414, relheight=0.117
                , relwidth=0.404, bordermode='ignore')
        self.lockZoomCheck.configure(activebackground="#d9d9d9")
        self.lockZoomCheck.configure(anchor='w')
        self.lockZoomCheck.configure(justify='left')
        self.lockZoomCheck.configure(text='''Lock zoom''')
        self.lockZoomCheck.configure(variable=self.lockZoom)

        self.lockViewsCheck = tk.Checkbutton(self.windowFrame)
        self.lockViewsCheck.place(relx=0.556, rely=0.207, relheight=0.117
                , relwidth=0.404, bordermode='ignore')
        self.lockViewsCheck.configure(activebackground="#d9d9d9")
        self.lockViewsCheck.configure(anchor='w')
        self.lockViewsCheck.configure(justify='left')
        self.lockViewsCheck.configure(text='''Lock views''')
        self.lockViewsCheck.configure(variable=self.lockViews)

        self.batchmodeCheck = tk.Checkbutton(self.windowFrame)
        self.batchmodeCheck.place(relx=0.037, rely=0.207, relheight=0.117
                , relwidth=0.478, bordermode='ignore')
        self.batchmodeCheck.configure(activebackground="#d9d9d9")
        self.batchmodeCheck.configure(anchor='w')
        self.batchmodeCheck.configure(justify='left')
        self.batchmodeCheck.configure(text='''Quiet Batchmode''')
        self.batchmodeCheck.configure(variable=self.batchmode)
        self.batchmodeCheck.configure(width=129)

        self.progressFrame = tk.LabelFrame(top)
        self.progressFrame.place(relx=0.532, rely=0.879, relheight=0.103
                , relwidth=0.449)
        self.progressFrame.configure(relief='groove')
        self.progressFrame.configure(text='''Progress''')
        self.progressFrame.configure(width=270)

        self.fileProgressbar = ttk.Progressbar(self.progressFrame)
        self.fileProgressbar.place(relx=0.556, rely=0.267, relwidth=0.407
                , relheight=0.0, height=19, bordermode='ignore')

        self.fileMessage = tk.Message(self.progressFrame)
        self.fileMessage.place(relx=0.037, rely=0.267, relheight=0.227
                , relwidth=0.485, bordermode='ignore')
        self.fileMessage.configure(anchor='w')
        self.fileMessage.configure(text='''Files: ''')
        self.fileMessage.configure(width=131)

        self.frameMessage = tk.Message(self.progressFrame)
        self.frameMessage.place(relx=0.037, rely=0.667, relheight=0.227
                , relwidth=0.559, bordermode='ignore')
        self.frameMessage.configure(anchor='w')
        self.frameMessage.configure(text='''Frames: ''')
        self.frameMessage.configure(width=151)

        self.frameProgressbar = ttk.Progressbar(self.progressFrame)
        self.frameProgressbar.place(relx=0.556, rely=0.667, relwidth=0.407
                , relheight=0.0, height=19, bordermode='ignore')

# The following code is added to facilitate the Scrolled widgets you specified.
class AutoScroll(object):
    '''Configure the scrollbars for a widget.'''

    def __init__(self, master):
        #  Rozen. Added the try-except clauses so that this class
        #  could be used for scrolled entry widget for which vertical
        #  scrolling is not supported. 5/7/14.
        try:
            vsb = ttk.Scrollbar(master, orient='vertical', command=self.yview)
        except:
            pass
        hsb = ttk.Scrollbar(master, orient='horizontal', command=self.xview)

        #self.configure(yscrollcommand=_autoscroll(vsb),
        #    xscrollcommand=_autoscroll(hsb))
        try:
            self.configure(yscrollcommand=self._autoscroll(vsb))
        except:
            pass
        self.configure(xscrollcommand=self._autoscroll(hsb))

        self.grid(column=0, row=0, sticky='nsew')
        try:
            vsb.grid(column=1, row=0, sticky='ns')
        except:
            pass
        hsb.grid(column=0, row=1, sticky='ew')

        master.grid_columnconfigure(0, weight=1)
        master.grid_rowconfigure(0, weight=1)

        # Copy geometry methods of master  (taken from ScrolledText.py)
        if py3:
            methods = tk.Pack.__dict__.keys() | tk.Grid.__dict__.keys() \
                  | tk.Place.__dict__.keys()
        else:
            methods = tk.Pack.__dict__.keys() + tk.Grid.__dict__.keys() \
                  + tk.Place.__dict__.keys()

        for meth in methods:
            if meth[0] != '_' and meth not in ('config', 'configure'):
                setattr(self, meth, getattr(master, meth))

    @staticmethod
    def _autoscroll(sbar):
        '''Hide and show scrollbar as needed.'''
        def wrapped(first, last):
            first, last = float(first), float(last)
            if first <= 0 and last >= 1:
                sbar.grid_remove()
            else:
                sbar.grid()
            sbar.set(first, last)
        return wrapped

    def __str__(self):
        return str(self.master)

def _create_container(func):
    '''Creates a ttk Frame with a given master, and use this new frame to
    place the scrollbars and the widget.'''
    def wrapped(cls, master, **kw):
        container = ttk.Frame(master)
        container.bind('<Enter>', lambda e: _bound_to_mousewheel(e, container))
        container.bind('<Leave>', lambda e: _unbound_to_mousewheel(e, container))
        return func(cls, container, **kw)
    return wrapped

class ScrolledListBox(AutoScroll, tk.Listbox):
    '''A standard Tkinter Text widget with scrollbars that will
    automatically show/hide as needed.'''
    @_create_container
    def __init__(self, master, **kw):
        tk.Listbox.__init__(self, master, **kw)
        AutoScroll.__init__(self, master)

import platform
def _bound_to_mousewheel(event, widget):
    child = widget.winfo_children()[0]
    if platform.system() == 'Windows' or platform.system() == 'Darwin':
        child.bind_all('<MouseWheel>', lambda e: _on_mousewheel(e, child))
        child.bind_all('<Shift-MouseWheel>', lambda e: _on_shiftmouse(e, child))
    else:
        child.bind_all('<Button-4>', lambda e: _on_mousewheel(e, child))
        child.bind_all('<Button-5>', lambda e: _on_mousewheel(e, child))
        child.bind_all('<Shift-Button-4>', lambda e: _on_shiftmouse(e, child))
        child.bind_all('<Shift-Button-5>', lambda e: _on_shiftmouse(e, child))

def _unbound_to_mousewheel(event, widget):
    if platform.system() == 'Windows' or platform.system() == 'Darwin':
        widget.unbind_all('<MouseWheel>')
        widget.unbind_all('<Shift-MouseWheel>')
    else:
        widget.unbind_all('<Button-4>')
        widget.unbind_all('<Button-5>')
        widget.unbind_all('<Shift-Button-4>')
        widget.unbind_all('<Shift-Button-5>')

def _on_mousewheel(event, widget):
    if platform.system() == 'Windows':
        widget.yview_scroll(-1*int(event.delta/120),'units')
    elif platform.system() == 'Darwin':
        widget.yview_scroll(-1*int(event.delta),'units')
    else:
        if event.num == 4:
            widget.yview_scroll(-1, 'units')
        elif event.num == 5:
            widget.yview_scroll(1, 'units')

def _on_shiftmouse(event, widget):
    if platform.system() == 'Windows':
        widget.xview_scroll(-1*int(event.delta/120), 'units')
    elif platform.system() == 'Darwin':
        widget.xview_scroll(-1*int(event.delta), 'units')
    else:
        if event.num == 4:
            widget.xview_scroll(-1, 'units')
        elif event.num == 5:
            widget.xview_scroll(1, 'units')

if __name__ == '__main__':
    vp_start_gui()





