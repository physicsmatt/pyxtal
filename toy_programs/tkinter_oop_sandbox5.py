import tkinter as tk
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_agg import FigureCanvasAgg
#https://matplotlib.org/gallery/user_interfaces/embedding_in_tk_canvas_sgskip.html
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)


#This class defines the main window.  It will read options (from file or command line) and allow changes.
#Then user hits go, and image files are read and processed, each in a new window.
class pyxtal_app_main(tk.Tk):
    def __init__(self, master):
        self.master = master
        self.master.protocol("WM_DELETE_WINDOW", self._kill_main_window)
        self.ctrl_frame = tk.Frame(self.master)
        self.lab1 = tk.Label(self.ctrl_frame,text="Hello again I am a Label").pack()
        self.button1 = tk.Button(self.ctrl_frame, text = 'New Window', width = 25, command = self.go_button_hit)
        self.button1.pack()
        self.spheresize = 6
        self.ctrl_frame.pack()

    #see http://effbot.org/tkinterbook/tkinter-events-and-bindings.htm .  This little bit turned out to be needed.
    #Without it, window would close but process would still be going, not returning control to terminal.
    def _kill_main_window(self):
        self.master.quit()
        self.master.destroy()

    def go_button_hit(self):
        #get list of filenames
        print("hello1")
        for i in range(0,2):
            #would need to hang on to handle if I want to influence the app again:
            self.app = pyxtal_win(tk.Toplevel(self.master), self)
            self.app.msgfromparent.config(text = "This text is from parent.")

#    def new_window(self):
#        self.newWindow = tk.Toplevel(self.master)
#        self.app = pyxtal_win(tk.Toplevel(self.master), self)


class pyxtal_win(tk.Tk):
    def __init__(self, master, parent):
        self.master = master

      #start of controls area.
        self.ctrl_frame = tk.Frame(self.master, background='red',padx=12,pady=5)

        #These lines demonstrate how to get data from parent window
        self.localspheresize = parent.spheresize
        self.msg = tk.Message(self.ctrl_frame, text = "sphere size = "+str(self.localspheresize))
        self.msg.pack()

        #These lines demonstrate how parent can control this child window
        self.msgfromparent = tk.Message(self.ctrl_frame, text = "no message from parent yet")
        self.msgfromparent.pack()

        #define variable used by radiobuttons
        self.animal = tk.IntVar(value = 1)
        self.catbut = tk.Radiobutton(self.ctrl_frame,text="cat",
                                     command=self.update_img,variable=self.animal,value=1).pack()
        self.dogbut = tk.Radiobutton(self.ctrl_frame,text="dog",
                                     command=self.update_img,variable=self.animal,value=2).pack()


        #define button to quit this window
        self.quitButton = tk.Button(self.ctrl_frame, text = 'Quit', width = 25, command = self.close_windows)
        self.quitButton.pack()
        self.ctrl_frame.pack()

      #start of image area.

        self.img_frame = tk.Frame(self.master, background='blue',padx=20,pady=10)
        self.msg2 = tk.Message(self.img_frame, text= "hopefully in img_frame")
        self.msg2.pack()
        self.imgcanv = tk.Canvas(self.img_frame)

        fig,ax = plt.subplots()
        x = range(300)
        ax.plot(x, x, '--', linewidth=5, color='firebrick',zorder=1)
        thecat = plt.imread("cat.jpg")
        thedog = plt.imread("dog.jpg")

        self.catimg = ax.imshow(thecat, extent=[0, 400, 0, 300],zorder=0)
        #aximg.set_data(thedog)

        x2 = range(20,200,10)
        self.dogimg = ax.imshow(thedog, extent=[0, 400, 0, 300],zorder=0.5)#, alpha=0.5) 
        blueplot = ax.scatter(x2, x2, color='blue',zorder=2)
        blueplot.set_visible(1)
        #ax.axis([100,200,110,210])
        ax.axis('off')
        fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
        self.dogimg.set_visible(0)

        #Here's where we actually put the plot on the tk canvas:
        self.img_frame.canvas = FigureCanvasTkAgg(fig, master=self.img_frame)
        self.img_frame.canvas.draw()
        self.img_frame.canvas.get_tk_widget().pack()
        self.img_frame.pack()

      #Now set callbacks for mouse events (& eventually keyboard events.)

        self.img_frame.canvas.get_tk_widget().bind("<Button-4>", self._on_mousewheel)
        self.img_frame.canvas.get_tk_widget().bind("<Button-5>", self._on_mousewheel)
        self.master.bind("<Key>", self._on_keypress)
        #Not sure what the hell this is supposed to do:
        #self.canvas.mpl_connect("key_press_event", on_key_press)
        
        #These lines would keep the usual toolbar from matplotlib plots, which I don't want.
        #toolbar = NavigationToolbar2Tk(self.canvas, self.master)
        #toolbar.update()
        #self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


    def _on_mousewheel(self, event):
        print(event.x,event.y,event.delta,event.num)

    def _on_keypress(self, event):
        print("keypress: ", event.x,event.y,event.char,event.keysym,event.keycode)
        if event.char in ('d', 'D'):
             print("dog")
             self.animal.set(2)
        if event.char in ('c', 'C'):
             print("cat")
             self.animal.set(1)
        self.update_img()

    def update_img(self):
        if self.animal.get() == 1:
           self.dogimg.set_visible(0)
        if self.animal.get() == 2:
           self.dogimg.set_visible(1)
        self.img_frame.canvas.draw()

    def close_windows(self):
        #self.master.quit() #This line causes whole program to close, not just one window.
        self.master.destroy()

#print(mpl.__version__)
root = tk.Tk()
app = pyxtal_app_main(root)
root.mainloop()


