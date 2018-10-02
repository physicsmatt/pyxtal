import tkinter as tk
import PIL, PIL.Image, PIL.ImageTk
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


class pyxtal_app_main(tk.Tk):
    def __init__(self, master):
        self.master = master
        self.ctrl_frame = tk.Frame(self.master)
        self.lab1 = tk.Label(self.ctrl_frame,text="Hello again I am a Label").pack()
        self.button1 = tk.Button(self.ctrl_frame, text = 'New Window', width = 25, command = self.go_button_hit)
        self.button1.pack()
        self.howmany = 527
        self.ctrl_frame.pack()

    def go_button_hit(self):
        #get list of filenames
        print("hello1")
        for i in range(0,2):
            self.app = pyxtal_win(tk.Toplevel(self.master), self)
            self.app.msg.config(text = "Now who da boss!")

#    def new_window(self):
#        self.newWindow = tk.Toplevel(self.master)
#        self.app = pyxtal_win(tk.Toplevel(self.master), self)


class pyxtal_win(tk.Tk):
    def __init__(self, master, parent):
        self.master = master
        self.ctrl_frame = tk.Frame(self.master)
        self.img_frame = tk.Frame(self.master)
        self.nowhowmany = parent.howmany
        self.msg = tk.Message(self.ctrl_frame, text = "msgtext = "+str(self.nowhowmany))
        self.msg.config(text =  "animal is " + str(parent.howmany))
        self.msg.pack()
        self.animal = tk.IntVar(value = 1)
        self.catbut = tk.Radiobutton(self.ctrl_frame,text="cat",
                                     command=self.update_img,variable=self.animal,value=1).pack()
        self.dogbut = tk.Radiobutton(self.ctrl_frame,text="dog",
                                     command=self.update_img,variable=self.animal,value=2).pack()
        self.quitButton = tk.Button(self.ctrl_frame, text = 'Quit', width = 25, command = self.close_windows)
        self.quitButton.pack()
        self.ctrl_frame.pack()

        self.imgcanv = tk.Canvas(self.img_frame)
        self.catimg = PIL.ImageTk.PhotoImage(PIL.Image.open("cat.jpg"))
        self.dogimg = PIL.ImageTk.PhotoImage(PIL.Image.open("dog.jpg"))
        self.img_on_canv = self.imgcanv.create_image(100,100,image=self.catimg)
        self.imgcanv.pack()
        self.img_frame.pack()


#        self.img = plt.imread("cat.jpg")
#        fig2, ax = plt.subplots()
#        ax.imshow(self.img)

#        fig = mpl.figure.Figure(figsize=(5, 4), dpi=100)

        fig,ax = plt.subplots()
        x = range(300)
        ax.plot(x, x, '--', linewidth=5, color='firebrick',zorder=1)
        thecat = plt.imread("cat.jpg")
        thedog = plt.imread("dog.jpg")

        aximg = ax.imshow(thecat, extent=[0, 400, 0, 300],zorder=0)
        #aximg.set_data(thedog)
        x2 = range(20,200,10)
        aximg2 = ax.imshow(thedog, extent=[0, 400, 0, 300],zorder=0.5, alpha=0.5) 
        blueplot = ax.scatter(x2, x2, color='blue',zorder=2)
        blueplot.set_visible(1)
        ax.axis([100,200,110,210])
        ax.axis('off')


#        ax.imshow(plt.imread("cat.jpg"))
#        t = np.arange(0, 3, .01)
#        fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
#        fig.add_subplot(111).plot(t, 2 * np.cos(2 * np.pi * t))


        

        self.canvas = FigureCanvasTkAgg(fig, master=self.master)  # A tk.DrawingArea.
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.imgcanv.bind("<Button-4>", self._on_mousewheel)
#        self.imgcanv.bind("<MouseWheel>", self._on_mousewheel)
#        self.imgcanv.bind("<Button-1>", self._on_mousewheel)

        
#        toolbar = NavigationToolbar2Tk(self.canvas, self.master)
#        toolbar.update()
#        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


#        self.canvas.mpl_connect("key_press_event", on_key_press)

    def _on_mousewheel(self, event):
        print("mouse wheel!")
        

    def update_img(self):
        if self.animal.get() == 1:
           self.imgcanv.itemconfig(self.img_on_canv, image=self.catimg)
        if self.animal.get() == 2:
           self.imgcanv.itemconfig(self.img_on_canv, image=self.dogimg)

    def close_windows(self):
        self.master.quit()
        self.master.destroy()

print(mpl.__version__)
root = tk.Tk()
app = pyxtal_app_main(root)
root.mainloop()

