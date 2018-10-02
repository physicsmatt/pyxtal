import tkinter as tk
import PIL, PIL.Image, PIL.ImageTk

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
        for i in range(0,1):
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

    def update_img(self):
        if self.animal.get() == 1:
           self.imgcanv.itemconfig(self.img_on_canv, image=self.catimg)
        if self.animal.get() == 2:
           self.imgcanv.itemconfig(self.img_on_canv, image=self.dogimg)

    def close_windows(self):
        self.master.destroy()

root = tk.Tk()
app = pyxtal_app_main(root)
root.mainloop()

