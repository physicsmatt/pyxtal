import tkinter as tk

class pyxtal_app_main(tk.Tk):
    def __init__(self, master):
        self.master = master
        self.ctrl_frame = tk.Frame(self.master)
        self.lab1 = tk.Label(self.ctrl_frame,text="Hello again I am a Label").pack()
        self.button1 = tk.Button(self.ctrl_frame, text = 'New Window', width = 25, command = self.go_button_hit)
        self.button1.pack()
        self.howmany = 527
        self.ctrl_frame.pack()
        self.animal = tk.IntVar()
        self.catbut = tk.Radiobutton(self.ctrl_frame,text="cat",variable=self.animal,value=1).pack()
        self.dogbut = tk.Radiobutton(self.ctrl_frame,text="dog",variable=self.animal,value=2).pack()

#    def go_button_hit(self):
#        #get list of filenames
#        for i in range(0,1):
#            self.new_window()
#
#    def new_window(self):
#        self.newWindow = tk.Toplevel(self.master)
#        self.app = pyxtal_win(self.newWindow, self)

    def go_button_hit(self):
        #get list of filenames
        print("hello1")
        for i in range(0,1):
            self.app = pyxtal_win(tk.Toplevel(self.master), self)
            self.app.msg.config(text = "Now who da boss!")


#    def new_window(self):
#        self.newWindow = tk.Toplevel(self.master)
#        self.app = pyxtal_win(tk.Toplevel(self.master), self)

#class pyxtal_win:
class pyxtal_win(tk.Tk):
    def __init__(self, master, parent):
        self.master = master
        self.ctrl_frame = tk.Frame(self.master)
        self.nowhowmany = parent.howmany
        self.msg = tk.Message(self.ctrl_frame, text = "msgtext = "+str(self.nowhowmany))
        self.msg.config(text =  "animal is " + str(parent.animal.get()))
        self.msg.pack()


        self.quitButton = tk.Button(self.ctrl_frame, text = 'Quit', width = 25, command = self.close_windows)
        self.quitButton.pack()
        self.ctrl_frame.pack()

    def close_windows(self):
        self.master.destroy()
       

root = tk.Tk()
app = pyxtal_app_main(root)
root.mainloop()

