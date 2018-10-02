import tkinter as tk

def create_pyxtal_app():
    paroot = tk.Tk()
    ctrl_frame = tk.Frame(paroot).pack()
    lab = tk.Label(ctrl_frame,text="Hello again I am a Label").pack()
    msgtext = "This is a message widget"
    msg = tk.Message(ctrl_frame, text = msgtext)
    msg.pack()
# Note about pack(): it looks like you can either use it as .pack on creation, or as a separate statement,
# but not both.
    paroot.animal=tk.IntVar()
    tk.Radiobutton(ctrl_frame,text="cat",variable=animal,command=refresh_pic(paroot),value=1).pack()
    tk.Radiobutton(ctrl_frame,text="dog",variable=animal,command=refresh_pic(paroot),value=2).pack()
    msg.config(text =  "animal is " + str(animal.get()))
    #msg.pack
    print,("random thought.")
    return(paroot)

def refresh_pic(root):
   root.Title=("its " + str(animal.get()))


#root = tk.Tk()

#v = tk.IntVar()

#tk.Label(root, 
#        text="""Choose a 
#programming language:""",
#        justify = tk.LEFT,
#        padx = 20).pack()
#tk.Radiobutton(root, 
#              text="Python",
#              padx = 20, 
#              variable=v, 
#              value=1).pack(anchor=tk.W)
#tk.Radiobutton(root, 
#              text="Perl",
#              padx = 20, 
#              variable=v, 
#              value=2).pack(anchor=tk.W)

pyxtal_app = create_pyxtal_app()

pyxtal_app.mainloop()


