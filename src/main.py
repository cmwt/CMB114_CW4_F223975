# The purpose of this file is to build the GUI application with the tkinter module,
# and to call on other files within the 'data' directory to use the functionalities of the
# application. 

elements = [
    "Carbon",
    "Hydrogen",
    "Oxygen",
    "Sulfur"
]

import tkinter as tk # Import the tkinter module in order to build the GUI
from tkinter import *

root = tk.Tk() # Initiate a window
root.geometry("1500x750")

frame = tk.Frame(root, width = 300, height = 750, borderwidth = 5, relief = tk.GROOVE)
frame.pack_propagate(False)
frame.pack(side="left")

frame1 = tk.Frame(root, width = 1200, height = 750, borderwidth = 5, relief = tk.GROOVE)
frame1.pack_propagate(False)
frame1.pack(side="right")

element = StringVar()
element.set("Element")
drop = OptionMenu(frame, element, *elements)
drop.pack()
         
# Sets different visual attributes about the window, such as title, size, etc.
root.title("CW4 Project - Molecular Dynamics")

root.mainloop()
