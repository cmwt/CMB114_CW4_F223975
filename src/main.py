# The purpose of this file is to build the GUI application with the tkinter module,
# and to call on other files within the 'data' directory to use the functionalities of the
# application. 

import tkinter as tk # Import the tkinter module in order to build the GUI
from tkinter import *

window = tk.Tk() # Initiate a window

def show():
    label.config( text = Element.get() )

options = [
    "Carbon",
    "Hydrogen",
    "Oxygen",
    "Sulfur"
]

element = StringVar()
element.set("Carbon")
drop = OptionMenu(window, element, *options)
drop.pack(side = "top")

# The following code formats the file menu at the top of the GUI window
menubar = Menu(window)

# Sets different visual attributes about the window, such as title, size, etc.
window.title("CW4 Project - Molecular Dynamics")
window.attributes('-fullscreen', True)

window.mainloop()
