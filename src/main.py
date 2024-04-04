# The purpose of this file is to build the GUI application with the tkinter module,
# and to call on other files within the 'data' directory to use the functionalities of the
# application. 

import tkinter as tk # Import the tkinter module in order to build the GUI
from tkinter import *

window = tk.Tk() # Initiate a window
toolbar = tk.Frame(window, background="light grey", height = 30) # Defines a toolbar at the top of the window
toolbar.pack(side="top", fill="x") # Inserts the toolbar into the window

# Sets different visual attributes about the window, such as title, size, etc.
window.title("CW4 Dr. Jolley Project")

window.mainloop()
