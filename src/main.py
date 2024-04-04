# The purpose of this file is to build the GUI application with the tkinter module,
# and to call on other files within the 'data' directory to use the functionalities of the
# application. 

import tkinter as tk # Import the tkinter module in order to build the GUI
from tkinter import *

window = tk.Tk() # Initiate a window

w = tk.Label (window, text="Hello, CW4")
w.pack()

# Sets different visual attributes about the window, such as title, size, etc.
window.title("CW4 Dr. Jolley Project")
window.geometry("400x350+500+500")

window.mainloop()
