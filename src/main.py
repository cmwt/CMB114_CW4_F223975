# The purpose of this file is to build the GUI application with the tkinter module,
# and to call on other files within the 'data' directory to use the functionalities of the
# application. 

import tkinter as tk # Import the tkinter module in order to build the GUI
from tkinter import *

window = tk.Tk() # Initiate a window

# The following code formats the file menu at the top of the GUI window
menubar = Menu(window)

filemenu = Menu(menubar,tearoff=0)

# Add commands to menu
filemenu.add_command(label="New File")
filemenu.add_command(label="Open")
filemenu.add_command(label="Save")
menubar.add_cascade(label="File", menu=filemenu)
window.config(menu=menubar) 

# Creates a space for a toolbar at the top of the GUI window.
toolbar = tk.Frame(window, background="light grey", height = 30) # Defines a toolbar at the top of the window
toolbar.pack(side="top", fill="x") # Inserts the toolbar into the window

# The following code assigns a format for the bulk of the GUI window
main = tk.PanedWindow(window, background="light grey") # Adds a section in the window
main.pack(side='top', fill='both', expand=True) # Packs the paned window

# Creates a paned section in the bulk of the GUI window
left_pane = tk.Frame(main, background="dark grey", width=200)
right_pane = tk.PanedWindow(main, background="dark grey", width=200)
main.add(left_pane)
main.add(right_pane)

# Sets different visual attributes about the window, such as title, size, etc.
window.title("CW4 Project - Molecular Dynamics")
window.attributes('-fullscreen', True)

window.mainloop()
