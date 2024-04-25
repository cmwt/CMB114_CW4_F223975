# The purpose of this file is to build the GUI application with the tkinter module,
# and to call on other files within the 'data' directory to use the functionalities of the
# application. 

molecules = [
    "Benzene",
    "Naphthalene",
    "Cyclohexane",
    "Acetone"
]

import tkinter as tk # Import the tkinter module in order to build the GUI
from tkinter import ttk

root = tk.Tk() # Initiate a window
root.geometry("1500x750")

# ****** Toolbar ******

toolbar = ttk.Frame(root, borderwidth = 5, relief = tk.GROOVE)
toolbar.pack(side="top", fill = "x")
select_but = ttk.Button(toolbar, text="Select")
select_but.pack(side="left", padx = 2, pady = 2)
nav_but = ttk.Button(toolbar, text = "Navigate")
nav_but.pack(side = "left", padx = 2, pady = 2)

# ****** Left Pane ******

frame = ttk.Frame(root, width = 300, height = 750, borderwidth = 5, relief = tk.GROOVE)
frame.pack_propagate(False)
frame.pack(side="left", fill = "y")

drop = ttk.Combobox(frame, value= molecules)
drop.place(x = 0, y = 0)

# ****** Right Frame ******

frame1 = ttk.Frame(root, width = 1200, height = 750, borderwidth = 5, relief = tk.GROOVE)
frame1.pack_propagate(False)
frame1.pack(side="left", fill = "both")

# Sets different visual attributes about the window, such as title, size, etc.
root.title("CW4 Project - Molecular Dynamics")

root.mainloop()
