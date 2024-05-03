# The purpose of this file is to build the GUI application with the tkinter module,
# and to call on other files within the 'data' directory to use the functionalities of the
# application. 

molecules = ["Benzene", "Naphthalene", "Cyclohexane", "Acetone"]

# Defining a function to turn atom number button on/off

global atom_on # Sets button to "On" automatically when the GUI starts
atom_on = True

def atom_switch():
    global atom_on
    if atom_on:
       atom_num.config(text = "Off") # Switches the button to "Off" when clicked
       atom_on = False 
    else:
        atom_num.config(text = "On") # Switches the button to "On" when clicked again
        atom_on = True


import tkinter as tk # Import the tkinter module in order to build the GUI
from tkinter import ttk
from PIL import Image, ImageTk
from rdkit import Chem
from rdkit.Chem import Draw

root = tk.Tk() # Initiate a window
root.geometry("1500x750")
root.configure( background = "white")

# ***** Molecules *****

# Benzene
benzene_SMILES = "c1ccccc1" # Assigns benzene a SIMLES moleucle format
benzene = Chem.MolFromSmiles(benzene_SMILES) # Creates the SMILES molecule
benz_img = Draw.MolToImage(benzene) # Creates the benzene image
benz_photo = ImageTk.PhotoImage(benz_img) # Puts the benzene image onto a tkinter PhotoImage

# Naphthalene
naph_SMILES = "c1c2ccccc2ccc1"
naph = Chem.MolFromSmiles(naph_SMILES) 
naph_img = Draw.MolToImage(naph) 
naph_photo = ImageTk.PhotoImage(naph_img)

# Cyclohexane
cyhex_SMILES = "C1CCCCC1"
cyhex = Chem.MolFromSmiles(cyhex_SMILES) 
cyhex_img = Draw.MolToImage(cyhex) 
cyhex_photo = ImageTk.PhotoImage(cyhex_img)

# Acetone
ace_SMILES = "CC(=O)C"
ace = Chem.MolFromSmiles(ace_SMILES) 
ace_img = Draw.MolToImage(ace) 
ace_photo = ImageTk.PhotoImage(ace_img)

# ****** Toolbar ******

toolbar = ttk.Frame(root, borderwidth = 5,  relief = tk.GROOVE)
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
drop.pack(pady = 5)

atom_num = tk.Button(frame, text = "On", command=atom_switch)
atom_label = tk.Label(frame, text= "Atom Numbers:")
atom_label.pack(pady=2)
atom_num.pack(pady = 2)

# ****** Right Frame ******

frame1 = ttk.Frame(root, width = 1200, height = 750, borderwidth = 5, relief = tk.GROOVE)
frame1.pack_propagate(False)
frame1.pack(side="left", fill = "both")

benz_canvas = tk.Canvas(frame1)
benz_canvas.create_image(0, 0, anchor=tk.NW, image=benz_photo)
benz_canvas.grid(row=1, column=1)

naph_canvas = tk.Canvas(frame1)
naph_canvas.create_image(0, 0, anchor=tk.NW, image=naph_photo)
naph_canvas.grid(row=1, column=2)

cyhex_canvas = tk.Canvas(frame1)
cyhex_canvas.create_image(0, 0, anchor=tk.NW, image=cyhex_photo)
cyhex_canvas.grid(row=1, column=3)

ace_canvas = tk.Canvas(frame1)
ace_canvas.create_image(0, 0, anchor=tk.NW, image=ace_photo)
ace_canvas.grid(row=2, column=1)

root.title("CW4 Project - Molecular Dynamics")
frame1.config(row=3, column=3)

root.mainloop()
