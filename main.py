# The purpose of this file is to build the GUI application with the tkinter module,
# and to call on other files within the 'data' directory to use the functionalities of the
# application.

# ********************
# *                  *
# *     Imports      *
# *                  *
# ********************

import tkinter as tk 
from tkinter import ttk
from PIL import Image, ImageTk
from rdkit import Chem
from rdkit.Chem import Draw
from lib.mol import *

# **********************
# *                    *
# *     Functions      *
# *                    *
# **********************

# Function to print properties text in left frame
def reference_text(x):
    reference_label.config(text = str(text[x]) )

# Function to turn atom number button on/off
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

# Function to draw correct molecule and text when selected form dropdown
def option_selected(event):
    if drop.get() == "Benzene":
        draw_benz()
        reference_text(drop.get())
    if drop.get() == "Naphthalene":
        draw_naph()
        reference_text(drop.get())
    if drop.get() == "Cyclohexane":
        draw_cyhex()
        reference_text(drop.get())
    if drop.get() == "Acetone":
        draw_ace()
        reference_text(drop.get())

# Pre-programmed molecule drawing functions
def draw_benz():
    benz_canvas = tk.Canvas(frame1, width = benz_photo.width(), height = benz_photo.height())
    benz_canvas.create_image(0, 0, anchor=tk.NW, image=benz_photo)
    benz_canvas.grid(row=1, column=1)

def draw_naph():
    naph_canvas = tk.Canvas(frame1, width=naph_photo.width(), height=naph_photo.height())
    naph_canvas.create_image(0, 0, anchor=tk.NW, image=naph_photo)
    naph_canvas.grid(row=1, column=2)

def draw_cyhex():
    cyhex_canvas = tk.Canvas(frame1, width=cyhex_photo.width(), height=cyhex_photo.height())
    cyhex_canvas.create_image(0, 0, anchor=tk.NW, image=cyhex_photo)
    cyhex_canvas.grid(row=1, column=3)

def draw_ace():
    ace_canvas = tk.Canvas(frame1, width=ace_photo.width(), height=ace_photo.height())
    ace_canvas.create_image(0, 0, anchor=tk.NW, image=ace_photo)
    ace_canvas.grid(row=2, column=1)

def draw_custom():
    mol = Chem.MolFromSmiles(str(custom_entry.get()))
    custom_img = Draw.MolToImage(mol)
    custom_photo = ImageTk.PhotoImage(custom_img)
    custom_canvas = tk.Canvas(frame1, width = custom_photo.width(), height = custom_photo.height())
    custom_canvas.create_image(0, 0, anchor=tk.NW, image=custom_photo)
    custom_canvas.grid(row=2, column=2)

# Clear function erase contents of right frame
def clear_frame():
    for widgets in frame1.winfo_children():
        widgets.destroy()

# ***************************
# *                         *
# *     Tkinter Window      *
# *                         *
# ***************************

root = tk.Tk()
root.geometry("570x750")
root.configure( background = "white")
root.title("CW4 Project - Molecular Dynamics")

molecules = ["Benzene", "Naphthalene", "Cyclohexane", "Acetone"]

# *** Left Frame ***

frame = ttk.Frame(root, width = 300, height = 750, borderwidth = 5, relief = tk.GROOVE)
frame.grid(row=1, column=1, sticky="NW")

drop_label = tk.Label(frame, text = "Pre-Programmed Molecules:")
drop_label.grid(row=1, column=1, sticky="W")
drop = ttk.Combobox(frame, value= molecules)
drop.grid(row=2, column=1, sticky="W")
drop.bind("<<ComboboxSelected>>", option_selected)

clear_but = tk.Button(frame, text="Clear", command=clear_frame)
clear_but.grid(row=2, column=2, sticky="W")

custom_label = tk.Label(frame, text="Custom SMILES code:")
custom_label.grid(row=3, column=1, sticky="W")

custom_entry = tk.Entry(frame)
custom_entry.grid(row=4, column=1, sticky="W")

draw_button = tk.Button(frame, text="Draw", command=lambda:[draw_custom()])
draw_button.grid(row=4, column=2, sticky="W")

atom_label = tk.Label(frame, text= "Atom Numbers:")
atom_label.grid(row=5, column=1, sticky="W")
atom_num = tk.Button(frame, text = "On", command=atom_switch)
atom_num.grid(row=5, column=2, sticky="W")

reference_label = tk.Label(frame)
reference_label.grid(row=6, column=1, sticky="W")


text = {"Benzene":"Molar Mass: 78.11 g/mol\nBoiling Point: 80.1°C\nMelting Point: 5.5°C",
        "Naphthalene":"Molar Mass: 128.17 g/mol",
                  }

# *** Right Frame ***

frame1 = ttk.Frame(root, width = 1200, height = 750, borderwidth = 5, relief = tk.GROOVE)
frame1.grid(row=1, column=2, sticky="NW")

# *** Molecule images ***

benz_img = Draw.MolToImage(benzene) # Creates the benzene image
benz_photo = ImageTk.PhotoImage(benz_img) # Puts the benzene image onto a tkinter PhotoImage

naph_img = Draw.MolToImage(naph) 
naph_photo = ImageTk.PhotoImage(naph_img)

cyhex_img = Draw.MolToImage(cyhex) 
cyhex_photo = ImageTk.PhotoImage(cyhex_img)

ace_img = Draw.MolToImage(ace) 
ace_photo = ImageTk.PhotoImage(ace_img)

root.mainloop()
