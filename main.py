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
from tkinter import messagebox
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


# Function to draw the molecule from SMILES
def draw_custom():
    smiles = custom_entry.get()  # Get the SMILES string from the entry widget
    try:
        mol = Chem.MolFromSmiles(smiles)  # Convert the SMILES string to a molecule object
        if mol is None:  # Checks if input is a valid SMILES code
            raise ValueError("Invalid SMILES string")
        
        img = Draw.MolToImage(mol)  # Draw the molecule and get a PIL image
        tk_img = ImageTk.PhotoImage(img)
        
        # Clears previous image from the canvas
        canvas.delete("all")
        canvas.create_image(150, 150, image=tk_img)  # Display the image on the canvas
        canvas.image = tk_img  # Keep a reference to avoid garbage collection
        
    except Exception as e:  # Catch any exceptions that occur
        messagebox.showerror("Error", str(e))  # Show an error message box
        

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

frame = ttk.Frame(root, width = 300, height = 1000, borderwidth = 5, relief = tk.GROOVE)
frame.grid_propagate(False)
frame.grid(row=1, column=1, sticky="NW")

drop_label = tk.Label(frame, text = "Pre-Programmed Molecules:")
drop_label.grid(row=1, column=1, sticky="W", pady=5)
drop = ttk.Combobox(frame, value= molecules)
drop.grid(row=2, column=1, sticky="W", pady=5)
drop.bind("<<ComboboxSelected>>", option_selected)

clear_but = tk.Button(frame, text="Clear", command=clear_frame)
clear_but.grid(row=2, column=2, sticky="W", pady=5, padx=5)

custom_label = tk.Label(frame, text="Custom SMILES code:")
custom_label.grid(row=3, column=1, sticky="W", pady=5)

custom_entry = tk.Entry(frame)
custom_entry.grid(row=4, column=1, sticky="W", pady=5)

draw_button = tk.Button(frame, text="Draw", command=lambda:[draw_custom()])
draw_button.grid(row=4, column=2, sticky="W", pady=5, padx=5)

atom_label = tk.Label(frame, text= "Atom Numbers:")
atom_label.grid(row=5, column=1, sticky="W", pady=5)
atom_num = tk.Button(frame, text = "On", command=atom_switch)
atom_num.grid(row=5, column=2, sticky="W", pady=5, padx=5)

prop_label = tk.Label(frame, text="Molecule Properties:")
prop_label.grid(row=6, column=1, sticky="W", pady=5)

reference_label = tk.Label(frame)
reference_label.grid(row=7, column=1, sticky="W", pady=5)


text = {"Benzene":"Benzene:\nMolar Mass: 78.11 g/mol\nBoiling Point: 80.1°C\nMelting Point: 5.5°C",
        "Naphthalene":"Naphthalene:\nMolar Mass: 128.17 g/mol\nBoiling Point: 218°C\nMelting Point: 80.3°C",
        "Cyclohexane":"Cyclohexane:\nMolar Mass: 84.16 g/mol\nBoiling Point: 80.8°C\nMelting Point: 6.5°C",
        "Acetone":"Acetone:\nMolar Mass: 58.1 g/mol\nBoiling Point: 56°C\nMelting Point: -95°C"
                  }

# *** Right Frame ***

frame1 = ttk.Frame(root, width = 1200, height = 1000, borderwidth = 5, relief = tk.GROOVE)
frame1.grid_propagate(False)
frame1.grid(row=1, column=2, sticky="NW")

canvas = tk.Canvas(frame1, width=300, height=300)
canvas.grid(row=2, column=2)

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
