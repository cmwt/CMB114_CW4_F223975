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
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolTransforms
import os
from lib.mol import *
from lib.orca import *

# ***************************
# *                         *
# *     Main Functions      *
# *                         *
# ***************************

# Function to print properties text in left frame
def reference_text(x):
    reference_label.config(text = str(text[x]) )

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
    reference_label.configure(text="")

# *****************
# *               *
# *     Orca      *
# *               *
# *****************

def orca_input_file(basis_set, calc_type, coordinates):

    filename = custom_entry.get() # Sets the filename to whatever the SMILES string is 
    
    input_dir = "ORCA_input"
    if not os.path.exists(input_dir): # If ORCA directory does not exist yet, create it
        os.makedirs(input_dir)

    filepath = os.path.join(input_dir, filename) # Allows a file to be opened within the new directory

    with open(f"{filepath}.txt", "w") as f:
        f.write("# CW4 - Molecular Dynamics generated ORCA input file\n")
        f.write(f"\n! {calc_type} {basis_set}\n")
        f.write("* xyz 0 1\n")
        f.write(coordinates)
        f.write("*\n")

def generate_orca():
    basis_set = basis_set_drop.get()
    calc_type = calc_type_drop.get()
    smiles = custom_entry.get()
    
    # Generate cartesian coordinates from SMILES
    mol2 = Chem.MolFromSmiles(smiles)
    if mol2 is None:
        status_label.config(text="Invalid SMILES code")
        return
    
    mol2 = Chem.AddHs(mol2)
    AllChem.EmbedMolecule(mol2, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol2)
    
    coordinates = ""
    conf = mol2.GetConformer()
    for atom in mol2.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        coordinates += f"  {atom.GetSymbol()}  {pos.x:.6f}  {pos.y:.6f}  {pos.z:.6f}\n"

    orca_input_file(basis_set, calc_type, coordinates)
    status_label.config(text="Input file created in ORCA_input")

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

frame = ttk.Frame(root, width = 400, height = 1000, borderwidth = 5, relief = tk.GROOVE)
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

prop_label = tk.Label(frame, text="Molecule Properties:")
prop_label.grid(row=6, column=1, sticky="W", pady=5)

reference_label = tk.Label(frame)
reference_label.grid(row=7, column=1, sticky="W", pady=5)

basis_label = tk.Label(frame, text="Select Basis Set:")
basis_label.grid(row=8, column=1, pady=5)
basis_set_drop = ttk.Combobox(frame, textvariable=basis_sets, values=basis_sets)
basis_set_drop.grid(row=8, column=2, padx=10, pady=5)

calc_label = tk.Label(frame, text="Select Calculation Type")
calc_label.grid(row=9, column=1, pady=5)
calc_type_drop = ttk.Combobox(frame, textvariable=calc_types, values=calc_types)
calc_type_drop.grid(row=9, column=2, padx=10, pady=5)

gen_but = tk.Button(frame, text="Generate Orca Input", command=generate_orca)
gen_but.grid(row=10, column=0, columnspan=2, pady=5)

status_label = tk.Label(frame, text="")
status_label.grid(row=11, column=0, columnspan=2, pady=5)


text = {"Benzene":"Benzene:\nMolar Mass: 78.11 g/mol\nBoiling Point: 80.1°C\nMelting Point: 5.5°C",
        "Naphthalene":"Naphthalene:\nMolar Mass: 128.17 g/mol\nBoiling Point: 218°C\nMelting Point: 80.3°C",
        "Cyclohexane":"Cyclohexane:\nMolar Mass: 84.16 g/mol\nBoiling Point: 80.8°C\nMelting Point: 6.5°C",
        "Acetone":"Acetone:\nMolar Mass: 58.1 g/mol\nBoiling Point: 56°C\nMelting Point: -95°C"
                  }

# *** Right Frame ***

frame1 = ttk.Frame(root, width = 2000, height = 1000, borderwidth = 5, relief = tk.GROOVE)
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
