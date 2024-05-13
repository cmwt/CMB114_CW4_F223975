# This file contains the data for each of the molecules that can be displayed
# on the tkinter window.
# All of the functions related to the molecules can also be found in here. For
# example, the draw_xxx() functions for each molecule.

from PIL import Image, ImageTk
from rdkit import Chem
from rdkit.Chem import Draw

# ***** Molecules *****

# Benzene
benzene_SMILES = "c1ccccc1" # Assigns benzene a SIMLES moleucle format
benzene = Chem.MolFromSmiles(benzene_SMILES) # Creates the SMILES molecule
# Naphthalene
naph_SMILES = "c1c2ccccc2ccc1"
naph = Chem.MolFromSmiles(naph_SMILES) 
# Cyclohexane
cyhex_SMILES = "C1CCCCC1"
cyhex = Chem.MolFromSmiles(cyhex_SMILES) 
# Acetone
ace_SMILES = "CC(=O)C"
ace = Chem.MolFromSmiles(ace_SMILES)
