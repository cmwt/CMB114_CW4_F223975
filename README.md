# CMB114_CW4_F223975

This repository houses all of the files and directories for Dr. Kenny Jolley's CW4 project for CMB114,
in which a software application surrounding molecular dynamics must be produced.

The requirements for the project are set out below:

- Must have a very good graphical user interface (GUI)
- Organise source files in an appropriate directory structure
- Create a README file (this file) that should contain
    - How to build the code
    - Purpose of the application
    - Tutorials for any features of the application
- Must be strongly related to molecular dynamics
- Use GitHub to develop the code, showing commit messages etc.

# How to build the code:

In order for the code to run, the user must first download the most up-to-date version of Python on their computer and then install the correct libraries. 

Python can be downloaded from the following link:
https://www.python.org/downloads/

When installing Python it is important to check the PATH tickbox when prompted since this allows the libraries to be downloaded to the operating system and be accessible to Python. 

In order to install the libraries required, type the following into the command prompt (or terminal on MacOS):

'pip install rdkit'

'pip install PIL'

RDKit is the library that allows the application to process molecular data.

PIL is the library that allows us to view the molecules as images and have them drawn to a tkinter canvas. 

tkinter is used for the GUI itself, however, this is a standard library of Python and, therefore, does not need to be specifically installed. 

# Purpose of the Application:

This GUI-based application takes inspiration from the Avogadro software in that it can be used to view molecules and properties of these molecules can be seen on screen when selected. 
The code for the application is written in such a way that allows for the addition of many more molecules at a later date, allowing for databases of molecules to be read in with their SMILES code and properties being displayed within the GUI window for the user to see. 

# Tutorials:
