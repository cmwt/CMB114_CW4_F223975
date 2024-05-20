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

In order for the code to run, the user must first download the most up-to-date version of Python and Git on their computer and then install the correct libraries. 

Python can be downloaded from the following link:
https://www.python.org/downloads/

Git can be downloaded here:
https://git-scm.com/downloads

When installing Python it is important to check the PATH tickbox when prompted since this allows the libraries to be downloaded to the operating system and be accessible to Python. 

In order to install the libraries required, type the following into the command prompt (or terminal on MacOS):

`pip install rdkit`

RDKit is the library that allows the application to process molecular data.

PIL is another library used in the code, that should come standard when downloading Python, that allows us to view the molecules as images and have them drawn to a tkinter canvas. 

tkinter is used for the GUI itself, however, this is a standard library of Python and, therefore, does not need to be specifically installed. 

To run the code itself, the user should open the command line and navigate to the directory within which they would like to clone the repository. Then, type the following:

`git clone `

After this has run (the user may be prompted to sign into git on their device) navigate into the 'CMB114_CW4_F223975' directory and run:

`python main.py` 

This should allow the user to use the software freely and unlock all of its functionality. 

For further help with running Orca calculations, please visit the Orca website at the following link:

https://www.orcasoftware.de/tutorials_orca/index.html

# Purpose of the Application:

This GUI-based application takes inspiration from the Avogadro software in that it can be used to view molecules and the properties of these pre-programmed molecules can be seen on screen when selected. As well as this feature, the user can write in their own custom molecules, as a SMILES string, for them to be drawn. Using the features of the GUI, the user can then generate an Orca input file that can be run on the command line, providing the Orca software is downloaded, to run molecular dynamics calculations on these molecules. 

The code for the application is written in such a way that allows for the addition of many more molecules at a later date, allowing for databases of molecules to be read in with their SMILES code and properties being displayed within the GUI window for the user to see. 

# Tutorials:

For tutorials on the features of the software, and testing surrounding its more complex features, please refer to 'tutorials_and_testing.docx' in the 'docs' directory. 
