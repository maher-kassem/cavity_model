# Protein Cavity Model Demonstration 
The purpose of this repository is to reproduce the core results in PAPER entitled TITLE.

## Data
Protein G  
    explain  
DMS  
    explain  
Guerois  
    explain  

## Code
dir x  
    explain  
file y  
    explain  

## Installation

### Step by step guide to get this running yourself (**Only tested on linux with Minoconda**)
0. Clone the repository and change directory  
`git clone https://github.com/mahermkassem/Cavity_Model_Demo.git`  
`cd Cavity_Model_Demp/`

1. Install and activate the provided exported conda environment.  
`conda env create -f environment.yaml`  
`conda activate cavity-model`

2. Install reduce. This program is used by my parser to add missing hydrogens to the proteins  
`git clone https://github.com/rlabduke/reduce.git`  
`cd reduce/`  
`make; make install` # This might give an error but provide the reduce executable in this directory

3. You should be able to run all the code in the `cavity_model.ipynb` notebook.
