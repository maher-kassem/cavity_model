# Reproducing Cavity Model results
The purpose of this repository is to reproduce the core results in PAPER entitled TITLE.

## Data
Protein G  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; explain  
DMS  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; explain  
Guerois  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; explain  

## Code
dir x  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; explain  
file y  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; explain  

## Installation

### Step by step guide to get this running yourself (**Only tested on linux with Minoconda**)
0. Clone the repository and change directory  
`git clone https://github.com/mahermkassem/Cavity_Model_Demo.git`  
`cd Cavity_Model_Demp/`

1. Install and activate the provided exported conda environment.  
`conda create --name cavity-model python=3.6`  
`conda activate cavity-model`  
`conda install notebook black nb_black pandas scipy numpy=1.17.3 pdbfixer=1.5 pytorch=1.2.0 biopython=1.72 openmm=7.3.1 matplotlib=3.1.1 -c omnia -c conda-forge -c anaconda -c defaults`  

2. Install reduce. This program is used by my parser to add missing hydrogens to the proteins  
`git clone https://github.com/rlabduke/reduce.git`  
`cd reduce/`  
`make; make install` # This might give an error but provide the reduce executable in this directory

3. You should be able to run all the code in the `cavity_model.ipynb` notebook.
