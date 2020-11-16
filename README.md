# Reproducing Cavity Model results
The purpose of this repository is to reproduce the core results in PAPER entitled TITLE.

## Data
Protein G  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `data/data_protein_g/`(https://doi.org/10.1101/484949)  
DMS  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `data/data_dms/`(https://dx.doi.org/10.1038%2Fs41588-018-0122-z)  
Guerois  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `data/data_guerois/`(https://doi.org/10.1016/s0022-2836(02)00442-4)  
Symmetric  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `data/data_symmetric/`(https://doi.org/10.1093/bioinformatics/bty348)  

## Code
`cavity_model.py` contains the main cavity model and downstream model classes and data loaders  
`cavity_model_pipeline.py` contains the full pipeline to reproduce the results of the paper  
`helpers.py` protected helper functions for the pipeline  
`pdb_parser_scripts/`contains a pdb cleaning script and a parser script to extract residue environments  
`simulation_script/` contains scripts for running explicit solvent MD and Rosetta  
`get_and_parse_pdbs_for_cavity_model.sh` bash script to parse PDBs for the cavity model  
`get_and_parse_pdbs_for_downstream_model.sh` Bash script to parse PDBs for the downstream model (i.e. PDBs of the ddG datasets)   

## Installation
### (**Only tested on linux with Minoconda**)
0. Clone the repository and change directory  
`git clone https://github.com/mahermkassem/Cavity_Model_Demo.git`  
`cd Cavity_Model_Demp/`

1. Install and activate the provided exported conda environment.  
`conda create --name cavity-model python=3.6`  
`conda activate cavity-model`  
`conda install notebook black nb_black pandas scipy numpy=1.17.3 pdbfixer=1.5 pytorch=1.2.0 biopython=1.72 openmm=7.3.1 matplotlib=3.1.1 openmm=7.3.1 -c omnia -c conda-forge -c anaconda -c defaults`  

2. Install reduce. This program is used by my parser to add missing hydrogens to the proteins  
`git clone https://github.com/rlabduke/reduce.git`  
`cd reduce/`  
`make; make install` # This might give an error but provide the reduce executable in this directory

3. You should be able to run all the code in the `cavity_model.ipynb` notebook.  
`cd ../`  
`jupyter notebook`
