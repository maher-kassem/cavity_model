# Protein Gap Model Demonstration 
This is  demo of a protein structure based amino acid classifier using 3D convolutions. The central idea is to look at the environment of a query amino acid, remove the atoms the query amino acid and have the classifier attempt to predict the label of the missing amino acid given the environment defined as by the sphere around the C-alpha atom of the query amino acid. 

## Foundation for Transfer Learning
Given a well-trained classifier, the accuracy on a homology reduced test should be around 60 %. In this [paper](http://papers.nips.cc/paper/6935-spherical-convolutions-and-their-application-in-molecular-modelling), they reason that the difference in the log-likelihood of a wildtype label and a mutant label must have some correlation with the change in stability upon that mutation which turned out to be true. Furthermore, they found that they could add a down stream model on top of the classifier's output log likelihoods and gain state of the art performance on the prediction of protein stability change upon mutation.

## Documentation

### Step by step guide to get this running yourself (on a unix machine with miniconda installed)
0. Clone the repository and change directory  
`git clone https://github.com/mahermkassem/Protein_Gap_Model_Demo.git`  
`cd Protein_Gap_Model_Demo`

1. Install the conda environment in the `conda_env` directory.

`conda env create -f conda_env/py3.6-gapmodeldemo.yml`  
`conda activate py3.6-gapmodeldemo`

2. Install reduce. This program is used by my parser to add missing hydrogens to the proteins  
`git clone https://github.com/rlabduke/reduce.git`  
`cd reduce`

`make; make install` # This should give an error but provide the reduce executable in this directory

3. Download and Parse the data using my bash script. It downloads the PDBs listed in `data/transfer_learning_data/pdbids_250.txt` and parses them. The final parsed files are saved structures numpy arrays in `data/transfer_learning_data/structural_environments/`.

`cd ../data`

`./download_and_process_data.sh`

Once it has downloaded, cleaned and parsed all the PDBs, you should be able to run all the code in the `prot_gap_model.ipynb` notebook.
