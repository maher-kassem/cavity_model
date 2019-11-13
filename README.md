# Protein Gap Model Demonstration 
This is  demo of a protein structure based amino acid classifier using 3D convolutions. The central idea is to look at the environment of a query amino acid, remove the atoms the query amino acid and have the classifier attempt to predict the label of the missing amino acid given the environment defined as by the sphere around the C-alpha atom of the query amino acid. 

## Foundation for Transfer Learning
Given a well-trained classifier, the accuracy on a homology reduced test should be around 60 %. In this [paper](http://papers.nips.cc/paper/6935-spherical-convolutions-and-their-application-in-molecular-modelling), they reason that the difference in the log-likelihood of a wildtype label and a mutant label must have some correlation with the change in stability upon that mutation which turned out to be true. Furthermore, they found that they could add a down stream model on top of the classifier's output log likelihoods and gain state of the art performance on the prediction of protein stability change upon mutation.

## Documentation

### Step by step guid to get this running yourself
1. Install the conda environment in the `conda_env` directory.
`conda env create -f conda_env/py3.6-gapmodeldemo.yml`
2

### Parsing the data
Parsing PDB structures can be tedious given the massive heterogeneity between PDBs despite being in the same format. 

### The model
The prot_gap_model.ipynb notebook is a step by step demonstration of how to train an amino acid classifier on a small set (for convenience) of PDB protein structures (N=250).
