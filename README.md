# Protein Cavity Model Demonstration 
This is  demo of a protein structure based amino acid classifier using 3D convolutions. The central idea is to look at the environment of a query amino acid, remove the atoms the query amino acid and have the classifier attempt to predict the label of the missing amino acid given the environment defined as by the sphere around the C-alpha atom of the query amino acid. 

## Predicting change in protein stability
Given a well-trained classifier, the accuracy on a homology reduced test should be around 60 %. In this [paper](http://papers.nips.cc/paper/6935-spherical-convolutions-and-their-application-in-molecular-modelling), they reason that the difference in the log-likelihood of a wildtype label and a mutant label must have some correlation with the change in stability upon that mutation which turned out to be true. Furthermore, they found that they could add a down stream model on top of the classifier's output log likelihoods and gain state of the art performance on the prediction of protein stability change upon mutation.

## Documentation

### Step by step guide to get this running yourself (**on a linux machine with Miniconda installed**)
0. Clone the repository and change directory  
`git clone https://github.com/mahermkassem/Cavity_Model_Demo.git`  
`cd Cavity_Model_Demp/`

1. Install and activate the provided exported conda environment.  
`conda env create -f environment.yaml`  
`conda activate cavity-model`

2. Install reduce. This program is used by my parser to add missing hydrogens to the proteins  
`git clone https://github.com/rlabduke/reduce.git`  
`cd reduce/`  
`make; make install` # This should give an error but provide the reduce executable in this directory

3. You should be able to run all the code in the `cavity_model.ipynb` notebook.
