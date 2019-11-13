#!/bin/bash
source activate py3.6-ddg


# Settings
pdbid_filename=transfer_learning_data/pdbids_250.txt
reduce_exe=../reduce/reduce
n_pdbs=$(cat $pdbid_filename | wc -l)

# Download the protein structures (.pdb files) from the protein data bank
counter=1
while read p;
do
    # Download
    wget -O transfer_learning_data/raw/$p.pdb 'http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId='$p &> /dev/null

    # Check that wget gives exit code 0. If not, remove potentially non-existant or partial PDB file.
    if [ $? -eq 0 ]
    then
	echo "Successfully downloaded $p.pdb. Finished $counter/$n_pdbs."
    else
	echo "Could not download $p.pdb. Attempting to remove potentially non-existant or potentially partial $p.pdb" >&2
	rm $pdb.pdb # Might give an error, depending on why wget failed.
    fi
    counter=$((counter+1))
done < $pdbid_filename

# Clean all successfully downloaded pdbs
counter=1
for pdb in transfer_learning_data/raw/*.pdb;
do
    python ../scripts/clean_pdb.py --pdb_file_in $pdb  --out_dir transfer_learning_data/cleaned/ --reduce_exe $reduce_exe &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
	echo "Successfully cleaned $pdb. Finished $counter/$n_pdbs."
    else
	echo "Error when cleaning $pdb. Skipping.." >&2
    fi
    counter=$((counter+1))
done

# Extract structural environments from cleaned pdb files
cd transfer_learning_data/structural_environments/
counter=1
for pdb_clean in ../cleaned/*pdb
do
    python ../../../scripts/extract_environments.py --pdb_in $pdb_clean &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
	base_pdb=$(basename $pdb_clean)
	echo "Successfully extracted environments from $base_pdb. Finished $counter/$n_pdbs."
    else
	echo "Error extracting $pdb_clean. Skipping.." >&2
    fi
    counter=$((counter+1))
done
