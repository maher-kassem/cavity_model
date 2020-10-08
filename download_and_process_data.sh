#!/bin/bash
source activate cavity-model

# Settings
pdbid_filename=$1
reduce_exe=reduce/reduce
n_pdbs=$(cat $pdbid_filename | wc -l)

# Create data directories
mkdir -p data/raw
mkdir -p data/cleaned
mkdir -p data/parsed

# Download the protein structures (.pdb files) from the protein data bank
counter=1
while read p;
do
    # Download
    wget -O data/raw/$p.pdb 'http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId='$p &> /dev/null

    # Check that wget gives exit code 0. If not, remove potentially non-existant or partial PDB file.
    if [ $? -eq 0 ]
    then
	echo "Successfully downloaded $p.pdb to data/raw/$p.pdb. $counter/$n_pdbs."
    else
	echo "Could not download $p.pdb. Attempting to remove potentially non-existant or potentially partial $p.pdb" >&2
	rm data/raw/$p.pdb 
    fi
    counter=$((counter+1))
done < $pdbid_filename

# Clean all successfully downloaded pdbs
counter=1
for pdb in data/raw/*.pdb;
do
    python scripts/clean_pdb.py --pdb_file_in $pdb  --out_dir data/cleaned/ --reduce_exe $reduce_exe &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
	echo "Successfully cleaned $pdb to data/cleaned/$pdb.pdb. $counter/$n_pdbs."
    else
	echo "Error when cleaning $pdb. Skipping.." >&2
    fi
    counter=$((counter+1))
done

# Parse pdbs and save in npz format
cd data/parsed/
counter=1
for pdb_clean in ../cleaned/*pdb
do
    python ../../scripts/extract_environments.py --pdb_in $pdb_clean &> /dev/null

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
