#!/bin/bash
conda_activate_path=$1
source $conda_activate_path cavity-model


# Settings
reduce_exe=reduce/reduce

# Parse ddG data PDBs
mkdir -p data/data_dms/pdbs_cleaned
mkdir -p data/data_dms/pdbs_parsed
mkdir -p data/data_guerois/pdbs_cleaned
mkdir -p data/data_guerois/pdbs_parsed
mkdir -p data/data_protein_g/pdbs_cleaned
mkdir -p data/data_protein_g/pdbs_parsed
mkdir -p data/data_symmetric/pdbs_cleaned
mkdir -p data/data_symmetric/pdbs_parsed

## DMS DATA
# Clean pdbs in data/dat_dms/pdbs_raw/*.pdb
counter=1
dms_pdbs=data/data_dms/pdbs_raw/*.pdb
n_pdbs=$(echo $dms_pdbs | wc -w)
for pdb in $dms_pdbs;
do
    python pdb_parser_scripts/clean_pdb.py --pdb_file_in $pdb  \
                                           --out_dir data/data_dms/pdbs_cleaned/ \
                                           --reduce_exe $reduce_exe &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    echo "Successfully cleaned $pdb and added it to data/data_dms/pdbs_cleaned/. $counter/$n_pdbs."
    else
    echo "Error when cleaning $pdb. Skipping.." >&2
    fi
    counter=$((counter+1))
done

# Parse pdbs and save in npz format
counter=1
out_dir=data/data_dms/pdbs_parsed
for pdb_clean in data/data_dms/pdbs_cleaned/*.pdb
do
    python pdb_parser_scripts/extract_environments.py --pdb_in $pdb_clean --out_dir $out_dir  &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    base_pdb=$(basename $pdb_clean)
    echo "Successfully parsed $pdb_clean and moved to $out_dir. \
Finished $counter/$n_pdbs."
    else
    echo "Error extracting $pdb_clean. Skipping.." >&2
    fi
    counter=$((counter+1))
done

## GUEROIS DATA
# Clean pdbs in data/dat_guerois/pdbs_raw/*.pdb
counter=1
guerois_pdbs=data/data_guerois/pdbs_raw/*.pdb
n_pdbs=$(echo $guerois_pdbs | wc -w)
for pdb in $guerois_pdbs;
do
    python pdb_parser_scripts/clean_pdb.py --pdb_file_in $pdb  \
                                           --out_dir data/data_guerois/pdbs_cleaned/ \
                                           --reduce_exe $reduce_exe &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    echo "Successfully cleaned $pdb and added it to data/data_guerois/pdbs_cleaned/. $counter/$n_pdbs."
    else
    echo "Error when cleaning $pdb. Skipping.." >&2
    fi
    counter=$((counter+1))
done

# Parse pdbs and save in npz format
counter=1
out_dir=data/data_guerois/pdbs_parsed
for pdb_clean in data/data_guerois/pdbs_cleaned/*.pdb
do
    python pdb_parser_scripts/extract_environments.py --pdb_in $pdb_clean --out_dir $out_dir  &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    base_pdb=$(basename $pdb_clean)
    echo "Successfully parsed $pdb_clean and moved to $out_dir. \
Finished $counter/$n_pdbs."
    else
    echo "Error extracting $pdb_clean. Skipping.." >&2
    fi
    counter=$((counter+1))
done

## PROTEIN G DATA
# Clean pdbs in data/dat_protein_g/pdbs_raw/*.pdb
counter=1
protein_g_pdbs=data/data_protein_g/pdbs_raw/*.pdb
n_pdbs=$(echo $protein_g_pdbs | wc -w)
for pdb in $protein_g_pdbs;
do
    python pdb_parser_scripts/clean_pdb.py --pdb_file_in $pdb  \
                                           --out_dir data/data_protein_g/pdbs_cleaned/ \
                                           --reduce_exe $reduce_exe &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    echo "Successfully cleaned $pdb and added it to data/data_protein_g/pdbs_cleaned/. $counter/$n_pdbs."
    else
    echo "Error when cleaning $pdb. Skipping.." >&2
    fi
    counter=$((counter+1))
done

# Parse pdbs and save in npz format
counter=1
out_dir=data/data_protein_g/pdbs_parsed
for pdb_clean in data/data_protein_g/pdbs_cleaned/*.pdb
do
    python pdb_parser_scripts/extract_environments.py --pdb_in $pdb_clean --out_dir $out_dir  &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    base_pdb=$(basename $pdb_clean)
    echo "Successfully parsed $pdb_clean and moved to $out_dir. \
Finished $counter/$n_pdbs."
    else
    echo "Error extracting $pdb_clean. Skipping.." >&2
    fi
    counter=$((counter+1))
done


## SYMMETRIC DATA
# Download the protein structures (.pdb files) from the protein data bank
counter=1
n_pdbs=$(cat data/data_symmetric/pdbids.txt | wc -l)
while read p;
do
    # Download
    wget -O data/data_symmetric/pdbs_raw/$p.pdb \
    'http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId='$p\
    &> /dev/null

    # Check that wget gives exit code 0. If not, remove potentially non-existant or partial PDB file.
    if [ $? -eq 0 ]
    then
        echo "Successfully downloaded $p.pdb to data/data_symmetric/pdbs_raw/$p.pdb. $counter/$n_pdbs."
    else
        echo "Could not download $p.pdb. Attempting to remove potentially non-existant or \
potentially partial $p.pdb" >&2
    
        rm data/pdbs/raw/$p.pdb 
    fi
    counter=$((counter+1))
done < data/data_symmetric/pdbids.txt


# Clean pdbs in data/data_symmetric/pdbs_raw/*.pdb
counter=1
symmetric_pdbs=data/data_symmetric/pdbs_raw/*.pdb
n_pdbs=$(echo $symmetric_pdbs | wc -w)
for pdb in $symmetric_pdbs;
do
    python pdb_parser_scripts/clean_pdb.py --pdb_file_in $pdb  \
                                           --out_dir data/data_symmetric/pdbs_cleaned/ \
                                           --reduce_exe $reduce_exe &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    echo "Successfully cleaned $pdb and added it to data/data_symmetric/pdbs_cleaned/. $counter/$n_pdbs."
    else
    echo "Error when cleaning $pdb. Skipping.." >&2
    fi
    counter=$((counter+1))
done

# Parse pdbs and save in npz format
counter=1
out_dir=data/data_symmetric/pdbs_parsed
for pdb_clean in data/data_symmetric/pdbs_cleaned/*.pdb
do
    python pdb_parser_scripts/extract_environments.py --pdb_in $pdb_clean --out_dir $out_dir  &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    base_pdb=$(basename $pdb_clean)
    echo "Successfully parsed $pdb_clean and moved to $out_dir. \
Finished $counter/$n_pdbs."
    else
    echo "Error extracting $pdb_clean. Skipping.." >&2
    fi
    counter=$((counter+1))
done
