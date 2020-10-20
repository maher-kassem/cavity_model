#!/bin/bash
source activate cavity-model


# Settings
reduce_exe=reduce/reduce


# Parse ddG data PDBs
mkdir -p data/data_dms/pdbs_cleaned
mkdir -p data/data_dms/pdbs_parsed
mkdir -p data/data_guerois/pdbs_cleaned
mkdir -p data/data_guerois/pdbs_parsed


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
