#!/bin/bash

pdb_relaxed=$1
mutfile=$2
cartesian_ddg_path=/home/chz526/software/rosetta_bin_linux_2018.33.60351_bundle/main/source/bin/cartesian_ddg.static.linuxgccrelease

$cartesian_ddg_path \
  -database /home/chz526/software/rosetta_bin_linux_2018.33.60351_bundle/main/database \
  -fa_max_dis 9.0 \
  -ddg::dump_pdbs true \
  -ddg:iterations 3 \
  -score:weights beta_nov16_cart \
  -missing_density_to_jump \
  -ddg:bbnbr 1 \
  -beta_cart \
  -ddg:mut_only \
  -ex1 \
  -ex2 \
  -s $pdb_relaxed \
  -ddg:mut_file $mutfile \
  -ddg::cartesian

