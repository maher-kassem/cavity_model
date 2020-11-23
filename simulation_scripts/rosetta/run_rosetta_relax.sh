#!/bin/bash

# Path to cleaned pdb and rosetta relax binary
pdb=$1
rosetta_relax_path=/home/chz526/software/rosetta_bin_linux_2018.33.60351_bundle/main/source/bin/relax.static.linuxgccrelease

# Run
$rosetta_relax_path \
  -database \
  /home/chz526/software/rosetta_bin_linux_2018.33.60351_bundle/main/database \
  -s \
  $pdb \
  -ex1 \
  -ex2 \
  -packing:linmem_ig 10 \
  -no_optH false \
  -flip_HNQ \
  -relax:dualspace true \
  -nonideal true \
  -missing_density_to_jump \
  -use_input_sc \
  -ignore_unrecognized_res \
  -beta \
  -score:weights beta_nov16_cart \
  -relax:constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
  -relax:min_type lbfgs_armijo_nonmonotone \
  -relax:ramp_constraints false \
  -ignore_zero_occupancy false 
