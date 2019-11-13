import sys
import os
import numpy as np
import Bio
import Bio.PDB
import Bio.PDB.Vector
import simtk
import simtk.openmm
import simtk.openmm.app
import simtk.unit
import enum
import argparse
basepath = os.path.dirname(os.path.realpath(__file__))
utilspath = basepath + "/utils"
sys.path.insert(1, utilspath)
import grid

def extract_atomic_features(pdb_filename):
	"""Extract atomic features from pdb"""

	# Parse structure with biopython
	pdb_parser = Bio.PDB.PDBParser()
	pdb_id = os.path.basename(pdb_filename).split('.')[0]
	structure = pdb_parser.get_structure(pdb_id, pdb_filename)
	first_model = structure.get_list()[0]
	ppb = Bio.PDB.PPBuilder()
	
	sequence = [] # Only used for an assertion / sanity check 
	sequence_onehot = []
	resids_pdb = []
	chain_ids = []
	for i, chain in enumerate(first_model):
		# Append chainid
		chain_ids.append(chain.id)
	
		# Sequence of residue names for this chain
		sequence_chain = []
		for res in chain.get_residues():
			resids_pdb.append(res.id[1])
			sequence_chain.append(res.resname.strip())
	
		# Add to global container for this protein
		sequence.append(sequence_chain)
	
		# Convert residue names to amino acid indices
		aa_indices = []
		for aa in sequence_chain:
			try:
				aa_index = Bio.PDB.Polypeptide.three_to_index(aa)
			except:
				aa_index = 20
			aa_indices.append(aa_index)
	
		# Convert to one-hot encoding
		aa_onehot_chain = np.zeros((len(aa_indices), 21))
		aa_onehot_chain[np.arange(len(aa_indices)), aa_indices] = 1
		sequence_onehot.append(aa_onehot_chain)
	
	
	# Keep track of boundaries of individual chains
	chain_boundary_indices = np.cumsum([0]+[len(entry) for entry in sequence_onehot])
	
	# Collapse all chain segments into one. The individual chains
	# will be recoverable through chain_boundary_indices
	sequence_onehot = np.concatenate(sequence_onehot)
	
	# Convert resids of pdb to arr
	resids_pdb = np.array(resids_pdb)
	
	# Extract positions using OpenMM (I could have used Bio, but easier to adapt Wouters code with simtk and opens future possibility of using mass and charge)
	pdb_simtk = simtk.openmm.app.PDBFile(pdb_filename)
	positions = pdb_simtk.getPositions()
		
	# Save features in a dictionary
	features = {}
	features["atom_names"] = [] # Atom names
	features["res_indices"] = [] # Global res indices (no reset across chains and starts from 0)
	features["resids_pdb"] = [] # Resids in the actual odbs
	features["x"] = []
	features["y"] = []
	features["z"] = []
	
	# Iterate over chain,residue,atoms and extract features
	for i,chain in enumerate(pdb_simtk.getTopology().chains()):
		chain_start_index = chain_boundary_indices[i]
		for j,residue in enumerate(chain.residues()):
				for atom in residue.atoms():
	
					# Extract atom features
					index = atom.index
					position = list(positions[index].value_in_unit(simtk.unit.angstrom))
					features["atom_names"].append(atom.name)
					features["res_indices"].append(residue.index)
					features["x"].append(position[0])
					features["y"].append(position[1])
					features["z"].append(position[2])
					
					# Sanity check
					residue_index_local = residue.index - chain_start_index
					assert(residue.name == sequence[i][residue_index_local])
	
	# Convert valid lists to numpy arrays (even convert atom_names since its simpler to mask with despite being str)
	features["atom_names"] = np.array(features["atom_names"], dtype = 'a5')
	features["res_indices"] = np.array(features["res_indices"], dtype = np.int)
	features["x"] = np.array(features["x"], dtype = np.float32)
	features["y"] = np.array(features["y"], dtype = np.float32)
	features["z"] = np.array(features["z"], dtype = np.float32)

	return features, sequence_onehot, chain_ids, chain_boundary_indices, resids_pdb

def extract_coordinates(features, max_radius, include_center):
	"""Extract environment coordinates within specifies crieteria"""

	# Extract coordinates as normal numpy array
	position_array = np.vstack([features["x"], features["y"], features["z"]]).T

	# Retrieve residue indices as numpy int array
	# This array has many repeats, since it follows the sequence of atoms, not residues
	# It counts globally across chains, i.e. no resets and starts from zero.
	res_indices_glob = features['res_indices']
	res_indices_uniq = np.unique(res_indices_glob) 	# Has length == number of total residues

	# Begin
	selector_list = []
	indices_list = []
	xyz_ref_origo_list = []	
	for residue_index in res_indices_uniq:
		# Extract origin mask
		if (np.logical_and(res_indices_glob == residue_index, features['atom_names'] == b"N").any() and
			np.logical_and(res_indices_glob == residue_index, features['atom_names'] == b"CA").any() and
			np.logical_and(res_indices_glob == residue_index, features['atom_names'] == b"C").any()):

			N_mask = np.logical_and(res_indices_glob == residue_index, features['atom_names'] == b"N")
			CA_mask = np.logical_and(res_indices_glob == residue_index, features['atom_names'] == b"CA")
			C_mask = np.logical_and(res_indices_glob == residue_index, features['atom_names'] == b"C")
		else:
			# Store None to maintain indices
			indices_list.append(None)
			selector_list.append(None)
			continue
		
		# Extract origin
		pos_N = np.array( [features['x'][N_mask], features['y'][N_mask], features['z'][N_mask]] ).squeeze()
		pos_CA = np.array( [features['x'][CA_mask], features['y'][CA_mask], features['z'][CA_mask]] ).squeeze()
		pos_C = np.array( [features['x'][C_mask], features['y'][C_mask], features['z'][C_mask]] ).squeeze()

		# Define Cartesian coordinate system
		coordinate_system = grid.CoordinateSystem["cartesian"]
		z_direction = grid.ZDirection[grid.ZDirection.outward.name]

		# Define local coordinate system
		rot_matrix = grid.define_coordinate_system(pos_N, pos_CA, pos_C, z_direction)

		# Calculate coordinates relative to origin
		xyz = position_array - pos_CA

		# Rotate to the local reference
		xyz = np.dot(rot_matrix, xyz.T).T
		
		# Calculate radius
		r = np.sqrt(xyz[:,0]**2 + xyz[:,1]**2 + xyz[:,2]**2)

		if include_center:
			selector = np.where(r < max_radius)[0]
		else:
			# ALSO exclude features from residue itself
			selector = np.where(np.logical_and(r < max_radius, res_indices_glob != residue_index))[0]

		xyz_ref_origo_list.append(xyz[selector])
		selector_list.append(selector)

	# Find max number of atoms in an environment with radias = max_radius
	max_selector = max([len(selector) for selector in selector_list if selector is not None])
	selector_array = np.full((len(selector_list), max_selector), -1, dtype=np.int32)
	for i, selector in enumerate(selector_list):
		if selector is not None:
			selector_array[i,:len(selector)] = selector.astype(np.int32)

	# Output actual coordinates and atom types rather than embedding in a grid.
	atom_type_list = ["C", "N", "O", "H", "S", "P"]
	atom_types_numeric = np.array([atom_type_list.index(x.decode("utf-8")[0]) for x in features['atom_names']]) # Zero refers to the first letter of atom name 

	# Save selected coordinates in array with shape defined by the max number of atoms in any of environment
	xyz_ref_origo_arr = np.full(shape = [len(selector_list), max_selector, 3],
								fill_value = [-99, -99, -99],
								dtype = np.float32)
	for i, xyz_selected in enumerate(xyz_ref_origo_list):
		if xyz_selected is not None:
			xyz_ref_origo_arr[i, :xyz_selected.shape[0], :] = xyz_selected

	return xyz_ref_origo_arr, atom_types_numeric, selector_array


if __name__ == "__main__":
	# Argument parser
	parser = argparse.ArgumentParser()
	parser.add_argument("--pdb_in", type = str)
	parser.add_argument("--max_radius", type = float, default = 9.0)
	parser.add_argument("--include_center", type = str, default = "False")
	args_dict = vars(parser.parse_args())

	# Settings
	pdb_filename = args_dict["pdb_in"]
	pdb_id = os.path.basename(pdb_filename).split('.')[0]
	max_radius = args_dict["max_radius"]
	if args_dict["include_center"] == "False":
		include_center = False
	else:
		include_center = True
	

	# Extract atomic features and other relevant info
	features, sequence_onehot, chain_ids, chain_boundary_indices, resids_pdb = extract_atomic_features(pdb_filename)
	
	# Extract relevant coordinates (already masked with selector and referenced), all atom types and the mask for each residue (selector_array) 
	xyz_ref_origo_arr, atom_types_numeric, selector_array = extract_coordinates(features, max_radius, include_center)

	# Save as .npz
	np.savez_compressed("{:s}_coordinate_features".format(pdb_id),
						atom_types_numeric = atom_types_numeric,
						positions = xyz_ref_origo_arr,
						selector = selector_array,
						aa_onehot = sequence_onehot,
						chain_boundary_indices = chain_boundary_indices,
						chain_ids = chain_ids,
						residue_numbers = resids_pdb)
