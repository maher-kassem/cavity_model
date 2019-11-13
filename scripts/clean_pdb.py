from io import StringIO, BytesIO
import tempfile
import subprocess
import os
import sys
import Bio.PDB
import Bio.PDB.Polypeptide
import Bio.SeqIO
import pdbfixer
import simtk
import simtk.openmm
import simtk.openmm.app
import argparse

class NonHetSelector(Bio.PDB.Select):
		""" Remove HET atoms and choose first conformation of disordered atoms"""
		def accept_residue(self, residue):
			norm_res_bool = residue.get_resname() in [pdbfixer.pdbfixer.substitutions[key] for key in pdbfixer.pdbfixer.substitutions]
			abnorm_res_bool = residue.get_resname() in [key for key in pdbfixer.pdbfixer.substitutions]
			return norm_res_bool or abnorm_res_bool # Accept abnorm since they are converted later
		def accept_atom(self, atom):
			return (not atom.is_disordered() or atom.get_altloc() == 'A' or atom.get_altloc() == '1') and atom.id[0] in ['C', 'H', 'N', 'O', 'S', 'P']
class PDBFixerResIdentifiabilityIssue(Exception):
	pass


if __name__ == "__main__":
	# Argument Parser
	parser = argparse.ArgumentParser()
	parser.add_argument("--pdb_file_in", type = str)
	parser.add_argument("--out_dir", type = str)
	parser.add_argument("--reduce_exe", type = str)

	# Parse arguments
	args_dict = vars(parser.parse_args())
	pdb_input_filename = args_dict["pdb_file_in"]
	out_dir = args_dict["out_dir"]
	reduce_executable = args_dict["reduce_exe"]
	
	# Test that all args are provided since we don't have default values
	if pdb_input_filename is None:
		sys.exit("Error, <pdb_file_in> not provided. Exitting...")
	if reduce_executable is None:
		sys.exit("Error, <reduce_exec> not provided. Exitting...")
	if out_dir is None:
		sys.exit("Error, <out_dir> not provided. Exitting...")
	
	
	with tempfile.NamedTemporaryFile(mode = "wt", delete=True) as temp1:
		# Add hydrogens using reduce program
		command = [reduce_executable, '-BUILD',
				   '-DB', os.path.join(os.path.dirname(reduce_executable), 'reduce_wwPDB_het_dict.txt'),
				   '-Quiet',
				   pdb_input_filename]
		error_code = subprocess.Popen(command, stdout=temp1).wait()
		temp1.flush()
	
		# Parse with Biopython to remove HET and disordered atoms
		pdbid = pdb_input_filename.split("/")[-1].split(".pdb")[0]
		pdb_parser = Bio.PDB.PDBParser(PERMISSIVE=0)
		structure = pdb_parser.get_structure(pdbid, temp1.name)
		
		first_model = structure[0]
		pdbio = Bio.PDB.PDBIO()
	
		# Save temp pdb with NonHetSelector filter
		with tempfile.NamedTemporaryFile(mode = "wt", delete=True) as temp2:
			pdbio.set_structure(first_model)
			pdbio.save(temp2, select = NonHetSelector())
			temp2.flush()
	
			# Load temp pdb and replace altloc chars to " "
			with tempfile.NamedTemporaryFile(mode = "wt", delete=True) as temp3:
				structure = pdb_parser.get_structure(temp2.name, temp2.name)
				first_model = structure[0]
				for chain in first_model:
					for res in chain:
						for atom in res:
							atom.set_altloc(' ')
				pdbio.set_structure(first_model)
				pdbio.save(temp3)
				temp3.flush()
	
				# Use PDBFixer to fix common PDB errors
				fixer = pdbfixer.PDBFixer(temp3.name)		
				fixer.findMissingResidues()
				fixer.findNonstandardResidues()
				fixer.replaceNonstandardResidues()
				fixer.findMissingAtoms()
				fixer.addMissingAtoms()
				fixer.addMissingHydrogens(7.0)
				with tempfile.NamedTemporaryFile(mode = "wt", delete=True) as temp4:
					simtk.openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, temp4, keepIds=False)
					temp4.flush()
					# Fix IDs manually since pdbfixer does not preserve insertion codes
					pdb_parser = Bio.PDB.PDBParser(PERMISSIVE=1)
					structure_before = pdb_parser.get_structure(temp3.name, temp3.name)
					structure_after = pdb_parser.get_structure(temp4.name, temp4.name)
					residues_before = []
					for chain in structure_before[0]:
						residues_before.append(chain.get_list())
					residues_after = []
					for chain in structure_after[0]:
						residues_after.append(chain.get_list())
					chain_counter = ""
					for i,chain in enumerate(structure_before[0]):
						try:
							if structure_after[0].get_list()[i].id != structure_before[0].get_list()[i].id:
								try:
									# HACK BECAUSE OF https://github.com/biopython/biopython/issues/1551
									# Essentially, a new change in biopython prevents you from changing the
									# id to an already existing id which broke this initial script.
									# Therefore, we now change the ids to "change_counter" which will never look
									# like a canonical chainid.
									structure_after[0][structure_before[0].get_list()[i].id].id = chain_counter  
									chain_counter += "KK"
								except KeyError:
									pass
								structure_after[0].get_list()[i].id = structure_before[0].get_list()[i].id
							if len(residues_before[i]) != len(residues_after[i]):
								raise PDBFixerResIdentifiabilityIssue()
	
						# When exceeding chainid Z, pdbfixer has discarded it, whereas biopython has not.
						# For simplicity, we just discard it as well and pretend it does not exist.
						# This is a very rare instance and will likely never be a problem unless you 
						# are extremely unlucky to work with huge proteins where you care about the
						# truncation.
						except IndexError: 
							continue
	
						counter = 99999 # A large residue number that will never exist in a pdb.
						for res1, res2 in zip(residues_before[i], residues_after[i]):
							assert(res1.get_resname().strip() == res2.get_resname().strip() or pdbfixer.pdbfixer.substitutions[res1.get_resname()].strip() == res2.get_resname().strip())
							if res2.id != res1.id:
								try:
									# Similar issue as previous hack https://github.com/biopython/biopython/issues/1551
									structure_after[0][chain.get_id()][res1.id].id = (' ', counter,' ') 
								except KeyError:
									pass
								res2.id = res1.id
								counter += 1
					with open(os.path.join(out_dir, pdbid+"_clean.pdb"), "w") as outpdb:
						pdbio.set_structure(structure_after[0])
						pdbio.save(outpdb)
