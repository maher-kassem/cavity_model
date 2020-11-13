import tempfile
from sys import argv, stdout
import argparse

from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *


def run_simulation(
    pdb_filename: str,
    n_nanoseconds: int,
    x_dim: float = 7.0,
    y_dim: float = 5.0,
    z_dim: float = 5.0,
    n_report_steps: int = 25000,
    output_traj_suffix: str = "_md_traj.pdb",
) -> None:

    """
    Function to run explicit solvent MD with the amber14 force field using tip3 model. 
    Default parameters were chosen where available. The integrator takes steps of 2 femtoseconds.

    Paramaters
    ----------
    pdb_filename: str
        Filename of pdb
    n_nanoseconds: int
        Number of nanoseconds to run the simulation for
    x_dim: float
        x dimension size of box
    y_dim: float
        y dimension size of box
    y_dim: float
        y dimension size of box
    n_report_steps: int
        Number of steps before logging a new sampled structure
    output_traj_suffix: str
        Name of suffix for output trajectory PDB file
    """

    box_size = [float(x_dim), float(y_dim), float(z_dim)]
    n_steps = n_nanoseconds * 1000000  # converting nano seconds to femto seconds

    # Load pdb
    pdb = PDBFile(pdb_filename)
    output_pdb = pdb_filename.strip(".pdb") + output_traj_suffix

    # Setup protein system
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(
        forcefield, boxSize=Vec3(box_size[0], box_size[1], box_size[2]) * nanometers
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * nanometer,
        constraints=HBonds,
    )

    # Simulation
    integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    simulation.reporters.append(PDBReporter(output_pdb, n_report_steps))
    simulation.reporters.append(
        StateDataReporter(
            stdout, n_report_steps, step=True, potentialEnergy=True, temperature=True
        )
    )
    simulation.step(n_steps)


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=str)
    parser.add_argument("--ns", type=int, default=20)  # nanoseconds
    parser.add_argument("--x", type=float, default=7.0)  # nanometers
    parser.add_argument("--y", type=float, default=5.0)  # nanometers
    parser.add_argument("--z", type=float, default=5.0)  # nanometers
    args_dict = vars(parser.parse_args())

    # Input arguments
    pdb_filename = args_dict["pdb"]
    n_nanoseconds = args_dict["ns"]
    x_dim = args_dict["x"]
    y_dim = args_dict["y"]
    z_dim = args_dict["z"]

    run_simulation(pdb_filename, n_nanoseconds, x_dim, y_dim, z_dim)
