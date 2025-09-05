#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
openmm-basic-md.py

This script runs a basic MD simulation starting with a solvated PDB file.

Input:
- Solvated PDB file

Output:
- DCD trajectory file

Dependencies:
- os
- sys
- argparse
- openmm

Usage:
- Run the script using Python 3.

Author: Brad Dallin
Date: Aug. 7, 2025
"""

########################################################################
## Imports
########################################################################
import os
import sys
import argparse
from openmm import app
from openmm import unit
from openmm import LangevinMiddleIntegrator


########################################################################
## Functions
########################################################################
def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        prog="qry_sqldb.py",
        description="Run a basic MD simulation starting with a PDB file.",
    )
    parser.add_argument(
        "input",
        action="store",
        type=str,
        dest="input",
        help="Input PDB file",
    )
    parser.add_argument(
        "-s", "-reporter_step_size",
        action="store",
        type=int,
        dest="reporter_step_size",
        default=1000,
        required=False,
        help="Report every N MD steps"
    )
    parser.add_argument(
        "-n", "-total_steps",
        action="store",
        type=int,
        dest="total_steps",
        default=10000,
        required=False,
        help="Total number of MD steps"
    )
    parser.add_argument(
        "-t", "-tail_file",
        action="store_true",
        dest="tail_file",
        required=False,
        help="Write log to stdout"
    )
    args = parser.parse_args(argv)
    return args


# Create simulation object
def create_simulation(input_file: str) -> app.Simulation:
    """Create a MD simulation object"""
    # Load molecular structure
    pdb = app.PDBFile(input_file)
    # Define force field parameters
    forcefield = app.ForceField(
        'amber19-all.xml', 'amber19/tip3pfb.xml'
    )
    # Create physical system with simulation parameters
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1*unit.nanometer,
        constraints=app.HBonds
    )
    # Set up Langevin integrator for constant temperature dynamics
    integrator = LangevinMiddleIntegrator(
        300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds
    )
    # Create simulation object and set initial positions
    simulation = app.Simulation(
        pdb.topology, system, integrator
    )
    simulation.context.setPositions(pdb.positions)
    return simulation


# Main function
def main(argv: list[str]) -> int:
    """Main function"""
    try:
        args = parse_args(argv)
        basename, _ = os.path.splitext(args.input)
        trj_file = basename + "_trj.dcd"
        csv_file = basename + "_state.csv"
        # Create simulation object
        simulation = create_simulation(args.input)
        # Minimize energy to remove bad contacts
        simulation.minimizeEnergy()
        # Set up trajectory and progress reporting
        simulation.reporters.append(
            app.DCDReporter(trj_file, args.reporter_step_size)
        )
        simulation.reporters.append(
            app.StateDataReporter(
                csv_file, args.reporter_step_size,
                step=True,
                potentialEnergy=True,
                temperature=True
            )
        )
        if args.tail_file is True:
            simulation.reporters.append(
                app.StateDataReporter(
                    sys.stdout, args.reporter_step_size,
                    step=True,
                    potentialEnergy=True,
                    temperature=True
                )
            )
        # Run molecular dynamics simulation
        simulation.step(args.total_steps)
        return 0
    except Exception as err:
        print(f"Error: {err}")
        return 1


########################################################################
## Run
########################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
