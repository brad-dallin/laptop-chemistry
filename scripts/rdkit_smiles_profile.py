#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rdkit_smiles_profile.py


Author: Brad Dallin
Date: April 6, 2025
"""


########################################################################
## Imports
########################################################################
import os
import sys
import argparse
import pandas as pd
import rdkit
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

# Suppress RDKit Output
rdkit.RDLogger.DisableLog('rdApp.*')

########################################################################
## Functions
########################################################################
def parse_args(argv: list) -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        prog = "rdkit_smiles_profile.py",
        description = "Load CSV file with a SMILES column, process SMILES and " \
        "calculate 2D RDKit descriptors, and return CSV file with descriptor columns"
    )
    parser.add_argument(
        "input",
        action = "store",
        type = str,
        help = "CSV file with SMILES column",
    )
    parser.add_argument(
        "output",
        action = "store",
        type = str,
        help = "CSV file with processed SMILES and descriptor columns",
    )
    parser.add_argument(
        "-c", "--smiles_column",
        action = "store",
        type = str,
        default = "smiles",
        help = "Column name containing SMILES",
    )
    parser.add_argument(
        "-d", "--descriptors",
        action = "store",
        type = str,
        default = "all",
        help = "Names of RDKit descriptors to calculate",
    )
    args = parser.parse_args(argv)
    return args


# Read molecule
def load_smiles(smi: str) -> rdkit.Chem.rdchem.Mol:
    """Load SMILES as RDKit Molecule"""
    mol = None
    if pd.notna(smi):
        try:
            mol = rdkit.Chem.MolFromSmiles(
                smi,
                sanitize = False
            )
            assert mol is not None
        except Exception as e:
            print(f"Error processing SMILES at index {ii}: {smi}")
            print(f"Exception: {e}")
    return mol


# Sanitize molecule
def sanitize_molecule(mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
    """Sanitize molecules with RDKit, keep largest fragment, and uncharge"""
    largest_frag_app = rdMolStandardize.LargestFragmentChooser()
    uncharge_app = rdMolStandardize.Uncharger()
    try:
        rdkit.Chem.SanitizeMol(mol)
        largest_frag_app.chooseInPlace(mol)
        uncharge_app.unchargeInPlace(mol)
    except Exception as e:
        print(f"Error processing SMILES: {rdkit.Chem.MolToSmiles(mol)}")
        print(f"Exception: {e}")
        return False
    return True


# Calculate descriptors
def calculate_descriptors(
        mol: rdkit.Chem.rdchem.Mol,
        descriptors: list = []
    ) -> dict:
    """Calculate RDKit descriptors"""
    rd_descriptors = {
        "SMILES": rdkit.Chem.MolToSmiles(mol),
    }
    for descriptor, fxn in Descriptors._descList:
        if descriptor in descriptors or descriptors.lower() == "all":
            try:
                value = fxn(mol)
            except Exception as e:
                value = None
            rd_descriptors[f"rdkit_{descriptor}"] = value
    return rd_descriptors

# Process molecules and calculate descriptors
def process_molecules(
        df: pd.DataFrame,
        smiles_column: str = "smiles",
        descriptors: list = []
    ) -> tuple:
    """Preprocess molecules with RDKit by checking compatibility"""
    mols = []
    skipped_indices = []
    for ii, smi in enumerate(df[smiles_column]):
        mol = load_smiles(smi)
        if mol is None:
            skipped_indices.append(ii)
            continue
        status = sanitize_molecule(mol)
        if status is False:
            skipped_indices.append(ii)
            continue
        rddescriptors = calculate_descriptors(
            mol,
            descriptors = descriptors
        )
        rddescriptors["idx"] = ii
        mols.append(rddescriptors)
    print(f"{len(mols)}/{df.shape[0]} molecules processed!")
    print(f"{df.shape[0]-len(mols)}/{df.shape[0]} molecules skipped!")
    return (mols, skipped_indices)


# Main function
def main(argv):
    """Main function"""
    args = parse_args(argv)
    if args.descriptors != "all":
        args.descriptors = [ dd.strip() for dd in args.descriptors.split(",") ]
    df = pd.read_csv(args.input)

    # Process SMILES
    mols = []
    skipped_indices = []
    for ii, smi in enumerate(df[args.smiles_column]):
        mol = load_smiles(smi)
        if mol is None:
            skipped_indices.append(ii)
            continue
        status = sanitize_molecule(mol)
        if status is False:
            skipped_indices.append(ii)
            continue
        rddescriptors = calculate_descriptors(
            mol,
            descriptors = args.descriptors
        )
        rddescriptors["idx"] = ii
        mols.append(rddescriptors)
    print(f"{len(mols)}/{df.shape[0]} molecules processed!")
    print(f"{df.shape[0]-len(mols)}/{df.shape[0]} molecules skipped!")

    # Update dataframe with RDKit descriptors
    mol_df = pd.DataFrame(mols)
    mol_df.set_index("idx", inplace=True)

    mol_df.columns = [ 
        col if col == "SMILES" else col.lower() 
        for col in mol_df.columns]

    df = pd.concat([df, mol_df], axis=1)
    cols = df.columns.tolist()
    cols.remove('SMILES')
    cols.insert(0, 'SMILES')
    df = df[cols]

    # Save to file
    df.to_csv(
        args.output,
        index=False,
        encoding="utf-8"
    )
    return 0


########################################################################
## Run
########################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
