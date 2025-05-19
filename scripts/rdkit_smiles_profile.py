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
from typing import List, Dict, Tuple, Optional
import rdkit
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

# Suppress RDKit Output
rdkit.RDLogger.DisableLog('rdApp.*')

########################################################################
## Functions
########################################################################
def parse_args(argv: List[str]) -> argparse.Namespace:
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
        "-c",
        "--smiles_column",
        action = "store",
        type = str,
        default = "smiles",
        help = "Column name containing SMILES",
    )
    parser.add_argument(
        "-d",
        "--descriptors",
        action = "store",
        type = str,
        default = "all",
        help = "Names of RDKit descriptors to calculate",
    )
    args=parser.parse_args(argv)
    return args


# Read molecule
def load_smiles(smi: str) -> rdkit.Chem.rdchem.Mol:
    """Load SMILES as RDKit Molecule"""
    mol = None
    if pd.notna(smi):
        try:
            mol = rdkit.Chem.MolFromSmiles(
                smi,
                sanitize=False
            )
            assert mol is not None
        except Exception as e:
            print(f"Error processing SMILES: {smi}")
            print(f"Exception: {e}")
    return mol


# Sanitize molecule
def sanitize_molecule(mol: rdkit.Chem.rdchem.Mol) -> bool:
    """Sanitize molecules with RDKit, keep largest fragment, and uncharge"""
    try:
        rdkit.Chem.SanitizeMol(mol)
        largest_frag_app = rdMolStandardize.LargestFragmentChooser()
        uncharge_app = rdMolStandardize.Uncharger()
        largest_frag_app.chooseInPlace(mol)
        uncharge_app.unchargeInPlace(mol)
        return True
    except Exception as e:
        print(f"Error processing SMILES: {rdkit.Chem.MolToSmiles(mol)}")
        print(f"Exception: {e}")
        return False


# Calculate descriptors
def calculate_descriptors(
        mol: rdkit.Chem.rdchem.Mol,
        descriptors: List[str] = ["all"]
    ) -> Dict[str, float]:
    """Calculate RDKit descriptors"""
    rd_descriptors = {
        "SMILES": rdkit.Chem.MolToSmiles(mol),
    }
    for descriptor, fxn in Descriptors._descList:
        if descriptor in descriptors or descriptors[0].lower() == "all":
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
        descriptors: List[str] = ["all"]
    ) -> Tuple[List[Dict[str, float]], List[int]]:
    """Preprocess molecules with RDKit by checking compatibility"""
    mols = []
    unique_smiles = []
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
        frag_smi = rdkit.Chem.MolToSmiles(mol)
        if frag_smi in unique_smiles:
            skipped_indices.append(ii)
            continue
        rddescriptors = calculate_descriptors(
            mol,
            descriptors = descriptors
        )
        rddescriptors["idx"] = ii
        mols.append(rddescriptors)
        unique_smiles.append(frag_smi)
    print(f"{len(mols)}/{df.shape[0]} molecules processed!")
    print(f"{df.shape[0]-len(mols)}/{df.shape[0]} molecules skipped!")
    return (mols, skipped_indices)


# Main function
def main(argv: List[str]) -> int:
    """Main function"""
    args = parse_args(argv)
    args.descriptors = [dd.strip() for dd in args.descriptors.split(",")]
    df = pd.read_csv(args.input)

    # Process molecules
    mols, skipped_indices = process_molecules(
        df,
        smiles_column=args.smiles_column,
        descriptors=args.descriptors
    )

    # Update dataframe with RDKit descriptors
    mol_df = pd.DataFrame(mols)
    mol_df.set_index("idx", inplace=True)

    mol_df.columns = [col if col == "SMILES" else col.lower()
                      for col in mol_df.columns]

    df = pd.concat([df, mol_df], axis=1)
    cols = df.columns.tolist()
    cols.remove('SMILES')
    cols.insert(0, 'SMILES')
    df = df[cols]

    # Drop skipped indices
    df = df.drop(index=skipped_indices)

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
