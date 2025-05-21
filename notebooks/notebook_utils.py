"""Utility functions for notebooks."""

import rdkit
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO


def show_acs1996(
    mol: rdkit.Chem.rdchem.Mol, 
    legend: str = "",
    width: int = -1,
    height: int = -1
) -> Image.Image:
    """
    Draw molecule in ACS1996 format and return as PIL Image.
    """
    rdkit.Chem.rdDepictor.Compute2DCoords(mol)
    rdkit.Chem.rdDepictor.StraightenDepiction(mol)
    d2d = Draw.MolDraw2DCairo(width, height)
    Draw.DrawMoleculeACS1996(d2d, mol, legend=legend)
    bio = BytesIO(d2d.GetDrawingText())
    return Image.open(bio)