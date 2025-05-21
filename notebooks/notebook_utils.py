"""Utility functions for notebooks."""

import rdkit
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO


def show_acs1996(
    mol: rdkit.Chem.rdchem.Mol,
    legend: str = "",
    width: int = -1,
    height: int = -1,
) -> Image.Image:
    """
    Draw molecule in ACS1996 format and return as PIL Image.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object to visualize
    legend : str, optional
        Text label for the molecule, by default ""
    width : int, optional
        Width of the image in pixels, by default -1
    height : int, optional
        Height of the image in pixels, by default -1

    Returns
    -------
    PIL.Image.Image
        PIL Image object containing the rendered molecule
    """
    # Compute 2D coordinates
    rdkit.Chem.rdDepictor.Compute2DCoords(mol)
    rdkit.Chem.rdDepictor.StraightenDepiction(mol)

    # Create drawer with specified dimensions
    d2d = Draw.MolDraw2DCairo(width, height)

    # Draw the molecule
    Draw.DrawMoleculeACS1996(d2d, mol, legend=legend)
    d2d.FinishDrawing()

    # Convert to PIL Image
    bio = BytesIO(d2d.GetDrawingText())
    return Image.open(bio)
