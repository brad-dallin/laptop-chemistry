"""Utility functions for notebooks."""

from io import BytesIO
from PIL import Image

import rdkit
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO


def show_acs1996(
    mol: rdkit.Chem.rdchem.Mol,
    legend: str = "",
    width: int = -1,
    height: int = -1,
    highlight_atoms: list = None,
    highlight_atom_colors: dict = None,
    highlight_color: tuple = (1, 1, 0.7),
) -> Image.Image:
    """
    Draw molecule in ACS1996 format and return as PIL Image.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object to visualize
    legend : str, optional
        Text label for the molecule, default: ""
    width : int, optional
        Width of the image in pixels, default: -1
    height : int, optional
        Height of the image in pixels, default: -1
    highlight_atoms : list, optional
        List of atom indices to highlight, default: None
    highlight_atom_colors : dict, optional
        Dictionary mapping atom indices to RGB colors, default: None
    highlight_color : tuple, optional
        RGB color for highlighting atoms, default: (1, 1, 0.7)

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

    # Set up atom highlighting
    highlight_atom_map = {}
    if highlight_atoms and not highlight_atom_colors:
        highlight_atom_map = {atom_idx: highlight_color for atom_idx in highlight_atoms}
    elif highlight_atom_colors:
        highlight_atom_map = highlight_atom_colors

    # Draw the molecule with highlighting
    Draw.DrawMoleculeACS1996(
        d2d,
        mol,
        legend=legend,
        highlightAtoms=highlight_atoms or [],
        highlightAtomColors=highlight_atom_map,
    )
    d2d.FinishDrawing()

    # Convert to PIL Image
    bio = BytesIO(d2d.GetDrawingText())
    return Image.open(bio)
