#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Visualization functions for molecular structures.

This module provides functions for visualizing molecular structures
using RDKit and other visualization libraries.
"""

import rdkit
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO
from typing import Optional, Union


def show_acs1996(
    mol: rdkit.Chem.rdchem.Mol, 
    legend: str = "",
    width: int = -1,
    height: int = -1
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
        Width of the image in pixels, by default -1 (auto)
    height : int, optional
        Height of the image in pixels, by default -1 (auto)
        
    Returns
    -------
    PIL.Image.Image
        PIL Image object containing the rendered molecule
        
    Examples
    --------
    >>> from rdkit import Chem
    >>> from cheminformatics.visualization import show_acs1996
    >>> mol = Chem.MolFromSmiles("CCO")
    >>> img = show_acs1996(mol, legend="Ethanol")
    >>> img.save("ethanol.png")
    """
    rdkit.Chem.rdDepictor.Compute2DCoords(mol)
    rdkit.Chem.rdDepictor.StraightenDepiction(mol)
    d2d = Draw.MolDraw2DCairo(width, height)
    Draw.DrawMoleculeACS1996(d2d, mol, legend=legend)
    bio = BytesIO(d2d.GetDrawingText())
    return Image.open(bio)