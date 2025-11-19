"""
Mol2Lewis Package
=================

A Python package for converting SMILES strings to LaTeX chemfig code with 
automatic lone pair electron placement based on 2D molecular geometry.

Based on mol2chemfig (https://github.com/Augus1999/mol2chemfigPy3)
with added support for automatic Lewis structure electron placement.

Main Functions:
    smiles_to_lewis(smiles, **kwargs): Convert SMILES to formatted chemfig code
    format_chemfig(code): Format chemfig code for human readability
"""

from .core import (
    smiles2lewis,
    formula2lewis,
    cid2lewis,
    iupac2lewis,
    inchi2lewis,
    inchikey2lewis,
)
from .formatting import format_chemfig

__version__ = "1.0.0"
__all__ = [
    'smiles2lewis',
    'formula2lewis',
    'cid2lewis',
    'iupac2lewis',
    'inchi2lewis',
    'inchikey2lewis',
    'format_chemfig',
]