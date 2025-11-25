"""
Mol2Lewis Package
=================

A Python package for converting various chemical identifiers (SMILES, CIDs, 
InChI, InChIKey, formulas, names) to LaTeX chemfig code with automatic lone 
pair electron placement based on 2D molecular geometry.

Based on mol2chemfigPy3 (https://github.com/Augus1999/mol2chemfigPy3)
with added support for automatic Lewis structure electron placement.

Main Functions:
    lewis(input, **kwargs): Convert any chemical identifier to chemfig with lone pairs
    formula2lewis(formula, **kwargs): Convert molecular formula (queries PubChem)
    format_chemfig(code): Format chemfig code for human readability
"""

from .core import (
    lewis,
)
from .formatting import format_chemfig

__version__ = "1.0.3"
__all__ = [
    'lewis',
    'format_chemfig',
]