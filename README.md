# mol2lewis

Convert various chemical identifiers (SMILES, formulas, CIDs, names, InChI, etc.) to ChemFig LaTeX code with Lewis structures, including lone pairs and draft mode annotations.

## Features

- Convert SMILES, molecular formulas, PubChem CIDs, IUPAC names, InChI, or InChIKey to ChemFig LaTeX with automatic lone pair placement
- Support for draft mode with visual markers for bonds and lone pairs
- Handles complex molecules with rings, heteroatoms, and various bond types
- Integrates with mol2chemfig for 2D structure generation and PubChem for data fetching

## Installation

### Prerequisites

- Python 3.8+
- RDKit
- mol2chemfig (command-line tool)

### Install

```bash
pip install mol2lewis
```

Note: You also need to install mol2chemfig separately, as it's a command-line tool.

## Usage

```python
from mol2lewis import smiles2lewis, formula2lewis, cid2lewis, iupac2lewis, inchi2lewis, inchikey2lewis

# From SMILES
chemfig_code = smiles2lewis("CCO")  # Ethanol

# From molecular formula (gets first isomer from PubChem)
formula_code = formula2lewis("C2H6O")  # Ethanol

# From formula with options (get all isomers)
all_isomers = formula2lewis("C2H6O", selection='all')  # List of ChemFig codes

# From PubChem CID
cid_code = cid2lewis(702)  # Ethanol

# From IUPAC name
iupac_code = iupac2lewis("ethanol")  # Ethanol

# From InChI
inchi_code = inchi2lewis("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")  # Ethanol

# From InChIKey
inchikey_code = inchikey2lewis("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")  # Ethanol

# Draft mode
draft_code = smiles2lewis("CCO", draft=True)

# With options
rotated_code = smiles2lewis("CCO", angle=45, aromatic_circles=True)
```
## API Reference

- smiles2lewis(smiles, **options): Convert SMILES string.
- formula2lewis(formula, selection='first', n=None, **options): Convert molecular formula. selection can be 'first', 'random', 'all', or 'first_n' (with n).
- cid2lewis(cid, **options): Convert PubChem CID.
- iupac2lewis(name, **options): Convert IUPAC name.
- inchi2lewis(inchi, **options): Convert InChI string.
- inchikey2lewis(inchikey, **options): Convert InChIKey.

All functions support mol2chemfig options (e.g., angle, draft, show_carbons) and return ChemFig LaTeX code.

## Dependencies

- RDKit: For molecular structure processing.
- mol2chemfig: For SMILES to ChemFig conversion (must be installed separately).
- pubchempy: For fetching data from PubChem.

## License

MIT License

### Notes
- **Error Handling**: Each function raises `ValueError` if PubChem lookup fails.
- **Multiple Results**: `formula2lewis` can return lists for 'all' or 'first_n'.
- **Testing**: I'd test with real examples like CID 702 (ethanol) or formula "C2H6O".
- **PubChem Support**: IUPAC names support multiple languages if PubChem has them; InChI/InChIKey are standard.

If you'd like me to adjust anything (e.g., add more options, change return formats, or update the README further), just let me know! This should make the package much more versatile. 😊