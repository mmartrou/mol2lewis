# mol2lewis

Convert SMILES strings to ChemFig LaTeX code with Lewis structures, including lone pairs and draft mode annotations.

## Features

- Convert SMILES to ChemFig LaTeX with automatic lone pair placement
- Support for draft mode with visual markers for bonds and lone pairs
- Handles complex molecules with rings, heteroatoms, and various bond types
- Integrates with mol2chemfig for 2D structure generation

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
from mol2lewis import smiles_to_lewis

# Basic usage
chemfig_code = smiles_to_lewis("CCO")
print(chemfig_code)

# Draft mode
draft_code = smiles_to_lewis("CCO", draft=True)
print(draft_code)

# With options
rotated_code = smiles_to_lewis("CCO", angle=45, aromatic_circles=True)
```

## Dependencies

- RDKit: For molecular structure processing
- mol2chemfig: For SMILES to ChemFig conversion (must be installed separately)

## License

MIT License