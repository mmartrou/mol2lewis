# mol2lewis

Convert various chemical identifiers (SMILES, formulas, CIDs, names, InChI, InChIKey) to ChemFig LaTeX code with Lewis structures, including lone pairs and draft mode annotations. 
It's based on [mol2chemfig package](https://ctan.org/pkg/mol2chemfig) and the [python3 version](https://github.com/Augus1999/mol2chemfigPy3). It only works for relatively easy molecules. It can't be working with molecules with hundreads of elements.

## Features

- **Universal Input**: Automatically recognizes SMILES, PubChem CIDs, InChI, InChIKey, IUPAC names, and file paths
- **Automatic Lone Pairs**: Adds proper lone pair placement based on 2D molecular geometry
- **Draft Mode**: Visual markers for bonds and lone pairs during development
- **Powered by mol2chemfig**: Uses the mol2chemfig command-line tool for structure generation

## Installation

### Prerequisites

- Python 3.8+
- RDKit
- mol2chemfig command-line tool (installed via mol2chemfigPy3 package)


### Install

```bash
pip install mol2lewis
```

This will automatically install RDKit, mol2chemfigPy3 (which includes the mol2chemfig command), and pubchempy.

## Usage

The lewis command returns a list of python dictionaries. One dictionary for each isomer.
So, if only one isomer is available, only one dictionary will appear in the list.

Each dictionary has two keys: 'normal' and 'draft' keys, corresponding to the two output modes.

Draft mode adds red visual markers to bonds for educational purposes.

### Simple Examples

```python
from mol2lewis import lewis

# From SMILES
code = lewis("CCO")  # Ethanol
print(code[0]['normal'])    # chemfig code for the lewis code
```

will give you
```bash
H
-[:343.9]C
(
    -[:73.9]H
)
(
    -[:253.9]H
)
-[:343.9]C
(
    -[:253.9]H
)
(
    -[:73.9]H
)
-[:343.9]\charge{224=\|,344=\|}{O}
-[:43.9]H
```

and 

```python
print(code[0]['draft'])    # chemfig code for the lewis code
```

will give you

```bash
\charge{344=\red}{H}
-[:343.9]\charge{344=\red,164=\red,74=\red,254=\red}{C}
(
    -[:73.9]\charge{254=\red}{H}
)
(
    -[:253.9]\charge{74=\red}{H}
)
-[:343.9]\charge{164=\red,344=\red,254=\red,74=\red}{C}
(
    -[:253.9]\charge{74=\red}{H}
)
(
    -[:73.9]\charge{254=\red}{H}
)
-[:343.9]\charge{164=\red,44=\red,224=\:,344=\:}{O}
-[:43.9]\charge{224=\red}{H}
```

The same `code` can be obtained from other input types:

```python
# From PubChem CID (numeric)
code = lewis(702)  # Ethanol

# From InChI
code = lewis("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")  # Ethanol

# From InChIKey
code = lewis("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")  # Ethanol

# From IUPAC name (via PubChem lookup)
code = lewis("ethanol")

# From file
code = lewis("molecule.mol")
code = lewis("structure.smi")

# From molecular formula (various isomers can be returned)
code = lewis("C2H6O")  # First isomer
codes = lewis("C2H6O", selection='all')  # All isomers
codes = lewis("C2H6O", selection='first_n', n=3)  # First 3 isomers
```


**LaTeX Requirements:**

If you objective is only to get the pdf of the molecule, you can just compile the LaTeX code below :
```latex
\documentclass{standalone}
\usepackage{chemfig}
\usepackage{mol2lewis}  % Required only for draft mode
\begin{document}
\chemfig{
    ...the molecule code here...
    }
\end{document}
```

### With Options

```python
# Rotate molecule
rotated = lewis("CCO", angle=45)

# Aromatic circles (default: False)
aromatic = lewis("c1ccccc1", aromatic_circles=True)

# Hide carbon atoms (default: False)
no_carbons = lewis("CCO", show_carbons=False)

# Combine options
code = lewis(702, angle=30, aromatic_circles=True, draft=True)
```


## Python API Reference

### Main Function

**`lewis(input, **options)`**

Universal converter that automatically recognizes input type (SMILES, CID, InChI, InChIKey, name, or file path).

**Parameters:**
- `input` (str or int): Chemical identifier
- `rotate` (float): Rotation angle (default: 0.0)
- `aromatic_circles` (bool): Draw aromatic circles (default: False)
- `relative_angle` (bool): Use relative angles (default: False)
- `show_carbons` (bool): Show carbon symbols (default: True)

**Returns:** A list of dictionaries with 'normal' and 'draft' keys. Each dictionary corresponds to an isomer.


## How It Works

1. **Input Recognition**: Uses mol2chemfigPy3's automatic input type detection
2. **ChemFig Generation**: Converts to basic ChemFig code
3. **Molecular Analysis**: Parses structure with RDKit for geometry
4. **Lone Pair Placement**: Calculates and adds `\charge{}` annotations
5. **Formatting**: Produces clean, readable LaTeX output

## Dependencies

- **RDKit**: Molecular structure processing and geometry
- **mol2chemfigPy3**: SMILES/InChI/CID to ChemFig conversion
- **pubchempy**: PubChem database queries (for formulas and names)

## Examples

### Basic Usage

```python
from mol2lewis import lewis

# Water
print(lewis("O"))

# Benzene with aromatic circle
print(lewis("c1ccccc1", aromatic=True))

# Ethanol from different sources
print(lewis("CCO"))           # SMILES
print(lewis(702))             # CID
print(lewis("ethanol"))       # Name
```

### Advanced: Multiple Isomers

To get all structural isomers for a given molecular formula (e.g. $C_2H_6O$), use the option `selection='all'`

```python
from mol2lewis import lewis

# List all isomers of C2H6O
isomers = lewis("C2H6O", selection="all")

for i, iso in enumerate(isomers, 1):
    print(f"Isomer {i}:")
    print("Normal:\n", iso["normal"], "\n")
    print("Draft:\n", iso["draft"], "\n")
```

Each element in the `isomers` list is a dictionary with two keys:
- `normal`: chemfig code ready to use in LaTeX
- `draft`: annotated version to visualize lone pairs and charges

You can copy and paste these codes into a LaTeX document using the chemfig package (see above for details).

## License

MIT License