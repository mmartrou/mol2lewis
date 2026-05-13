# mol2lewis

[![GitHub Repository](https://img.shields.io/badge/GitHub-mmartrou/mol2lewis-blue)](https://github.com/mmartrou/mol2lewis)

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

Each dictionary has four keys: 'normal', 'draft', 'iupac_name', and 'name', corresponding to the ChemFig code in normal/draft modes, the IUPAC name, and the common name (from PubChem).

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
\chemfig{
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
}
```

and 

```python
print(code[0]['draft'])    # chemfig code for the lewis code
```

will give you

```bash
\chemfig{
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
}
```

```python
isomers = lewis("C2H6O", selection='all')  # All isomers of C2H6O
print(isomers[0]['iupac_name'])  # IUPAC name: ethanol
print(isomers[0]['name'])       # Common name: ethanol
print(isomers[1]['iupac_name'])  # IUPAC name: dimethyl ether
print(isomers[1]['name'])       # Common name: methoxymethane
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

...the chemfig molecule code here...

\end{document}
```

The `mol2lewis.sty` LaTeX package is automatically included in the Python package installation for convenience.

**LaTeX Package Installation:**

After installing `mol2lewis` via pip, you can retrieve the `mol2lewis.sty` file from the Python installation:

```bash
# Find where pip installed mol2lewis
pip show mol2lewis

# Copy the .sty file to your LaTeX directory
cp $(python -c "import mol2lewis; import os; print(os.path.dirname(mol2lewis.__file__))")/mol2lewis.sty /your/latex/directory/
```

Alternatively, simply download the file from the GitHub repository and place it in your LaTeX directory.

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
- `angle` (float): Rotation angle (default: 0.0)
- `aromatic_circles` (bool): Draw aromatic circles (default: False)
- `show_carbons` (bool): Show carbon symbols (default: True)
- `show_methyls` (bool): Show methyl groups (default: False)
- `flip` (bool): Horizontal flip (default: False)
- `flop` (bool): Vertical flip (default: False)
- `selection` (str): For formulas, one of `first`, `random`, `all`, `first_n` (default: `first`)
- `n` (int): Number of molecules when `selection='first_n'`
- `enumerate_stereo` (bool): Enumerate stereoisomers locally with RDKit (default: False)
- `wrap_chemfig` (bool): Wrap output in `\\chemfig{...}` (default: True)

**Returns:** A list of dictionaries with 'normal', 'draft', 'iupac_name', and 'name' keys. Each dictionary corresponds to an isomer. Returns an empty list on failure.


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
print(lewis("c1ccccc1", aromatic_circles=True))

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
    print(f"Isomer {i}: {iso['iupac_name']} ({iso['name']})")
    print("Normal:\n", iso["normal"], "\n")
    print("Draft:\n", iso["draft"], "\n")
```

Each element in the `isomers` list is a dictionary with four keys:
- `normal`: chemfig code ready to use in LaTeX
- `draft`: annotated version to visualize lone pairs and charges
- `iupac_name`: IUPAC name of the molecule
- `name`: common name of the molecule

You can copy and paste these codes into a LaTeX document using the chemfig package (see above for details).

For draft mode, make sure to include the `mol2lewis.sty` package in your LaTeX preamble. It's a short file with only this inside (but draft mode can't work without it):

```latex
\ProvidesPackage{mol2lewis}[2025/01/22 v0.2.0 Commands for mol2lewis draft mode]
\RequirePackage{xcolor}

% Draft mode markers for bonds
\newcommand{\red}{\.[{.style={draw=red,fill=red}}]} % red dot for single bond in draft mode
\newcommand{\redd}{\:[{.style={draw=red,fill=red}}]}% red dots for double bond in draft mode
\newcommand{\draftlp}{\:}                           % two dots for lone pair in draft mode

\endinput
```

## License

MIT License