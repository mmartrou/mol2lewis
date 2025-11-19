"""
Core Module
===========

Main conversion logic from SMILES to chemfig code.
"""

import subprocess
import re
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
import random
import pubchempy as pcp

from .geometry import calculate_charge_entries
from .formatting import format_chemfig


def smiles2lewis(smiles, **mol2chemfig_options):
    """
    Convert a SMILES string to formatted chemfig LaTeX code with lone pairs.
    
    This function:
    1. Parses the SMILES string using RDKit
    2. Adds explicit hydrogens
    3. Generates 2D coordinates
    4. Converts to chemfig via mol2chemfig
    5. Adds \\charge{} annotations for lone pairs based on geometry
    6. Formats the output for readability
    
    Args:
        smiles (str): SMILES string representation of the molecule
        **mol2chemfig_options: Any mol2chemfig command-line options:
            - angle (float): Rotate molecule counterclockwise by this angle (default: 0.0)
            - show_carbons (bool): Show element symbol for carbon atoms (default: True)
            - show_methyls (bool): Show element symbols for methyl groups (default: False)
            - relative_angles (bool): Use relative bond angles (default: False)
            - flip (bool): Flip structure horizontally (default: False)
            - flop (bool): Flip structure vertically (default: False)
            - aromatic_circles (bool): Draw circles instead of double bonds in rings (default: False)
            - fancy_bonds (bool): Draw fancier double and triple bonds (default: False)
            - markers (str): Give atoms/bonds unique markers for arrows (default: None)
            - atom_numbers (bool): Show molfile atom numbers (default: False)
            - bond_scale (str): 'keep', 'scale', or 'normalize' (default: 'normalize')
            - bond_stretch (float): Scaling factor or average for bond lengths (default: 1.0)
            - recalculate_coordinates (bool): Discard existing coords (default: False)
            - hydrogens (str): 'keep', 'add', or 'delete' (default: 'keep')
        
    Returns:
        str: Formatted chemfig code ready to use in \\chemfig{...}, or None on error
        
    Example:
        >>> code = smiles2lewis('CCO')
        >>> code_rotated = smiles2lewis('CCO', angle=45, aromatic_circles=True)
    """
    draft_mode = mol2chemfig_options.pop('draft', False)
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # Add explicit hydrogens and compute 2D coordinates
    mol_with_h = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol_with_h)
    # Kekulize to get explicit single/double bond assignments for aromatics
    # This matches what mol2chemfig does internally
    Chem.Kekulize(mol_with_h, clearAromaticFlags=False)
    # Export to MOL file (temporary)
    Chem.MolToMolFile(mol_with_h, 'temp_molecule.mol')
    # Build mol2chemfig command with options
    cmd = ['mol2chemfig']
    # Map Python kwargs to mol2chemfig command-line options
    option_map = {
        'angle': '--angle',
        'show_carbons': '--show-carbons',
        'show_methyls': '--show-methyls',
        'relative_angles': '--relative-angles',
        'flip': '--flip',
        'flop': '--flop',
        'aromatic_circles': '--aromatic-circles',
        'fancy_bonds': '--fancy-bonds',
        'markers': '--markers',
        'atom_numbers': '--atom-numbers',
        'bond_scale': '--bond-scale',
        'bond_stretch': '--bond-stretch',
        'recalculate_coordinates': '--recalculate-coordinates',
        'hydrogens': '--hydrogens',
    }
    # Set default for show_carbons if not specified
    if 'show_carbons' not in mol2chemfig_options:
        mol2chemfig_options['show_carbons'] = True
    # Add options to command
    for key, value in mol2chemfig_options.items():
        if key in option_map:
            flag = option_map[key]
            if isinstance(value, bool):
                if value:  # Only add flag if True
                    cmd.append(flag)
            else:
                # Add flag with value
                cmd.extend([flag, str(value)])
    cmd.append('temp_molecule.mol')
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return None
    tex_content = result.stdout
    # Remove any deprecated \\lewis macros
    # Remove legacy \lewis macros left by mol2chemfig
    # Two-arg form: \lewis{...}{Se} -> Se
    tex_content = re.sub(r"\\lewis\{[^}]*\}\{([A-Za-z][a-z]?)\}", r"\1", tex_content)
    # One-arg with comma form: \lewis{2.,Se} -> Se
    tex_content = re.sub(r"\\lewis\{[^,}]*,([A-Za-z][a-z]?)\}", r"\1", tex_content)
    # Build charge entries for atoms
    conf = mol_with_h.GetConformer()
    pt = Chem.GetPeriodicTable()
    
    # Get transformation parameters for lone pair adjustment
    rotation_angle = mol2chemfig_options.get('angle', 0)
    flip = mol2chemfig_options.get('flip', False)
    flop = mol2chemfig_options.get('flop', False)
    aromatic_circles_opt = mol2chemfig_options.get('aromatic_circles', False)
    
    # In draft mode, annotate ALL atoms. In normal mode, only heteroatoms.
    if draft_mode:
        INTEREST = {'H', 'C', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Si', 'Ge', 'As', 'Se'}
    else:
        INTEREST = {'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Si', 'Ge', 'As', 'Se'}
    
    charge_dict = {}  # Map atom index to charge string
    sym_to_indices = {}
    for atom in mol_with_h.GetAtoms():
        sym = atom.GetSymbol()
        if sym not in INTEREST:
            continue
        sym_to_indices.setdefault(sym, []).append(atom.GetIdx())
        charge_entry = calculate_charge_entries(atom, conf, pt, rotation_angle, flip, flop, draft=draft_mode)
        if charge_entry:  # Only add non-empty entries
            charge_dict[atom.GetIdx()] = charge_entry
    
    # Replace using atom index from comments (% N)
    if draft_mode:
        # Match ALL atoms including additional main group elements
        # IMPORTANT: Put 2-letter symbols BEFORE 1-letter symbols to avoid partial matches
        pattern = re.compile(
            r'(?P<prefix>-\[:[^\]]+\])?(?P<atom>Cl|Br|Si|Ge|As|Se|I|F|O|N|S|P|C|H|B)(?P<bracks>\[[^\]]*\])?(?P<cmt>%\s*(?P<idx>\d+))?'
        )
    else:
        # Match heteroatoms, expanded set
        # IMPORTANT: Put 2-letter symbols BEFORE 1-letter symbols to avoid partial matches
        pattern = re.compile(
            r'(?P<prefix>-\[:[^\]]+\])?(?P<atom>Cl|Br|Si|Ge|As|Se|I|F|O|N|S|P|B)(?P<bracks>\[[^\]]*\])?(?P<cmt>%\s*(?P<idx>\d+))?'
        )
    
    used_indices = set()

    def _repl(m):
        sym = m.group('atom')
        idx_str = m.group('idx')
        
        # Get atom index from comment
        if idx_str is None:
            # Fallback: if this symbol appears exactly once, use that atom
            indices = sym_to_indices.get(sym, [])
            cand = [i for i in indices if i not in used_indices]
            if len(indices) == 1 and cand:
                atom_idx = cand[0]
            else:
                return m.group(0)
        else:
            # mol2chemfig uses 1-based indexing, RDKit uses 0-based
            atom_idx = int(idx_str) - 1
        charge_str = charge_dict.get(atom_idx, '')
        
        if not charge_str:
            return m.group(0)
        
        atom_token = sym + (m.group('bracks') or '')
        new_core = f"\\charge{{{charge_str}}}{{{atom_token}}}"
        used_indices.add(atom_idx)
        return f"{m.group('prefix') or ''}{new_core}{m.group('cmt') or ''}"
    
    tex_content = pattern.sub(_repl, tex_content)
    # Cleanup: remove accidental charge wrappers around bond-only tokens
    tex_content = re.sub(r"\\charge\{[^}]*\}\{(-\[:[^\]]+\])\}", r"\1", tex_content)
    # Sanitize for chemfig: remove all comments
    no_comments = re.sub(r"%.*", "", tex_content)
    sanitized = no_comments.strip()
    # DRAFT MODE POST-PROCESSING: already handled in geometry.py
    # Lone pairs are annotated with \draftlp\red when draft=True is passed
    # No further post-processing needed here
    # Format for readability
    formatted = format_chemfig(sanitized)
    return formatted


def formula2lewis(formula, selection='first', n=None, lang='en', **options):
    """
    Convert a molecular formula to ChemFig LaTeX with Lewis structures.
    
    Args:
        formula (str): Molecular formula (e.g., 'C2H6O').
        selection (str): 'first' (default), 'random', 'all', or 'first_n'.
        n (int): Number of results for 'first_n' (ignored otherwise).
        lang (str): Language for name ('en' for IUPAC, 'fr' for French, etc.).
        **options: Passed to smiles_to_lewis (e.g., draft=True).
    
    Returns:
        str or list: ChemFig code(s) for the selected SMILES.
    """
    try:
        compounds = pcp.get_compounds(formula, 'formula')
        if not compounds:
            raise ValueError(f"No compounds found for formula '{formula}'")
        
        # Filter out isotopically labeled compounds (e.g., [13C], [2H])
        import re
        smiles_list = []
        seen = set()
        for c in compounds:
            if c.smiles and not re.search(r'\[\d+', c.smiles):
                if c.smiles not in seen:
                    smiles_list.append(c.smiles)
                    seen.add(c.smiles)
        if not smiles_list:
            raise ValueError(f"No SMILES found for formula '{formula}'")
        
        if selection == 'first':
            selected_smiles = smiles_list[0]
        elif selection == 'random':
            selected_smiles = random.choice(smiles_list)
        elif selection == 'all':
            selected_smiles = smiles_list
        elif selection == 'first_n':
            if n is None:
                n = 1
            selected_smiles = smiles_list[:n]
        else:
            raise ValueError("Invalid selection. Use 'first', 'random', 'all', or 'first_n'.")
        
        if isinstance(selected_smiles, list):
            return [smiles2lewis(smiles, **options) for smiles in selected_smiles]
        else:
            return smiles2lewis(selected_smiles, **options)
    except Exception as e:
        raise ValueError(f"Error fetching data for formula '{formula}': {e}")

def cid2lewis(cid, **options):
    """
    Convert a PubChem CID to ChemFig LaTeX with Lewis structures.
    
    Args:
        cid (int or str): PubChem Compound ID.
        **options: Passed to smiles2lewis.
    
    Returns:
        str: ChemFig code.
    """
    try:
        compound = pcp.Compound.from_cid(cid)
        if not compound or not compound.smiles:
            raise ValueError(f"No SMILES found for CID '{cid}'")
        return smiles2lewis(compound.smiles, **options)
    except Exception as e:
        raise ValueError(f"Error fetching data for CID '{cid}': {e}")

def iupac2lewis(name, **options):
    """
    Convert an IUPAC name to ChemFig LaTeX with Lewis structures.
    
    Args:
        name (str): IUPAC name (supports multiple languages if available in PubChem).
        **options: Passed to smiles2lewis.
    
    Returns:
        str: ChemFig code.
    """
    try:
        compounds = pcp.get_compounds(name, 'name')
        if not compounds or not compounds[0].smiles:
            raise ValueError(f"No SMILES found for IUPAC name '{name}'")
        return smiles2lewis(compounds[0].smiles, **options)
    except Exception as e:
        raise ValueError(f"Error fetching data for IUPAC name '{name}': {e}")

def inchi2lewis(inchi, **options):
    """
    Convert an InChI string to ChemFig LaTeX with Lewis structures.
    
    Args:
        inchi (str): InChI string.
        **options: Passed to smiles2lewis.
    
    Returns:
        str: ChemFig code.
    """
    try:
        compounds = pcp.get_compounds(inchi, 'inchi')
        if not compounds or not compounds[0].smiles:
            raise ValueError(f"No SMILES found for InChI '{inchi}'")
        return smiles2lewis(compounds[0].smiles, **options)
    except Exception as e:
        raise ValueError(f"Error fetching data for InChI '{inchi}': {e}")

def inchikey2lewis(inchikey, **options):
    """
    Convert an InChIKey to ChemFig LaTeX with Lewis structures.
    
    Args:
        inchikey (str): InChIKey string.
        **options: Passed to smiles2lewis.
    
    Returns:
        str: ChemFig code.
    """
    try:
        compounds = pcp.get_compounds(inchikey, 'inchikey')
        if not compounds or not compounds[0].smiles:
            raise ValueError(f"No SMILES found for InChIKey '{inchikey}'")
        return smiles2lewis(compounds[0].smiles, **options)
    except Exception as e:
        raise ValueError(f"Error fetching data for InChIKey '{inchikey}': {e}")