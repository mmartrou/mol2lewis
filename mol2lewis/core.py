"""
Core Module
===========

Main conversion logic from various chemical identifiers to chemfig code with Lewis structures.
"""

import re
import os
import subprocess
import tempfile
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)
from rdkit.Chem import rdmolops
import random
import pubchempy as pcp

from .geometry import calculate_charge_entries
from .formatting import format_chemfig


def _is_chemical_formula(s):
    return bool(re.match(r'^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*$', s))


def _generate_chemfig(smiles, **options):
    """
    Generate chemfig code for a SMILES string, returning both normal and draft versions.
    """
    # Parse SMILES with RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Add explicit hydrogens and compute 2D coordinates
    mol_with_h = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol_with_h)
    Chem.Kekulize(mol_with_h, clearAromaticFlags=False)
    
    # Save to temporary MOL file for mol2chemfig
    with tempfile.NamedTemporaryFile(mode='w', suffix='.mol', delete=False) as f:
        temp_molfile = f.name
        Chem.MolToMolFile(mol_with_h, temp_molfile)
    
    try:
        # Build mol2chemfig command
        cmd = ['mol2chemfig']

        # Add atom numbers for matching
        cmd.append('--atom-numbers')

        # Map options to mol2chemfig flags
        if options.get('angle', 0) != 0:
            cmd.extend(['--angle', str(options['angle'])])
        if options.get('show_carbons', True):
            cmd.append('--show-carbons')
        if options.get('show_methyls', False):
            cmd.append('--show-methyls')
        if options.get('flip', False):
            cmd.append('--flip')
        if options.get('flop', False):
            cmd.append('--flop')
        if options.get('aromatic_circles', False):
            cmd.append('--aromatic-circles')
        # Ajout de l'option --wrap-chemfig (True par défaut)
        if options.get('wrap_chemfig', True):
            cmd.append('--wrap-chemfig')

        cmd.append(temp_molfile)

        # Run mol2chemfig
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return None

        tex_content = result.stdout.strip()

    finally:
        # Clean up temp file
        if os.path.exists(temp_molfile):
            os.unlink(temp_molfile)
    
    # Remove any deprecated \lewis macros
    tex_content = re.sub(r"\\lewis\{[^}]*\}\{([A-Za-z][a-z]?)\}", r"\1", tex_content)
    tex_content = re.sub(r"\\lewis\{[^,}]*,([A-Za-z][a-z]?)\}", r"\1", tex_content)
    
    # Convert \mcfatomno{N} macros to % N comments for matching
    # Pattern: \mcf..{...\mcfatomno{N}...}{ATOM} -> ATOM % N
    tex_content = re.sub(
        r'\\mcf(?:left|right|above|below)\{[^}]*\\mcfatomno\{(\d+)\}[^}]*\}\{([^}]+)\}',
        lambda m: f"{m.group(2)} % {m.group(1)}",
        tex_content
    )
    # Pattern: \mcf..{ATOM}{...\mcfatomno{N}...} -> ATOM % N
    tex_content = re.sub(
        r'\\mcf(?:left|right|above|below)\{([^}]+)\}\{[^}]*\\mcfatomno\{(\d+)\}[^}]*\}',
        lambda m: f"{m.group(1)} % {m.group(2)}",
        tex_content
    )
    # Standalone \mcfatomno{N}
    tex_content = re.sub(r'\\mcfatomno\{(\d+)\}', r'% \1', tex_content)
    
    # Build charge entries for atoms
    conf = mol_with_h.GetConformer()
    pt = Chem.GetPeriodicTable()
    
    # Get transformation parameters
    rotation_angle = options.get('angle', 0)
    flip = options.get('flip', False)
    flop = options.get('flop', False)
    aromatic_circles = options.get('aromatic_circles', False)
    
    # Generate normal version (draft=False)
    INTEREST_normal = {'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Si', 'Ge', 'As', 'Se'}
    
    charge_dict_normal = {}
    sym_to_indices_normal = {}
    for atom in mol_with_h.GetAtoms():
        sym = atom.GetSymbol()
        if sym not in INTEREST_normal:
            continue
        sym_to_indices_normal.setdefault(sym, []).append(atom.GetIdx())
        charge_entry = calculate_charge_entries(atom, conf, pt, rotation_angle, flip, flop, draft=False, aromatic_circles=aromatic_circles)
        if charge_entry:
            charge_dict_normal[atom.GetIdx()] = charge_entry
    
    # Match atoms for normal
    pattern_normal = re.compile(
        r'(?P<prefix>-\[:[^\]]+\])?(?P<atom>Cl|Br|Si|Ge|As|Se|I|F|O|N|S|P|B)(?P<bracks>\[[^\]]*\])?(?P<cmt>\s*%\s*(?P<idx>\d+))?'
    )
    
    used_indices_normal = set()

    def _repl_normal(m):
        sym = m.group('atom')
        idx_str = m.group('idx')
        
        if idx_str is None:
            indices = sym_to_indices_normal.get(sym, [])
            cand = [i for i in indices if i not in used_indices_normal]
            if len(indices) == 1 and cand:
                atom_idx = cand[0]
            else:
                return m.group(0)
        else:
            atom_idx = int(idx_str) - 1
        
        charge_str = charge_dict_normal.get(atom_idx, '')
        if not charge_str:
            return m.group(0)
        
        atom_token = sym + (m.group('bracks') or '')
        new_core = f"\\charge{{{charge_str}}}{{{atom_token}}}"
        used_indices_normal.add(atom_idx)
        return f"{m.group('prefix') or ''}{new_core}{m.group('cmt') or ''}"
    
    tex_normal = pattern_normal.sub(_repl_normal, tex_content)
    
    # Cleanup for normal
    tex_normal = re.sub(r"\\charge\{[^}]*\}\{(-\[:[^\]]+\])\}", r"\1", tex_normal)
    if not options.get('aromatic_circles', False):
        tex_normal = re.sub(r'\(\s*\n\s*-\[[^\]]*draw=[^]]*\][^\\]*\\mcfcringle\{[^}]*\}\s*\n\s*\)\s*\n', '', tex_normal)
    
    sanitized_normal = re.sub(r"\s*%\s*[^\n]*", "", tex_normal).strip()
    formatted_normal = format_chemfig(sanitized_normal)
    
    # Generate draft version (draft=True)
    INTEREST_draft = {'H', 'C', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Si', 'Ge', 'As', 'Se'}
    
    charge_dict_draft = {}
    sym_to_indices_draft = {}
    for atom in mol_with_h.GetAtoms():
        sym = atom.GetSymbol()
        if sym not in INTEREST_draft:
            continue
        sym_to_indices_draft.setdefault(sym, []).append(atom.GetIdx())
        charge_entry = calculate_charge_entries(atom, conf, pt, rotation_angle, flip, flop, draft=True, aromatic_circles=aromatic_circles)
        if charge_entry:
            charge_dict_draft[atom.GetIdx()] = charge_entry
    
    # Match atoms for draft
    pattern_draft = re.compile(
        r'(?P<prefix>-\[:[^\]]+\])?(?P<atom>Cl|Br|Si|Ge|As|Se|I|F|O|N|S|P|C|H|B)(?P<bracks>\[[^\]]*\])?(?P<cmt>\s*%\s*(?P<idx>\d+))?'
    )
    
    used_indices_draft = set()

    def _repl_draft(m):
        sym = m.group('atom')
        idx_str = m.group('idx')
        
        if idx_str is None:
            indices = sym_to_indices_draft.get(sym, [])
            cand = [i for i in indices if i not in used_indices_draft]
            if len(indices) == 1 and cand:
                atom_idx = cand[0]
            else:
                return m.group(0)
        else:
            atom_idx = int(idx_str) - 1
        
        charge_str = charge_dict_draft.get(atom_idx, '')
        if not charge_str:
            return m.group(0)
        
        atom_token = sym + (m.group('bracks') or '')
        new_core = f"\\charge{{{charge_str}}}{{{atom_token}}}"
        used_indices_draft.add(atom_idx)
        return f"{m.group('prefix') or ''}{new_core}{m.group('cmt') or ''}"
    
    tex_draft = pattern_draft.sub(_repl_draft, tex_content)
    
    # Cleanup for draft
    tex_draft = re.sub(r"\\charge\{[^}]*\}\{(-\[:[^\]]+\])\}", r"\1", tex_draft)
    
    # In draft mode with aromatic_circles, mark aromatic bonds and circles in red
    if options.get('aromatic_circles', False):
        aromatic_atoms = set()
        for atom in mol_with_h.GetAtoms():
            if atom.GetIsAromatic():
                aromatic_atoms.add(atom.GetIdx())
        
        # Color mcfcringle circles red
        tex_draft = re.sub(
            r'(-\[[^\]]*),draw=none([^\]]*\].*?\\mcfcringle)',
            r'\1,draw=red\2',
            tex_draft
        )
        tex_draft = re.sub(
            r'\\mcfcringle\{([^}]+)\}',
            r'\\color{red}{\\mcfcringle{\1}}',
            tex_draft
        )
    
    sanitized_draft = re.sub(r"\s*%\s*[^\n]*", "", tex_draft).strip()
    formatted_draft = format_chemfig(sanitized_draft)
    
    return {'normal': formatted_normal, 'draft': formatted_draft}


def _enumerate_stereo_via_rdkit(smiles):
    """Enumerate all stereoisomers for a given SMILES locally with RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    # Work on a stereo-free copy so RDKit can enumerate every possibility.
    mol_copy = Chem.Mol(mol)
    Chem.RemoveStereochemistry(mol_copy)
    base_smiles = Chem.MolToSmiles(mol_copy, isomericSmiles=False)
    base_mol = Chem.MolFromSmiles(base_smiles)
    if base_mol is None:
        return []

    options = StereoEnumerationOptions(onlyUnassigned=False, unique=True)
    try:
        isomers = EnumerateStereoisomers(base_mol, options=options)
    except Exception:
        return []

    enumerated = []
    for isomer in isomers:
        enumerated.append(Chem.MolToSmiles(isomer, isomericSmiles=True))

    # Fallback: return original SMILES if enumeration produced nothing.
    return enumerated or [smiles]


def lewis(input_string, **options):
    """
    Convert any chemical identifier to formatted chemfig LaTeX code with lone pairs.

    Automatically recognizes input type:
    - SMILES strings (e.g., 'CCO')
    - PubChem CID (numeric, e.g., 702)
    - InChI (e.g., 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')
    - InChIKey (e.g., 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
    - Chemical formulas (e.g., 'C2H4O2')
    - File paths (e.g., 'molecule.mol', 'structure.smi')

    Args:
        input_string (str or int): Chemical identifier
        **options: Options:
            - angle (float): Rotation angle (default: 0.0)
            - show_carbons (bool): Show carbon symbols (default: True)
            - show_methyls (bool): Show methyl symbols (default: False)
            - flip (bool): Flip horizontally (default: False)
            - flop (bool): Flip vertically (default: False)
            - aromatic_circles (bool): Draw aromatic circles (default: False)
            - selection (str): For formulas - 'first', 'random', 'all', or 'first_n' (default: 'first')
            - n (int): For formulas - number of results for 'first_n' (default: None)
            - enumerate_stereo (bool): For formulas - enumerate stereoisomers locally with RDKit
            - wrap_chemfig (bool): Wrap output in \\chemfig{...} (default: True)

    Returns:
        list: List of dicts with 'normal' and 'draft' chemfig codes, or None on error
    """
    # Convert input to SMILES
    smiles = None
    
    # Try to interpret the input - priority order: CID, InChI, InChIKey, File, SMILES, then Formula
    try:
        cid = int(input_string)
        # It's a CID - fetch from PubChem
        compound = pcp.Compound.from_cid(cid)
        smiles = compound.smiles if compound else None
    except (ValueError, TypeError):
        # String input - check for special formats first
        input_str = str(input_string)
        
        if input_str.startswith('InChI='):
            # InChI format
            try:
                mol = Chem.MolFromInchi(input_str)
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
            except:
                pass
        elif '-' in input_str and len(input_str) == 27 and input_str.count('-') == 2:
            # Likely an InChIKey (format: XXXXXXXXXXXXXX-XXXXXXXXXX-X)
            try:
                compounds = pcp.get_compounds(input_str, 'inchikey')
                if compounds and len(compounds) > 0:
                    smiles = compounds[0].smiles
            except:
                pass
        elif os.path.isfile(input_str):
            # File path - read molecular structure from file
            try:
                if input_str.endswith('.mol'):
                    mol = Chem.MolFromMolFile(input_str)
                elif input_str.endswith('.smi'):
                    with open(input_str, 'r') as f:
                        smiles_str = f.read().strip().split()[0]  # First word is SMILES
                        mol = Chem.MolFromSmiles(smiles_str)
                elif input_str.endswith('.sdf'):
                    suppl = Chem.SDMolSupplier(input_str)
                    mol = next(suppl) if suppl else None
                else:
                    mol = None
                
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
            except:
                pass
        else:
            # Try as SMILES first (even if it looks like a formula)
            # This prioritizes SMILES interpretation (e.g., "CCO" = ethanol, not C2O)
            RDLogger.DisableLog('rdApp.*')
            mol = Chem.MolFromSmiles(input_str)
            RDLogger.EnableLog('rdApp.*')
            if mol is not None:
                smiles = input_str
            elif _is_chemical_formula(input_str):
                # If SMILES parsing failed and it looks like a formula, use formula path
                # Will be handled below
                pass
    
    # If we found a valid SMILES, use it
    if smiles is not None:
        result = _generate_chemfig(smiles, **options)
        return [result] if result else None
    
    # Last resort: try as chemical formula (e.g., 'C2H4O2')
    if _is_chemical_formula(input_string):
        selection = options.get('selection', 'first')
        n = options.get('n', None)
        enumerate_stereo = bool(options.get('enumerate_stereo', False))
        try:
            compounds = pcp.get_compounds(input_string, 'formula')
            if not compounds:
                # Return empty list so caller can safely iterate
                return []
            
            # Filter out isotopically labeled compounds (e.g., [13C], [2H])
            smiles_list = []
            seen = set()
            for c in compounds:
                if c.smiles and not re.search(r'\[\d+', c.smiles):
                    if c.smiles not in seen:
                        smiles_list.append(c.smiles)
                        seen.add(c.smiles)
            if not smiles_list:
                return []
            
            # Optional: reduce to constitutional isomers only (ignore stereochem & isotopes)
            if selection == 'all':
                unique_map = {}
                for smi in smiles_list:
                    mol_u = Chem.MolFromSmiles(smi)
                    if mol_u is None:
                        continue

                    # Skip multi-fragment or charged species when looking for structural isomers.
                    if '.' in Chem.MolToSmiles(mol_u) or len(Chem.GetMolFrags(mol_u)) > 1:
                        continue
                    if rdmolops.GetFormalCharge(mol_u) != 0:
                        continue
                    if any(atom.GetFormalCharge() != 0 for atom in mol_u.GetAtoms()):
                        continue

                    # Remove isotopes
                    for a in mol_u.GetAtoms():
                        a.SetIsotope(0)
                    # Remove stereochemistry
                    Chem.RemoveStereochemistry(mol_u)
                    base = Chem.MolToSmiles(mol_u, isomericSmiles=False)
                    if base not in unique_map:
                        unique_map[base] = smi  # keep first original stereochemical representative
                smiles_list_reduced = list(unique_map.values())
            else:
                smiles_list_reduced = smiles_list

            if selection == 'first':
                selected_smiles = [smiles_list_reduced[0]]
            elif selection == 'random':
                selected_smiles = [random.choice(smiles_list_reduced)]
            elif selection == 'all':
                selected_smiles = smiles_list_reduced
            elif selection == 'first_n':
                if n is None:
                    n = 1
                selected_smiles = smiles_list_reduced[:n]
            else:
                return []

            if enumerate_stereo and selected_smiles:
                expanded = []
                seen_expanded = set()
                for base_smiles in selected_smiles:
                    # Enumerate locally to avoid downloading every stereoisomer from PubChem.
                    for stereo_smiles in _enumerate_stereo_via_rdkit(base_smiles):
                        if stereo_smiles not in seen_expanded:
                            expanded.append(stereo_smiles)
                            seen_expanded.add(stereo_smiles)
                if expanded:
                    selected_smiles = expanded
            
            generated = [_generate_chemfig(smiles, **options) for smiles in selected_smiles]
            generated = [g for g in generated if g is not None]
            return generated
        except Exception:
            return []
    else:
        # Try as PubChem name if nothing else matched
        try:
            compounds = pcp.get_compounds(input_str, 'name')
            if compounds and compounds[0].smiles:
                smiles = compounds[0].smiles
                result = _generate_chemfig(smiles, **options)
                return [result] if result else None
        except Exception:
            pass
    # Nothing worked
    return []





