"""
Geometry Module
===============

Handles angle calculations and lone pair placement based on 2D molecular geometry.
Contains all the geometric logic for determining where to place lone pairs
on heteroatoms based on bond angles from RDKit 2D coordinates.
"""

import math


def _norm_deg(a):
    """Normalize angle to [0, 360) range."""
    return (a + 360) % 360


def _min_angular_distance(deg, others):
    """Return minimal cyclic angular distance between deg and any in others."""
    if not others:
        return 180
    dists = [min(abs(_norm_deg(deg - o)), 360 - abs(_norm_deg(deg - o))) for o in others]
    return min(dists)


def _choose_max_spaced_candidates(required, candidates, occupied):
    """
    Choose `required` candidate angles from `candidates` so that the minimum 
    angular distance to `occupied` angles + already chosen ones is maximized.
    Uses greedy selection.
    """
    chosen = []
    avail = [c for c in candidates if c not in occupied]
    for _ in range(required):
        best = None
        best_min = -1
        for c in avail:
            others = occupied + chosen
            m = _min_angular_distance(c, others)
            if m > best_min:
                best_min = m
                best = c
        if best is None:
            break
        chosen.append(best)
        avail.remove(best)
    return chosen


def _angle_diff(a, b):
    """Calculate the angular difference between two angles."""
    d = abs(_norm_deg(a - b))
    return d if d <= 180 else 360 - d


def _circ_mean(angles):
    """Calculate the circular mean of a list of angles."""
    if not angles:
        return 0
    vx = sum(math.cos(math.radians(a)) for a in angles)
    vy = sum(math.sin(math.radians(a)) for a in angles)
    if abs(vx) < 1e-9 and abs(vy) < 1e-9:
        return 0
    return _norm_deg(math.degrees(math.atan2(vy, vx)))


# Candidate angular slots (eight positions, every 45 degrees)
CANDIDATE_ANGLES = [i * 45 for i in range(8)]


# Lone pair and charge symbols - different for normal vs draft mode
PAIR_SYMBOL_NORMAL = r'\|'   # Normal mode: full lone pair symbol
PAIR_SYMBOL_DRAFT = r'\:'    # Draft mode: spacing marker
SINGLE_SYMBOL = r'\.'


def calculate_charge_entries(atom, conf, pt, rotation_angle=0, flip=False, flop=False, draft=False, aromatic_circles=False):
    r"""
    Calculate the \\charge{} entries for a given atom based on its geometry.
    
    Args:
        atom: RDKit atom object
        conf: RDKit conformer with 2D coordinates
        pt: RDKit periodic table
        rotation_angle: Angle in degrees to rotate all lone pairs (default: 0)
        flip: Horizontal flip - angle becomes 180-angle (default: False)
        flop: Vertical flip - angle becomes -angle (default: False)
        draft: Draft mode - annotate bonds and lone pairs with visual markers (default: False)
        aromatic_circles: If True with draft mode, aromatic bonds use only \\red (default: False)
        
    Returns:
        str: Comma-separated charge entries (e.g., "0=\\|,180=\\|")
    """
    idx = atom.GetIdx()
    sym = atom.GetSymbol()

    # Calculate bond angles around the atom from 2D coordinates
    bond_angles = []
    single_bond_angles = []
    double_bond_angles = []
    triple_bond_angles = []
    aromatic_bond_angles = []
    bond_angle_to_order = {}  # Map angle to bond order
    pos0 = conf.GetAtomPosition(idx)
    bond_order_sum = 0

    for b in atom.GetBonds():
        nbr = b.GetOtherAtom(atom)
        posn = conf.GetAtomPosition(nbr.GetIdx())
        dx = posn.x - pos0.x
        dy = posn.y - pos0.y
        ang = math.degrees(math.atan2(dy, dx))
        ang_n = _norm_deg(ang)
        bond_angles.append(ang_n)
        
        bo = int(round(b.GetBondTypeAsDouble()))
        bond_angle_to_order[ang_n] = bo
        
        bond_order_sum += bo
        
        if b.GetIsAromatic():
            aromatic_bond_angles.append(ang_n)
        
        if bo >= 3:
            triple_bond_angles.append(ang_n)
        elif bo >= 2:
            double_bond_angles.append(ang_n)
        elif bo == 1:
            single_bond_angles.append(ang_n)
    
    # DRAFT MODE: Add markers at bond connection points for ALL atoms
    if draft:
        # Apply transformations to bond angles
        transformed_angles = bond_angles[:]
        if flip:
            transformed_angles = [_norm_deg(180 - a) for a in transformed_angles]
        if flop:
            transformed_angles = [_norm_deg(-a) for a in transformed_angles]
        if rotation_angle != 0:
            transformed_angles = [_norm_deg(a + rotation_angle) for a in transformed_angles]
        
        # Build entries based on bond order
        entries = []
        for i, orig_angle in enumerate(bond_angles):
            trans_angle = transformed_angles[i]
            bo = bond_angle_to_order[orig_angle]
            angle_int = int(round(trans_angle))
            
            # If aromatic_circles is True, use only \red for aromatic bonds
            is_aromatic = aromatic_bond_angles and orig_angle in aromatic_bond_angles
            
            if aromatic_circles and is_aromatic:
                # Aromatic bonds: only \red (pas d'alternance)
                entries.append(f"{angle_int}=\\red")
            elif bo >= 3:  # Triple bond: \red then \redd
                entries.append(f"{angle_int}=\\red")
                entries.append(f"{angle_int}=\\redd")
            elif bo >= 2:  # Double bond: \redd
                entries.append(f"{angle_int}=\\redd")
            else:  # Single bond: \red
                entries.append(f"{angle_int}=\\red")
        
        # In draft mode, still need to calculate and add lone pairs with \: marker
        # Don't return yet, continue to calculate lone pairs below
        # But for H atoms in draft mode, return bond markers only
        if sym == 'H':
            return ','.join(entries)
    
    # NORMAL MODE: Only process heteroatoms for lone pairs
    if sym == 'H':
        return ''

    # Valence electrons based count
    valence_e = pt.GetNOuterElecs(atom.GetAtomicNum())
    lone_elec = max(0, int(round(valence_e - bond_order_sum - atom.GetFormalCharge())))
    lone_pairs = lone_elec // 2

    # Cap pairs by element
    max_lp_by_atom = {
        'H': 0, 'B': 1, 'C': 0, 'N': 2, 'O': 2, 'S': 3,
        'P': 3, 'F': 3, 'Cl': 3, 'Br': 3, 'I': 3,
    }
    lone_pairs = min(lone_pairs, max_lp_by_atom.get(sym, 3))

    occupied = list(bond_angles)
    lp_angles = []

    # Case 1: oxygen, sulfur, selenium with two singles and two LPs
    if sym in ['O', 'S', 'Se'] and lone_pairs >= 2 and len(double_bond_angles) == 0 and len(single_bond_angles) == 2:
        a, b = single_bond_angles[0], single_bond_angles[1]
        sep = _angle_diff(a, b)
        bis = _circ_mean([a, b])
        opp = _norm_deg(bis + 180)
        half_sep = sep / 2.0
        lp_angles = [_norm_deg(opp - half_sep), _norm_deg(opp + half_sep)]

    # Case 2: at least one double bond
    elif double_bond_angles and lone_pairs > 0:
        # Special case: element with 5 valence electrons (N, P, As) with 1 double + 1 single bond
        if sym in ['N', 'P', 'As'] and len(double_bond_angles) == 1 and len(single_bond_angles) == 1:
            # Place lone pair opposite to the middle between double and single bond
            double_angle = double_bond_angles[0]
            single_angle = single_bond_angles[0]
            mid_angle = _circ_mean([double_angle, single_angle])
            lp_angle = _norm_deg(mid_angle + 180)
            lp_angles = [lp_angle]
        else:
            # General case: place opposite to double bond
            avg_double = _circ_mean(double_bond_angles)
            opp = _norm_deg(avg_double + 180)
            if lone_pairs >= 2:
                lp_angles = [_norm_deg(opp - 60), _norm_deg(opp + 60)]
            else:
                lp_angles = [opp]

    # Case 3: Nitrogen with 3 singles and 1 LP: place opposite first bond (main branch)
    elif sym == 'N' and lone_pairs == 1 and len(single_bond_angles) >= 2:
        # Place LP opposite the first bond (main branch connection)
        first_bond = single_bond_angles[0]
        lp_angle = _norm_deg(first_bond + 180)
        lp_angles = [lp_angle]

    # Fallback
    else:
        lp_angles = _choose_max_spaced_candidates(lone_pairs, CANDIDATE_ANGLES, occupied)

    # Apply transformations to all lone pair angles
    # Order: flip -> flop -> rotation
    if flip:
        lp_angles = [_norm_deg(180 - a) for a in lp_angles]
    if flop:
        lp_angles = [_norm_deg(-a) for a in lp_angles]
    if rotation_angle != 0:
        lp_angles = [_norm_deg(a + rotation_angle) for a in lp_angles]

    # Build entries for lone pairs - different symbols for normal vs draft
    pair_symbol = PAIR_SYMBOL_DRAFT if draft else PAIR_SYMBOL_NORMAL
    lp_entries = [f"{int(round(a))}={pair_symbol}" for a in lp_angles]

    fc = atom.GetFormalCharge()
    if fc != 0:
        sign = '+' if fc > 0 else '-'
        lp_entries.append(f"45={sign}{'' if abs(fc)==1 else abs(fc)}")

    # In draft mode, combine bond markers with lone pair markers
    if draft:
        # 'entries' already has bond markers from above
        all_entries = entries + lp_entries
        return ','.join(all_entries)
    else:
        # Normal mode: only return lone pair entries
        return ','.join(lp_entries)
