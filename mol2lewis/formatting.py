"""
Formatting Module
=================

Handles formatting of chemfig code for human readability.
"""


def format_chemfig(code):
    """
    Format chemfig code for human readability with proper indentation.
    
    mol2chemfig already provides basic structure with newlines and indentation.
    We enhance it by:
    - Normalizing indentation (4 spaces per level)
    - Ensuring one bond-atom per line
    - Keeping \\charge{} with its atom on the same line
    - Clear branch structure
    
    Args:
        code (str): Raw chemfig code string
        
    Returns:
        str: Formatted chemfig code with proper indentation
    """
    lines = code.split('\n')
    result = []
    indent_level = 0
    indent_str = '    '  # 4 spaces per level
    
    for line in lines:
        # Strip existing indentation
        stripped = line.strip()
        if not stripped:
            continue
        
        # Handle opening parenthesis
        if stripped == '(':
            result.append(indent_str * indent_level + stripped)
            indent_level += 1
            continue
        
        # Handle closing parenthesis
        if stripped == ')':
            indent_level = max(0, indent_level - 1)
            result.append(indent_str * indent_level + stripped)
            continue
        
        # Check if line ends with opening paren (e.g., "atom(")
        if stripped.endswith('(') and not stripped.startswith('('):
            # Split the line at the opening paren
            main_part = stripped[:-1].strip()
            if main_part:
                result.append(indent_str * indent_level + main_part)
            result.append(indent_str * indent_level + '(')
            indent_level += 1
            continue
        
        # Check if line ends with closing paren (e.g., "atom)")
        if stripped.endswith(')') and not stripped.startswith(')'):
            # Split the line at the closing paren
            main_part = stripped[:-1].strip()
            indent_level = max(0, indent_level - 1)
            if main_part:
                result.append(indent_str * (indent_level + 1) + main_part)
            result.append(indent_str * indent_level + ')')
            continue
        
        # Regular line - just add with current indentation
        result.append(indent_str * indent_level + stripped)
    
    return '\n'.join(result)
