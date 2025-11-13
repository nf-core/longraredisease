#!/usr/bin/env python3
import sys
import os

def fix_nextflow_formatting(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    fixed_lines = []
    for line in lines:
        # Remove trailing whitespace and add newline
        cleaned = line.rstrip()
        if cleaned:  # Non-empty line
            # Count leading spaces
            leading_spaces = len(line) - len(line.lstrip())
            if leading_spaces > 0:
                # Round to nearest multiple of 4
                new_indent = ((leading_spaces + 2) // 4) * 4
                fixed_line = ' ' * new_indent + line.lstrip().rstrip() + '\n'
            else:
                fixed_line = cleaned + '\n'
        else:  # Empty line
            fixed_line = '\n'

        fixed_lines.append(fixed_line)

    # Ensure file ends with newline
    if fixed_lines and not fixed_lines[-1].endswith('\n'):
        fixed_lines[-1] += '\n'

    with open(filepath, 'w') as f:
        f.writelines(fixed_lines)

# Usage
if len(sys.argv) > 1:
    fix_nextflow_formatting(sys.argv[1])
