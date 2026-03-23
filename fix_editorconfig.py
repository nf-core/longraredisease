#!/usr/bin/env python3
"""
Comprehensive fix for all editorconfig-checker issues.
Fixes: trailing whitespace, indentation (multiples of 4), and final newlines.
"""
import sys
import os

def fix_nextflow_formatting(filepath):
    """Fix formatting issues in a file according to editorconfig rules."""
    print(f"Processing: {filepath}")

    try:
        with open(filepath, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"  ✗ File not found: {filepath}")
        return False
    except Exception as e:
        print(f"  ✗ Error reading {filepath}: {e}")
        return False

    # Split into lines but keep track of original line endings
    lines = content.splitlines(keepends=True)
    if not lines:
        lines = []

    fixed_lines = []
    changes_made = False

    for line_num, line in enumerate(lines, 1):
        original_line = line

        # Remove any trailing whitespace and line ending
        cleaned = line.rstrip('\r\n\t ')

        if cleaned:  # Non-empty line
            # Count leading spaces (not tabs)
            leading_spaces = 0
            for char in line:
                if char == ' ':
                    leading_spaces += 1
                elif char == '\t':
                    # Convert tabs to spaces (assuming tab = 4 spaces)
                    leading_spaces += 4
                else:
                    break

            if leading_spaces > 0:
                # Round to nearest multiple of 4
                new_indent = ((leading_spaces + 2) // 4) * 4
                fixed_line = ' ' * new_indent + cleaned.lstrip() + '\n'

                if new_indent != leading_spaces or line != fixed_line:
                    changes_made = True
                    print(f"  Line {line_num}: indent {leading_spaces} → {new_indent}")
            else:
                fixed_line = cleaned + '\n'
                if line.rstrip('\r\n') != cleaned or not line.endswith('\n'):
                    changes_made = True
        else:  # Empty line
            fixed_line = '\n'
            if line != '\n':
                changes_made = True

        fixed_lines.append(fixed_line)

    # Ensure file ends with exactly one newline (LF)
    if not fixed_lines or not fixed_lines[-1].endswith('\n'):
        if fixed_lines:
            fixed_lines[-1] += '\n'
        else:
            fixed_lines = ['\n']
        changes_made = True
        print(f"  Added final newline")

    # Remove any extra trailing newlines (keep only one)
    while len(fixed_lines) > 1 and fixed_lines[-1] == '\n' and fixed_lines[-2].endswith('\n\n'):
        fixed_lines[-2] = fixed_lines[-2].rstrip('\n') + '\n'
        fixed_lines.pop()
        changes_made = True

    if changes_made:
        try:
            with open(filepath, 'w', newline='') as f:
                for line in fixed_lines:
                    f.write(line)
            print(f"  ✓ Fixed: {filepath}\n")
            return True
        except Exception as e:
            print(f"  ✗ Error writing {filepath}: {e}\n")
            return False
    else:
        print(f"  ○ No changes needed\n")
        return True


# List of all files with issues from your error report
files_to_fix = [
    'assets/snpeff_db.txt',
    'nextflow.config',
    'conf/modules.config',
    'modules/local/spectre/cnvcaller/main.nf',
    'subworkflows/local/call_str/tests/main.nf.test',
    'subworkflows/local/annotate_sv/main.nf',
    'subworkflows/local/align/main.nf',
    'subworkflows/local/annotsv_db/main.nf',
    'subworkflows/local/call_snv/main.nf',
    'subworkflows/local/call_str/main.nf',
    'subworkflows/local/haplotag_bam/main.nf',
    'subworkflows/local/longphase_variants/main.nf',
    'subworkflows/local/utils_nfcore_longraredisease_pipeline/main.nf',
    'subworkflows/local/merge_sv/main.nf',
    'workflows/longraredisease.nf',
    'tests/default.nf.test',
    'conf/test.config',
    'main.nf'
]

if __name__ == '__main__':
    print("=" * 70)
    print("Fixing editorconfig issues")
    print("=" * 70)
    print()

    # If specific files provided as arguments, use those
    if len(sys.argv) > 1:
        files_to_fix = sys.argv[1:]

    success_count = 0
    total_count = len(files_to_fix)

    for filepath in files_to_fix:
        if fix_nextflow_formatting(filepath):
            success_count += 1

    print("=" * 70)
    print(f"Summary: {success_count}/{total_count} files processed")
    print("=" * 70)
    print()

    if success_count == total_count:
        print("✅ All files processed successfully!")
        print()
        print("Now run: pre-commit run editorconfig-checker --all-files")
    else:
        print(f"⚠️  Some files had errors ({total_count - success_count} failed)")
    print()
