#!/bin/bash

echo "Fixing line endings and final newlines..."
find . -name "*.json" -o -name "*.txt" -o -name "*.bed" -o -name "*.html" -o -name "*.nf.test" | while read file; do
    # Fix line endings
    sed -i 's/\r$//' "$file"
    # Add final newline if missing
    tail -c1 "$file" | read -r _ || echo >> "$file"
done

echo "Removing trailing whitespace..."
find . -name "*.nf.test" -o -name "*.html" -o -name "*.txt" -o -name "*.config" | xargs sed -i 's/[[:space:]]*$//'

echo "Converting tabs to spaces..."
sed -i 's/\t/    /g' modules/nf-core/survivor/merge/tests/main.nf.test

echo "Fixing basic indentation issues..."
# Fix 2-space indentation to 4-space
find . -name "*.html" -o -name "*.txt" | xargs sed -i 's/^  /    /g'

echo "Done! Now commit the changes."
