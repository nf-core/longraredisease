#!/bin/bash
# Run the Python fixer and verify with pre-commit

echo "Running editorconfig fixer..."
echo ""

python3 fix_editorconfig.py

echo ""
echo "Verifying with pre-commit..."
echo ""

pre-commit run editorconfig-checker --all-files

exit_code=$?

echo ""
if [ $exit_code -eq 0 ]; then
    echo "✅ SUCCESS! All editorconfig issues resolved!"
else
    echo "⚠️  Some issues remain. Running fixer one more time..."
    python3 fix_editorconfig.py
    echo ""
    echo "Verifying again..."
    pre-commit run editorconfig-checker --all-files
fi
