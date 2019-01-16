#!/bin/bash

set -e
errors=0

# Run unit tests
python tcga_rppa_sig/tcga_rppa_sig_test.py || {
    echo "'python python/tcga_rppa_sig/tcga_rppa_sig_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E tcga_rppa_sig/*.py || {
    echo 'pylint -E tcga_rppa_sig/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"
