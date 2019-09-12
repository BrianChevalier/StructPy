# /bin/bash
clear

# This runs all available unit tests
printf "====== Running Unit Tests ======\n\n"
python3 -W ignore::DeprecationWarning -m unittest discover -s 'Unit Tests' -p '*_tests.py'

# This runs all available doctests
printf "\n\n\n\n====== Running Doc Tests ======\n\n"
python3 -m doctest StructPy/*.py