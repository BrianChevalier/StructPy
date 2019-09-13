#!/bin/bash
clear

rm Unit\ Tests/*.log # Remove old logs.

# This runs all available unit tests
printf "====== Running Unit Tests ======\n\n"
#python3 -W ignore::DeprecationWarning -m unittest discover -s 'Unit Tests' -p '*_tests.py'
python3 -W ignore::DeprecationWarning -m pytest Unit\ Tests/*_pytest.py
status1=$?

# This runs all available doctests
printf "\n\n\n\n====== Running Doc Tests ======\n\n"
python3 -m doctest StructPy/*.py
status2=$?


# Exit and return status.
echo "Exiting with status code: $((status1>status2 ? status1 : status2))"
exit $((status1>status2 ? status1 : status2))
