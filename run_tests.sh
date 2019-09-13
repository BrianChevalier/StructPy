# /bin/bash
clear

# This runs all available unit tests
printf "====== Running Unit Tests ======\n\n"
python3 -W ignore::DeprecationWarning -m unittest discover -s 'Unit Tests' -p '*_tests.py'
status1=$?

# This runs all available doctests
printf "\n\n\n\n====== Running Doc Tests ======\n\n"
python3 -m doctest StructPy/*.py
status2=$?

echo "Exiting with status code: $((status1>status2 ? status1 : status2))"
exit $((status1>status2 ? status1 : status2))