#!/bin/bash
# Before running this script you must build with scons TYPE=coverage
set -e

# This script assumes that either it is run from the ziran base directory
# or that the ZIRAN variable is set to the location of the base directory
ZIRAN="${ZIRAN:-$(pwd -P)}"

# reset execution counts to zero
lcov --directory "$ZIRAN/Tests" --zerocounters
# capture baseline data (so that unrun functions report zero coverage)
# uninstantiated templates will not be included
lcov --initial --capture --directory "$ZIRAN/Tests" --output-file coverage.initial
# run the tests
./Tests/tests
# capture run counts
lcov --capture --directory "$ZIRAN/Tests" --output-file coverage.run
# combine initial and run
lcov --add-tracefile coverage.initial --add-tracefile coverage.run --output-file coverage.total
# remove run counts for non-ziran files
lcov --extract coverage.total "$ZIRAN/*" --output-file coverage.ziran
# generate html 
genhtml coverage.ziran --output-directory coverage
rm coverage.initial coverage.run coverage.total coverage.ziran
