# Setup source script for the Stacks Test Suite.
#
# This file is sourced by the individual test files.  It expects the following
# globals to be previously defined by the test file:
#       $test_path      - i.e. stacks_source/tests/
#       $test_data_path - e.g. stacks_source/tests/process_radtags/
# The 'preamble' section of test files define these variables.

export source_path=$(cd $test_path/.. && pwd)
export out_path=${test_data_path}_out
export PATH=$PATH:$source_path

source $test_path/libstacks.sh
