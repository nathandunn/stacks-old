#!/usr/bin/env bash

# Preamble
test_path=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
test_data_path="$test_path/"$(basename "${BASH_SOURCE[0]}" | sed -e 's@\.t$@@')
source $test_path/setup.sh
freq_in=$test_data_path/data_files/f0_female
freq_in2=$test_data_path/data_files/f0_male
freq_in3=$test_data_path/data_files/progeny_002
 
plan 6

ok_ 'initialize catalog' \
    000_init_cat \
    "cstacks -s $freq_in -o %out -b 1"

ok_ 'Add on to an already existing catalog' \
    001_addon_cat \
    "cstacks --catalog $test_data_path/000_init_cat/eout/batch_1.catalog -s $freq_in2 -o %out -b 1"

#Note: if either of the previous two tests fail, but not this one (below) there is an error in the catalog initiation or addition. If this one and both the previous two all fail, likely there is a problem witht the bath ID number assignment, although there may still be a problem with catalog initiation/addition.
ok_ 'Set sql ID/sample ID to 1' \
    002_sqlid \
    "cstacks -s $freq_in -b 2 -o %out"

ok_ 'Include tags that match multiple entries' \
    003_mm \
    "cstacks -s $freq_in -s $freq_in2 -s $freq_in3 -m -o %out -b 1"

ok_ 'Report tags that match multiple entries' \
    004_mmrpt \
    "cstacks -s $freq_in -s $freq_in2 -s $freq_in3 --report_mmatches -o %out -b 1"

ok_ 'Allow 1 mismatched base between sample tags' \
    005_num_mis \
    "cstacks --catalog $test_data_path/000_init_cat/eout/batch_1.catalog -s $freq_in2 -o %out -n 1 -b 1"

finish