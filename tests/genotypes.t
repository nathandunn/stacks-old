#!/usr/bin/env bash

# Preamble
test_path=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
test_data_path="$test_path/"$(basename "${BASH_SOURCE[0]}" | sed -e 's@\.t$@@')
source $test_path/setup.sh
data_files=$test_data_path/data_files/

plan 9

ok_ -i $data_files \
    'Genetic cross map' \
    000_gen \
    "genotypes -P %in -b 1 -t gen"

ok_ -i $data_files \
    'Minimum progeny required to report marker' \
    001_minpro \
    "genotypes -P %in -b 1 -r 5"

skip_ 'Autocorrect results data' \
    002_autocor \
    "genotypes -P $data_files -b 1 -t"

ok_ -i $data_files \
    'CP cross map' \
    003_cp \
    "genotypes -P %in -b 1 -t cp"

ok_ -i $data_files \
    'CP cross map, output onemap' \
    004_onemap \
    "genotypes -P %in -b 1 -t cp -o onemap"

ok_ -i $data_files \
    'Set minimum stack depth for reporting' \
    005_mindep \
    "genotypes -P %in -b 1 -m 10"

ok_ -i $data_files \
    'create file for importing results to MySQL' \
    006_inmysql \
    "genotypes -P %in -b 1 -s"

ok_ -i $data_files \
    'Specify marker 1 as being blacklisted' \
    007_black \
    "genotypes -P %in -b 1 -B $data_files/blacklist.txt"

ok_ -i $data_files \
    'Specify markers 1,2,3 and 5 as being whitelisted' \
    008_white \
    "genotypes -P %in -b 1 -W $data_files/whitelist.txt"

finish