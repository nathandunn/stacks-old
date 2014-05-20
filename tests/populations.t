#!/usr/bin/env bash

# Preamble
test_path=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
test_data_path="$test_path/"$(basename "${BASH_SOURCE[0]}" | sed -e 's@\.t$@@')
source $test_path/setup.sh
data_files=$test_data_path/data_files/

plan 12

ok_ -i $data_files \
    'Basic run using batch num 1' \
    000_basic \
    "populations -P %in -b 1"

ok_ -i $data_files \
    'Use a population map' \
    001_popmap \
    "populations -P %in -b 1 -M $data_files/map.txt"

ok_ -i $data_files \
    'Create file to input results into MySQL' \
    002_insql \
    "populations -P %in -b 1 -s -M $data_files/map.txt"

ok_ -i $data_files \
    'Black list loci 5' \
    003_black \
    "populations -P %in -b 1 -B $data_files/blacklist.txt -M $data_files/map.txt"

ok_ -i $data_files \
    'White list loci 1,2,4,5' \
    004_white \
    "populations -P %in -b 1 -W $data_files/whitelist.txt -M $data_files/map.txt"

skip_ -i $data_files \
    'Minimum percent of population a loci must be present in to report' \
    005_minperc \
    "populations -P %in -b 1 -M $data_files/map.txt"

skip_ -i $data_files \
    'Minimum number of population a loci must be present in to report' \
    006_minnum \
    "populations -P %in -b 1 -M $data_files/map.txt"

ok_ -i $data_files \
    'minimum stack depth to report a loci' \
    007_mindep \
    "populations -P %in -b 1 -m 10 -M $data_files/map.txt"

ok_ -i $data_files \
    'Minimum minor allele frequency required before calculating Fst' \
    008_qfreq \
    "populations -P %in -b 1 -a 25 -M $data_files/map.txt"

ok_ -i $data_files \
    'Fst correction with p-value' \
    009_fst_p \
    "populations -P %in -b 1 -f p_value -M $data_files/map.txt"

ok_ -i $data_files \
    'Fst value correction with Bonferroni win' \
    010_fst_bw \
    "populations -P %in -b 1 -f bonferroni_win -M $data_files/map.txt"

ok_ -i $data_files \
    'Fst value correction with Bonferroni gen' \
    011_fst_bg \
    "populations -P %in -b 1 -f bonferroni_gen -M $data_files/map.txt"

#skip_ '' \
#    014_maxhet \
#    "genotypes -P $data_files -b 0 -t"

#skip_ '' \
#    015_cors \
#    "genotypes -P $data_files -b 0 -t"

finish