#!/usr/bin/env bash

# Preamble
test_path=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
test_data_path="$test_path/"$(basename "${BASH_SOURCE[0]}" | sed -e 's@\.t$@@')
source $test_path/setup.sh
freq_in=$test_data_path/data_files

plan 5

ok_ 'Compare progeny to full family, input catalog' \
    000_sib-fam_incat \
    "sstacks -c $freq_in/batch_2 -s $freq_in/progeny_002 -o %out"

ok_ 'Compare progeny to single parent, input single sample TSV files' \
    001_sib-f0m_infile \
    "sstacks -r $freq_in/f0_male -s $freq_in/progeny_002 -o %out"

ok_ 'Compare individual progeny to parents, similar to genotypes, input catalog' \
    002_sib-f0_incat \
    "sstacks -c $freq_in/batch_1 -s $freq_in/progeny_002 -o %out"

skip_ 'Use genomic location as the basis for matching' \
    003_base_gen \
    ""

skip_ 'Do not verify haplotype of matching locus' \
    004_nohap \
    ""

finish