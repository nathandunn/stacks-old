#!/usr/bin/env bash

# Preamble
test_path=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
test_data_path="$test_path/"$(basename "${BASH_SOURCE[0]}" | sed -e 's@\.t$@@')
source $test_path/setup.sh

# Setup
barcodes=$test_data_path/frequent_data/Barcodes.txt
freq_in=$test_data_path/frequent_data/in.fastq.gz 

plan 13

# # Example libtap tests.  Uncomment to run.
# ok "This test will pass" true
# ok "This test will fail" ls -al /this/file/does/not/exist
# diag 'I just love word plays ...'
# ok "This test is expected to fail# TODO fix this" ls -al /neither/does/this/file
# skip "This command doesn't make sense:" more /dev/null

# process_radtags tests
ok_ 'example failure # TODO: This is expected to fail' \
    000_example_failing_test \
    "process_radtags -i fastq -p %in -o %out -b $barcodes -E phred33 -e sbfI -r -c -q"

ok_ 'input gzfastq' \
    001_input_gzfastq \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI"

ok_ 'input fastq' \
    002_input_fastq \
    "process_radtags -i fastq -p %in -o %out -b $barcodes -E phred33 -e sbfI"

skip_ 'input bustard - not yet implemented' \
    003_input_bustard \
    "process_radtags -i bustard -p %in -o %out -b $barcodes -E phred33 -e sbfI"

diag 'FIXME: Input files for this test are NOT actaully phred64 encoded! This is just an example test...'

ok_ 'input phred64' \
    004_input_phred64 \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred64 -e sbfI"

ok_ 'clean' \
    005_clean_data \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -c"

ok_ 'discarded reads' \
    006_discarded_reads \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -D"

ok_ 'fasta output' \
    007_output_fasta \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -y fasta"

ok_ 'discard low quality reads' \
    008_discard_lq \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -q"

ok_ 'rescue barcodes and radtags' \
    009_rescue_bcrt \
    "process_radtags -i gzfastq -p %in -o %out -b $barcodes -E phred33 -e sbfI -r"

ok_ 'truncate final read length' \
    010_truncate \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -t 50"

ok_ 'set window length' \
    011_winlen \
    "process_radtags -i gzfastq -p %in -o %out -b $barcodes -E phred33 -e sbfI -q -w .12"

ok_ 'minimum window score' \
    012_minscore \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -q -s 15"

# I'm not sure yet what finish() does.
finish
