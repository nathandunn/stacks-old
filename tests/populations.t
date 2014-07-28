#!/usr/bin/env bash

# Preamble
test_path=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
test_data_path="$test_path/"$(basename "${BASH_SOURCE[0]}" | sed -e 's@\.t$@@')
source $test_path/setup.sh
data_files=$test_data_path/data_files/

plan 29

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

ok_ -i '%in' \
    'Minimum percent of population a loci must be present in to report' \
    005_minper \
    "populations -P %in/ -b 1 -M $data_files/map.txt -r 0.5"

ok_ -i '%in' \
    'Minimum number of population a loci must be present in to report' \
    006_minnum \
    "populations -P %in/ -b 1 -M $data_files/map.txt -p 3"

ok_ -i $data_files \
    'minimum stack depth to report a loci' \
    007_mindep \
    "populations -P %in/ -b 1 -m 10 -M $data_files/map.txt"

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

ok_ -i $data_files \
    'specify p-value required to keep Fst' \
    012_pcut \
    "populations -P %in -b 1 -M $data_files/map.txt --p_value_cutoff 0.03 -f p_value"

ok_ -i $data_files \
    'create genomic output file' \
    013_genomic \
    "populations -P %in -b 1 -M $data_files/map.txt --genomic"

ok_ -i $data_files \
    'create fasta output file' \
    014_fasta \
    "populations -P %in -b 1 -M $data_files/map.txt --fasta"

skip_ -i $data_files \
    'VCF format output' \
    015_vcf \
    "populations -P %in -b 1 -M $data_files/map.txt --vcf"

skip_ -i $data_files \
    'GenePop format output' \
    016_genepop \
    "populations -P %in -b 1 -M $data_files/map.txt --genepop"

skip_ -i $data_files \
    'Structure format output' \
    017_structure \
    "populations -P %in -b 1 -M $data_files/map.txt --structure"

ok_ -i $data_files \
    'PHASE/fastPHASE format output' \
    018_phase \
    "populations -P %in -b 1 -M $data_files/map.txt --phase"

skip_ -i $data_files \
    'Beagle output format' \
    019_beagle \
    "populations -P %in -b 1 -M $data_files/map.txt --beagle"

skip_ -i $data_files \
    'Plink output format' \
    020_plink \
    "populations -P %in -b 1 -M $data_files/map.txt --plink"

skip_ -i '%in' \
    'Phylip output format' \
    021_phylip \
    "populations -P %in -b 1 -M $data_files/map.txt --phylip"

skip_ -i $data_files \
    'Include variable sites in phylip output' \
    022_phylip_var \
    "populations -P %in -b 1 -M $data_files/map.txt --phylip --phylip_var"

skip_ -i '%in' \
    'Write only the first snp per locus in GenePop output' \
    023_wss \
    "populations -P %in -b 1 -M $data_files/map.txt --structure --write_single_snp"

ok_ -i $data_files \
    'Enable kernel smoothing' \
    024_ks \
    "populations -P %in -b 1 -M $data_files/map.txt -k"

skip_ -i $data_files \
    'set window size for smoothing' \
    025_winsize \
    "populations -P %in -b 1 -M $data_files/map.txt -k"

ok_ -i $data_files \
    'Bootstrap 100 reps' \
    026_btstrp \
    "populations -P %in -b 1 -M $data_files/map.txt --bootstrap"

ok_ -i $data_files \
    'Bootstrap 10 reps' \
    027_btstrp_reps \
    "populations -P %in -b 1 -M $data_files/map.txt --bootstrap --bootstrap_reps 10"

ok_ -i $data_files \
    'Bootstrap for whitelisted loci only' \
    028_btstrp_wl \
    "populations -P %in -b 1 -M $data_files/map.txt --bootstrap_wl $data_files/whitelist.txt -k"

finish