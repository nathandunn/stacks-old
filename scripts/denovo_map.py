#!/usr/bin/env python3

import sys, os, argparse, datetime, re, subprocess, unittest
version = '_VERSION_'
install_prefix = '_INSTALLPREFIX_'.rstrip('/')

# Function to run a command and save its output.
def run_command(command, log_file, args):
    assert type(command) == list
    for word in command:
        assert type(word) == str
    log_file.write(' '.join(['"{}"'.format(word) if ' ' in word else word for word in command]))
    log_file.write('\n')
    if args.dry_run:
        return
    # Run the command.
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        encoding='utf8')
    for line in process.stdout:
        log_file.write(line)
    # Check that the command succeeded.
    process.wait()
    if process.returncode != 0:
        for f in [log_file, sys.stderr]:
            print('Command failed, aborting.', file=f)
        sys.exit(1)
    return

# Get the (formatted) current time.
def get_current_time():
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    return current_time

# Create command line options for executing stacks pipeline
parser = argparse.ArgumentParser(
    description='execute stacks pipeline by running each components of the stacks individually')
parser.add_argument('-s','--samples', required=True, help='path to the directory of samples')
parser.add_argument('-p','--popmap', required=True, help='path to a population map file')
parser.add_argument('-d', '--dry-run', action='store_true')
parser.add_argument('-o', '--output', required=True, help='path to write pipline output files')
parser.add_argument('-M', '--ustacks',
    type=int, help='number of mismatches allowed between stacks within individuals (for ustacks)')
parser.add_argument('-n', '--cstacks',
    type=int,help='number of mismatches allowed between stacks between individuals (for cstacks)')
parser.add_argument('-T', '--threads', help='the number of threads/CPUs to use')
parser.add_argument('-X', '--program', action='append',
    help='additional options for specific pipeline components')
parser.add_argument('--paired', action='store_true',
    help='assemble contigs for each locus from paired-end reads')


def main():
    args = parser.parse_args()
    args.samples = args.samples.rstrip('/')

    # load command line options in a dictionary.
    dictionary_command_option = {}
    for x_options in args.program:
        x_option = x_options.split(":")
        cmd_opt_list = x_option[1].strip(' ').strip('\n').split(' ')
        dictionary_command_option[x_option[0].strip(' ')]= cmd_opt_list
    print(dictionary_command_option)
    # Create a log file, write a standard header.
    log_file = open('{}{}'.format(args.output,'/denovo_map.log'), 'w')
    log_file.write('denovo_map.py version {} started at {}\n'
        .format(version, get_current_time()))
    log_file.write(' '.join(sys.argv))
    log_file.write('\n')

    # Load the list of samples names
    list_of_samples_names = []
    with open(args.popmap, 'r') as popmap_file:
        for line in popmap_file:
            field = line.rstrip('\n').split('\t')
            sample_name = field[0]
            list_of_samples_names.append(sample_name)
    # ustacks
    # ==========
    log_file.write('\nustacks\n==========\n')
    list_of_sample_reads_paths = []
    # Find the extension of the reads files.
    known_extensions = ['.1.fa.gz', '.1.fq.gz']
    extension = None
    filenames = set(os.listdir(args.samples+'/'))
    first_sample = list_of_samples_names[0]
    for ext in known_extensions:
        first_sample_filename = first_sample + ext
        print('Trying extension {}; the file to find is {}'.format(ext, first_sample_filename))
        if first_sample_filename in filenames:
            print('File {} exists!'.format(first_sample_filename))
            extension = ext
            break
    if extension is None:
        print('Panic! I didnt find it!')
        sys.exit(1)

    # Create a list of all the input file paths.
    for sample_name in list_of_samples_names:
        list_of_sample_reads_paths.append('{}/{}{}'.format(args.samples, sample_name, extension))
    # For each sample, create and run the ustacks command.
    for sample_index, sample_name in enumerate(list_of_samples_names):
        input_file_path = list_of_sample_reads_paths[sample_index]
        log_file.write(
            "\nsample {} of {} '{}'\n----------\n"
            .format(sample_index + 1, len(list_of_samples_names), sample_name))
        ustacks_command = [
            '{}/bin/ustacks'.format(install_prefix),
            '-f', input_file_path,
            '-o', args.output,
            '-i', str(sample_index + 1),
            '--name', sample_name]
        if args.ustacks is not None:
            ustacks_command.append('-M')
            ustacks_command.append(str(args.ustacks))
        if args.threads is not None:
            ustacks_command.append('-p')
            ustacks_command.append(str(args.threads))
        if 'ustacks' in dictionary_command_option:
            ustacks_command.append(dictionary_command_option['ustacks'])
        run_command(ustacks_command, log_file, args)

    # Write command line code for cstacks, sstacks, tsv2bam, gstacks, and populations in a log file.
    # cstacks
    # ==========
    log_file.write('\ncstacks\n==========\n')
    cstacks_command = [
        '{}/bin/cstacks'.format(install_prefix),
        '-P', args.output,
        '-M', args.popmap]
    if args.cstacks is not None:
        cstacks_command.append('-n')
        cstacks_command.append(str(args.cstacks))
    if args.threads is not None:
        cstacks_command.append('-p')
        cstacks_command.append(str(args.threads))
    if 'cstacks' in dictionary_command_option:
        cstacks_command.append(dictionary_command_option['cstacks'])
    run_command(cstacks_command, log_file, args)
    # sstacks
    # ==========
    log_file.write('\nsstacks\n==========\n')
    sstacks_command = [
        '{}/bin/sstacks'.format(install_prefix),
        '-P', args.output,
        '-M', args.popmap]
    if args.threads is not None:
        sstacks_command.append('-p')
        sstacks_command.append(str(args.threads))
    if 'sstacks' in dictionary_command_option:
        sstacks_command.append(dictionary_command_option['sstacks'])
    run_command(sstacks_command, log_file, args)
    # tsv2bam
    # ==========
    log_file.write('\ntsv2bam\n==========\n')
    tsv2bam_command = [
        '{}/bin/tsv2bam'.format(install_prefix),
        '-P', args.output,
        '-M', args.popmap]
    if args.threads is not None:
        tsv2bam_command.append('-t')
        tsv2bam_command.append(str(args.threads))
    if 'tsv2bam' in dictionary_command_option:
        tsv2bam_command.append(dictionary_command_option['tsv2bam'])
    if args.paired:
        tsv2bam_command.append('-R')
        tsv2bam_command.append(args.samples+'/')
    run_command(tsv2bam_command, log_file, args)
    # gstacks
    # ==========
    log_file.write('\ngstacks\n==========\n')
    gstacks_command = [
        '{}/bin/gstacks'.format(install_prefix),
        '-P', args.output,
        '-M', args.popmap]
    if args.threads is not None:
        gstacks_command.append('-t')
        gstacks_command.append(str(args.threads))
    if 'gstacks' in dictionary_command_option:
        list_of_options = dictionary_command_option['gstacks']
        for element in list_of_options:
            gstacks_command.append(element)
    run_command(gstacks_command, log_file, args)
    # populations
    # ==========
    log_file.write('\npopulations\n==========\n')
    populations_command = [
        '{}/bin/populations'.format(install_prefix),
        '-P', args.output,
        '-M', args.popmap]
    if args.threads is not None:
        populations_command.append('-t')
        populations_command.append(str(args.threads))
    if 'populations' in dictionary_command_option:
        list_of_options = dictionary_command_option['populations']
        for element in list_of_options:
            populations_command.append(element)
    run_command(populations_command, log_file, args)

    log_file.write('\ndenovo_map.py is done.\n')
    log_file.write('\ndenovo_map.py completed at {}\n'.format(get_current_time()))

if __name__ == '__main__':
    main()
