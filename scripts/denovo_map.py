#!/usr/bin/env python3

import sys, os, argparse, datetime, re, subprocess, unittest
version = '_VERSION_'
install_prefix = '_INSTALLPREFIX_'.rstrip('/')

# Function to run a command and save its output.
def run_command(command, log_file, args):
    assert type(command) == list
    for word in command:
        assert type(word) == str
    log_file.write(' '.join(
        ['\'{}\''.format(word) if ' ' in word else word for word in command]))
    log_file.write('\n')
    if args.dry_run:
        return
    if args.time_components:
        command = ['/usr/bin/time'] + command
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
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--samples', required=True)
parser.add_argument('-p','--popmap', required=True)
parser.add_argument('-o', '--outdir', required=True)
parser.add_argument('-M', '--ustacks-M', type=int)
parser.add_argument('-n', '--cstacks-n', type=int)
parser.add_argument('--paired', action='store_true')
parser.add_argument('--model', type=float)
parser.add_argument('--var-alpha', type=float)
parser.add_argument('--gt-alpha', type=float)
parser.add_argument('-X', action='append')
parser.add_argument('-T', '--threads', type=int)
parser.add_argument('-d', '--dry-run', action='store_true')
parser.add_argument('--time-components', action='store_true')

# Overwrite the help/usage behavior.
parser.format_usage = lambda : '''\
{prog} {version}
{prog} --samples dir --popmap path -o dir (assembly options) [-X "prog:opts" ...]

  Input/Output files:
    --samples: path to the directory containing the samples reads files.
    --popmap: path to a population map file (format is "<name> TAB <pop>", one sample per line).
    --outdir: path to an output directory.

  General options:
    X: additional options for specific pipeline components, e.g. -X "populations: -p 3 -r 0.50".
    T: the number of threads/CPUs to use (default: 1).
    --dry-run: Dry run. Do not actually execute anything, just print the commands
               that would be executed.

  Assembly options:
    --ustacks-M: number of mismatches allowed between alleles within individuals.
    --cstacks-n: number of mismatches allowed between alleles between individuals (suggested: set to ustacks -M).
    --paired: assemble forward reads into RAD loci, then assemble a mini-contig
              using the paired-end reads of each locus.

  SNP model options (for gstacks):
    --model: model to use for calling SNPs and genotypes.
    --var-alpha: significance level at which to call SNPs.
    --gt-alpha: significance level at which to call genotypes.

'''.format(prog=os.path.basename(__file__), version=version)
parser.format_help = parser.format_usage

def main(args):
    # Post-process the parsed arguments.
    # ==========
    args.samples = args.samples.rstrip('/')
    # load command line options in a dictionary.
    x_options = args.X
    args.X = {
        'ustacks': [], 'cstacks': [], 'sstacks': [],
        'tsv2bam': [], 'gstacks': [], 'populations': [] }
    if x_options is not None:
        for x_option in x_options:
            if not ':' in x_option:
                print('Panic!', x_option)
                sys.exit(1)
            program, arguments = x_option.split(':', 1)
            program = program.strip()
            arguments = arguments.split()
            if not program in args.X:
                print('Panic!', program)
                sys.exit(1)
            args.X[program] += arguments
    # Load the samples names from the population map.
    # ==========
    list_of_samples_names = []
    with open(args.popmap, 'r') as popmap_file:
        for line in popmap_file:
            field = line.rstrip('\n').split('\t')
            sample_name = field[0]
            list_of_samples_names.append(sample_name)
    # Find the reads files.
    # ==========
    known_extensions = ['.1.fa.gz', '.1.fq.gz']
    extension = None
    filenames = set(os.listdir(args.samples+'/'))
    first_sample = list_of_samples_names[0]
    for ext in known_extensions:
        first_sample_filename = first_sample + ext
        if first_sample_filename in filenames:
            extension = ext
            break
    if extension is None:
        print('Panic! I didnt find it!')
        sys.exit(1)
    # Open the log file, write a standard header.
    # ==========
    log_file = open('{}{}'.format(args.outdir,'/denovo_map.log'), 'w')
    log_file.write('denovo_map.py version {} started at {}\n'
        .format(version, get_current_time()))
    log_file.write(' '.join(
        ['\'{}\''.format(word) if ' ' in word else word for word in sys.argv]))
    log_file.write('\n')
    # ustacks
    # ==========
    log_file.write('\nustacks\n==========\n')
    # For each sample, create and run the ustacks command.
    for sample_index, sample_name in enumerate(list_of_samples_names):
        input_file_path = '{}/{}{}'.format(args.samples, sample_name, extension)
        log_file.write(
            "\nsample {} of {} '{}'\n----------\n"
            .format(sample_index + 1, len(list_of_samples_names), sample_name))
        ustacks_command = [
            '{}/bin/ustacks'.format(install_prefix),
            '-f', input_file_path,
            '-o', args.outdir,
            '-i', str(sample_index + 1),
            '--name', sample_name]
        if args.ustacks_M is not None:
            ustacks_command += ['-M', str(args.ustacks_M)]
        if args.threads is not None:
            ustacks_command += ['-p', str(args.threads)]
        ustacks_command += args.X['ustacks']
        run_command(ustacks_command, log_file, args)
    # cstacks
    # ==========
    log_file.write('\ncstacks\n==========\n')
    cstacks_command = [
        '{}/bin/cstacks'.format(install_prefix),
        '-P', args.outdir,
        '-M', args.popmap]
    if args.cstacks_n is not None:
        cstacks_command += ['-n', str(args.cstacks_n)]
    if args.threads is not None:
        cstacks_command += ['-p', str(args.threads)]
    cstacks_command += args.X['cstacks']
    run_command(cstacks_command, log_file, args)
    # sstacks
    # ==========
    log_file.write('\nsstacks\n==========\n')
    sstacks_command = [
        '{}/bin/sstacks'.format(install_prefix),
        '-P', args.outdir,
        '-M', args.popmap]
    if args.threads is not None:
        sstacks_command += ['-p', str(args.threads)]
    sstacks_command += args.X['sstacks']
    run_command(sstacks_command, log_file, args)
    # tsv2bam
    # ==========
    log_file.write('\ntsv2bam\n==========\n')
    tsv2bam_command = [
        '{}/bin/tsv2bam'.format(install_prefix),
        '-P', args.outdir,
        '-M', args.popmap]
    if args.threads is not None:
        tsv2bam_command += ['-t', str(args.threads)]
    if args.paired:
        tsv2bam_command += ['-R', args.samples]
    tsv2bam_command += args.X['tsv2bam']
    run_command(tsv2bam_command, log_file, args)
    # gstacks
    # ==========
    log_file.write('\ngstacks\n==========\n')
    gstacks_command = [
        '{}/bin/gstacks'.format(install_prefix),
        '-P', args.outdir,
        '-M', args.popmap]
    if args.threads is not None:
        gstacks_command += ['-t', str(args.threads)]
    gstacks_command += args.X['gstacks']
    run_command(gstacks_command, log_file, args)
    # populations
    # ==========
    log_file.write('\npopulations\n==========\n')
    populations_command = [
        '{}/bin/populations'.format(install_prefix),
        '-P', args.outdir,
        '-M', args.popmap]
    if args.threads is not None:
        populations_command += ['-t', str(args.threads)]
    populations_command += args.X['populations']
    run_command(populations_command, log_file, args)
    # Finish.
    # ==========
    log_file.write('\ndenovo_map.py is done.\n')
    log_file.write('\ndenovo_map.py completed at {}\n'.format(get_current_time()))

if __name__ == '__main__':
    args = parser.parse_args()
    try:
        sys.exit(main(args))
    except Exception as e:
        if type(e) in (FileNotFoundError, PermissionError, NotADirectoryError):
            print(e)
            sys.exit(1)
        else:
            raise
