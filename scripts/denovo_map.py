#!/usr/bin/env python3

import sys, os, signal, argparse, datetime, re, subprocess
version = '_VERSION_'
install_prefix = '_INSTALLPREFIX_'.rstrip('/')
program_name = os.path.basename(__file__)

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
        log_file.write('\n{}:\n'.format(program_name))
        log_file.write('Error: Aborting because the previous command failed. ')
        log_file.write(str(process.returncode))
        if process.returncode < 0:
            sig = -process.returncode
            if sig in (signal.SIGINT, signal.SIGTERM, signal.SIGKILL):
                log_file.write('/killed')
            elif sig == signal.SIGABRT:
                log_file.write('/SIGABRT')
            elif sig == signal.SIGSEGV:
                log_file.write('/segmentation fault')
        log_file.write('\n')
        sys.exit('{}: Aborted, see log file for details.'.format(program_name))

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
parser.add_argument('--quiet', action='store_true')

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

'''.format(prog=program_name, version=version)
parser.format_help = parser.format_usage

def main(args):
    # Post-process the parsed arguments.
    # ==========
    args.samples = args.samples.rstrip('/')
    x_options = args.X
    args.X = {
        'ustacks': [], 'cstacks': [], 'sstacks': [],
        'tsv2bam': [], 'gstacks': [], 'populations': [] }
    if x_options is not None:
        for x_option in x_options:
            if not re.search(r'^ *\w+ *:', x_option):
                sys.exit('Error: Illegal -X option value \'{}\''.format(x_option))
            program, arguments = x_option.split(':', 1)
            program = program.strip()
            arguments = arguments.split()
            if not program in args.X:
                sys.exit('Error: Unknown program \'{}\' in -X\'{}\'.'
                    .format(program, x_option))
            args.X[program] += arguments
    # Load the samples names from the population map.
    # ==========
    samples_names = []
    with open(args.popmap, 'r') as popmap_file:
        for i, line in enumerate(popmap_file):
            if ' ' in line:
                print('WARNING: The population map contains spaces.', file=sys.stderr)
            fields = line.rstrip('\n').split('\t')
            if not 2 <= len(fields) <= 3:
                sys.exit('ERROR: Bad number of fields at line {} of the'
                    ' population map, \'{}\''
                    .format(i+1, line.rstrip('\n')))
            samples_names.append(fields[0])
    # Find the reads files.
    # ==========
    known_extensions = [
        '.1{}{}'.format(ext, gz)
        for ext in ['.fq', '.fastq', '.fa', '.fasta'] for gz in ['.gz', '']]
    if not args.paired:
        known_extensions += [e[2:] for e in known_extensions]
    extension = None
    first_sample = samples_names[0]
    filenames = os.listdir(args.samples)
    filenames = set([f for f in filenames if f.startswith(first_sample)])
    if not filenames:
        sys.exit('Error: No files matching \'{}/{}.*\''.format(args.samples, first_sample))
    for ext in known_extensions:
        first_sample_filename = first_sample + ext
        if first_sample_filename in filenames:
            extension = ext
            break
    if extension is None:
        sys.exit('Error: Could not find read files in \'{}/\'; '
            'expected \'{}{}.fq|fastq|fa|fasta(.gz)\'.'
            .format(args.samples, first_sample, '.1' if args.paired else '(.1)'))
    if args.paired:
        reads2_path = '{}/{}.2.{}'.format(args.samples, first_sample, extension[3:])
        if not os.path.exists(reads2_path):
            sys.exit('Error: Found \'{}/{}\' but \'{}\' does not exist (and --paired '
                'was specified).'
                .format(args.samples, first_sample + extension, reads2_path))
    # Build all commands.
    # ==========
    # ustacks
    ustacks = {}
    for sample_index, sample_name in enumerate(samples_names):
        input_file_path = '{}/{}{}'.format(args.samples, sample_name, extension)
        cmd = [
            '{}/bin/ustacks'.format(install_prefix),
            '-f', input_file_path,
            '-o', args.outdir,
            '-i', str(sample_index + 1),
            '--name', sample_name]
        if args.ustacks_M is not None: cmd += ['-M', str(args.ustacks_M)]
        if args.threads is not None: cmd += ['-p', str(args.threads)]
        cmd += args.X['ustacks']
        ustacks[sample_name] = cmd
    # cstacks
    cstacks = [
        '{}/bin/cstacks'.format(install_prefix),
        '-P', args.outdir,
        '-M', args.popmap]
    if args.cstacks_n is not None: cstacks += ['-n', str(args.cstacks_n)]
    if args.threads is not None: cstacks += ['-p', str(args.threads)]
    cstacks += args.X['cstacks']
    # sstacks
    sstacks = [
        '{}/bin/sstacks'.format(install_prefix),
        '-P', args.outdir,
        '-M', args.popmap]
    if args.threads is not None: sstacks += ['-p', str(args.threads)]
    sstacks += args.X['sstacks']
    # tsv2bam
    tsv2bam = [
        '{}/bin/tsv2bam'.format(install_prefix),
        '-P', args.outdir,
        '-M', args.popmap]
    if args.threads is not None: tsv2bam += ['-t', str(args.threads)]
    if args.paired: tsv2bam += ['-R', args.samples]
    tsv2bam += args.X['tsv2bam']
    # gstacks
    gstacks = [
        '{}/bin/gstacks'.format(install_prefix),
        '-P', args.outdir,
        '-M', args.popmap]
    if args.threads is not None: gstacks += ['-t', str(args.threads)]
    gstacks += args.X['gstacks']
    # populations
    populations = [
        '{}/bin/populations'.format(install_prefix),
        '-P', args.outdir,
        '-M', args.popmap]
    if args.threads is not None: populations += ['-t', str(args.threads)]
    populations += args.X['populations']
    # Open the log file & run the pipeline.
    # ==========
    os.listdir(args.outdir) # (Check for a directory.)
    log_file = open('{}/denovo_map.log'.format(args.outdir), 'w')
    if not args.quiet:
        print('Logging to \'{}/denovo_map.log\'...'.format(args.outdir))
    try:
        log_file.write('denovo_map.py version {} started at {}\n'
            .format(version, get_current_time()))
        log_file.write(' '.join(
            ['\'{}\''.format(word) if ' ' in word else word for word in sys.argv]))
        log_file.write('\n')
        # Run all commands.
        log_file.write('\nustacks\n==========\n')
        for sample_index, sample_name in enumerate(samples_names):
            log_file.write(
                "\nsample {} of {} '{}'\n----------\n"
                .format(sample_index + 1, len(samples_names), sample_name))
            run_command(ustacks[sample_name], log_file, args)
        log_file.write('\ncstacks\n==========\n')
        run_command(cstacks, log_file, args)
        log_file.write('\nsstacks\n==========\n')
        run_command(sstacks, log_file, args)
        log_file.write('\ntsv2bam\n==========\n')
        run_command(tsv2bam, log_file, args)
        log_file.write('\ngstacks\n==========\n')
        run_command(gstacks, log_file, args)
        log_file.write('\npopulations\n==========\n')
        run_command(populations, log_file, args)
        log_file.write('\ndenovo_map.py is done.\n')
        log_file.write('\ndenovo_map.py completed at {}\n'.format(get_current_time()))
    except Exception as e:
        print(type(e).__name__, ': ', e, '\nAborted.', sep='', file=log_file)
        raise

if __name__ == '__main__':
    args = parser.parse_args()
    try:
        main(args)
    except (FileNotFoundError, PermissionError, NotADirectoryError) as e:
            print(e, file=sys.stderr)
            sys.exit(1)
