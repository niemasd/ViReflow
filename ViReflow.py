#! /usr/bin/env python3
'''
ViReflow: An elastically-scaling parallelized AWS pipeline for viral consensus sequence generation
'''

# imports
from csv import reader
from json import load as jload
from os import chdir,remove
from sys import stderr
from time import localtime, strftime
from urllib.request import urlopen
import argparse
import sys

# useful constants
VERSION = '1.0.17'
RELEASES_URL = 'https://api.github.com/repos/niemasd/ViReflow/tags'
RUN_ID_ALPHABET = set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_.')
READ_TRIMMERS = {
    'reads': {
        'fastq': {'fastp', 'prinseq', 'ptrimmer'}, # trimmers that work on raw reads (FASTQ)
        'bam':   {'ivar'},                         # trimmers that work on mapped reads (BAM)
    },
    'primers': {
        'bed':   {'ivar', 'ptrimmer'},             # trimmers that use BED primers
        'fasta': {'fastp'},                        # trimmers that use FASTA primers
    }
}
READ_TRIMMERS_ALL = {k for i in READ_TRIMMERS for j in READ_TRIMMERS[i] for k in READ_TRIMMERS[i][j]}
READ_MAPPERS = {'bowtie2', 'bwa', 'hisat2', 'minimap2'}
VARIANT_CALLERS = {'freebayes', 'ivar', 'lofreq'}
VIRSTRAIN_DBS = {'H1N1', 'HIV', 'SCOV2'}
INSTANCE_INFO = {
    'docker_image': 'niemasd/vireflow:%s' % VERSION,
    'mem':          '1*GiB',
    'disk':         '3*GiB',
}
DATE_COMMAND_BASH = 'echo "[$(date +"%Y-%m-%d %T")]"'

# defaults
DEFAULT_THREADS = 1
DEFAULT_COMPRESS = 1
DEFAULT_MIN_ALT_FREQ = 0.5
DEFAULT_READ_MAPPER = 'minimap2'
DEFAULT_READ_TRIMMER = 'ivar'
DEFAULT_VARIANT_CALLER = 'lofreq'

# help text
HELP_TEXT_CORONASPADES = "Run SPAdes in coronaSPAdes mode (optional)"
HELP_TEXT_METAVIRALSPADES = "Run SPAdes in metaviralSPAdes mode (optional)"
HELP_TEXT_MEGAHIT = "Run MEGAHIT (optional)"
HELP_TEXT_MINIA = "Run minia (optional)"
HELP_TEXT_PANGOLIN = "Run Pangolin (optional)"
HELP_TEXT_PI = "Compute pi diversity metric (optional)"
HELP_TEXT_RNAVIRALSPADES = "Run SPAdes in rnaviralSPAdes mode (optional)"
HELP_TEXT_UNICYCLER = "Run Unicycler (optional)"
HELP_TEXT_VIRSTRAIN = "Run VirStrain (optional) (select DB: %s)" % ', '.join(sorted(VIRSTRAIN_DBS))

# clear argv
def clear_argv(keep_first_arg=True):
    tmp = sys.argv[0]; sys.argv.clear()
    if keep_first_arg:
        sys.argv.append(tmp)

# convert a ViReflow version string to a tuple of integers
def parse_version(s):
    return tuple(int(v) for v in s.split('.'))

# check if string can be parsed as int
def check_int(s):
    try:
        int(s); return True
    except:
        return False

# check if string can be parsed as float
def check_float(s):
    try:
        float(s); return True
    except:
        return False

# update ViReflow to the newest version
def update_vireflow():
    tags = jload(urlopen(RELEASES_URL))
    newest = max(tags, key=lambda x: parse_version(x['name']))
    if parse_version(newest['name']) <= parse_version(VERSION):
        print("ViReflow is already at the newest version (%s)" % VERSION); exit(0)
    old_version = VERSION; new_version = newest['name']
    url = 'https://raw.githubusercontent.com/niemasd/ViReflow/%s/ViReflow.py' % newest['commit']['sha']
    content = urlopen(url).read()
    try:
        with open(__file__, 'wb') as f:
            f.write(content)
    except PermissionError:
        print("ERROR: Received a permission error when updating ViReflow. Perhaps try running as root?", file=stderr); exit(1)
    print("Successfully updated ViReflow %s --> %s" % (old_version, new_version)); exit(0)

# parse user args
def parse_args():
    # check if user wants to update ViralMSA
    if '-u' in sys.argv or '--update' in sys.argv:
        update_vireflow()

    # use argparse to parse user arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-id', '--run_id', required=False, type=str, default=strftime('%Y-%m-%d.%H-%M-%S', localtime()), help="Unique Run Identifier (for output file naming)")
    parser.add_argument('-d', '--destination', required=True, type=str, help="Destination for Results (s3 folder)")
    parser.add_argument('-rf', '--reference_fasta', required=True, type=str, help="Reference Genome Sequence (s3/http/https/ftp to FASTA)")
    parser.add_argument('-rg', '--reference_gff', required=True, type=str, help="Reference Genome Annotation (s3/http/https/ftp to GFF3)")
    parser.add_argument('-p', '--primer_bed', required=True, type=str, help="Primer (s3/http/https/ftp to BED)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Reflow File (rf)")
    parser.add_argument('-t', '--threads', required=False, type=int, default=DEFAULT_THREADS, help="Number of Threads")
    parser.add_argument('-cl', '--compression_level', required=False, type=int, default=DEFAULT_COMPRESS, help="Compression Level (1 = fastest, 9 = best)")
    parser.add_argument('--mapped_read_cap', required=False, type=int, default=None, help="Successfully-Mapped Read Cap")
    parser.add_argument('--min_alt_freq', required=False, type=float, default=DEFAULT_MIN_ALT_FREQ, help="Minimum Alt Allele Frequency for consensus sequence")
    parser.add_argument('--read_mapper', required=False, type=str, default=DEFAULT_READ_MAPPER, help="Read Mapper (options: %s)" % ', '.join(sorted(READ_MAPPERS)))
    parser.add_argument('--read_trimmer', required=False, type=str, default=DEFAULT_READ_TRIMMER, help="Read Trimmer (options: %s)" % ', '.join(sorted(READ_TRIMMERS_ALL)))
    parser.add_argument('--variant_caller', required=False, type=str, default=DEFAULT_VARIANT_CALLER, help="Variant Caller (options: %s)" % ', '.join(sorted(VARIANT_CALLERS)))
    parser.add_argument('--optional_pangolin', action='store_true', help=HELP_TEXT_PANGOLIN)
    parser.add_argument('--optional_virstrain', required=False, type=str, default=None, help=HELP_TEXT_VIRSTRAIN)
    parser.add_argument('--optional_pi_metric', action='store_true', help=HELP_TEXT_PI)
    parser.add_argument('--optional_spades_coronaspades', action='store_true', help=HELP_TEXT_CORONASPADES)
    parser.add_argument('--optional_spades_metaviralspades', action='store_true', help=HELP_TEXT_METAVIRALSPADES)
    parser.add_argument('--optional_spades_rnaviralspades', action='store_true', help=HELP_TEXT_RNAVIRALSPADES)
    parser.add_argument('--optional_unicycler', action='store_true', help=HELP_TEXT_UNICYCLER)
    parser.add_argument('--optional_megahit', action='store_true', help=HELP_TEXT_MEGAHIT)
    parser.add_argument('--optional_minia', action='store_true', help=HELP_TEXT_MINIA)
    parser.add_argument('-u', '--update', action='store_true', help="Update ViReflow (current version: %s)" % VERSION)
    parser.add_argument('fastq_files', metavar='FQ', type=str, nargs='+', help="Input FASTQ Files (s3/http/https/ftp; single biological sample)")
    args = parser.parse_args()

    # check run ID is valid
    if len(args.run_id) == 0:
        stderr.write("Run ID cannot be empty string\n"); exit(1)
    for c in args.run_id:
        if c not in RUN_ID_ALPHABET:
            stderr.write("Invalid symbol in Run ID: %s\n" % c); exit(1)

    # check number of threads is valid
    if args.threads < 1:
        stderr.write("Invalid number of threads: %s\n" % args.threads); exit(1)

    # check compression level is valid
    if args.compression_level < 1 or args.compression_level > 9:
        stderr.write("Invalid compression level: %s\n" % args.compression_level); exit(1)

    # check read cap
    if args.mapped_read_cap is not None and args.mapped_read_cap < 1:
        stderr.write("Invalid mapped read cap: %s\n" % args.mapped_read_cap); exit(1)

    # check read mapper selection is valid
    args.read_mapper = args.read_mapper.lower()
    if args.read_mapper not in READ_MAPPERS:
        stderr.write("Invalid read mapper: %s\n" % args.read_mapper); exit(1)

    # check read trimmer selection is valid
    args.read_trimmer = args.read_trimmer.lower()
    if args.read_trimmer not in READ_TRIMMERS['reads']['fastq'] and args.read_trimmer not in READ_TRIMMERS['reads']['bam']:
        stderr.write("Invalid read trimmer: %s\n" % args.read_trimmer); exit(1)

    # check variant caller selection is valid
    args.variant_caller = args.variant_caller.lower()
    if args.variant_caller not in VARIANT_CALLERS:
        stderr.write("Invalid variant caller: %s\n" % args.variant_caller); exit(1)

    # check VirStrain DB selection is valid
    if args.optional_virstrain is not None:
        args.optional_virstrain = args.optional_virstrain.upper()
        if args.optional_virstrain not in VIRSTRAIN_DBS:
            stderr.write("Invalid VirStrain DB: %s (options: %s)\n" % (args.optional_virstrain, ', '.join(sorted(VIRSTRAIN_DBS)))); exit(1)

    # user args are valid, so return
    return args

# main function
def main():
    # parse user args
    args = parse_args()
    INSTANCE_INFO['cpu'] = args.threads

    # if CSV instead of FASTQ, run main() for each CSV
    if args.fastq_files[0].lower().endswith('.csv'):
        if len(args.fastq_files) != 1:
            stderr.write("If running in CSV mode, must only provide one positional argument\n"); exit(1)
        argv_minus_csv = [s for s in sys.argv if s != args.fastq_files[0]]
        id_ind = -1; o_ind = -1
        if '-id' in argv_minus_csv:
            id_ind = argv_minus_csv.index('-id')
        elif '--run_id' in argv_minus_csv:
            id_ind = argv_minus_csv.index('--run_id')
        if id_ind != -1:
            argv_minus_csv.pop(id_ind); argv_minus_csv.pop(id_ind)
        if '-o' in argv_minus_csv:
            o_ind = argv_minus_csv.index('-o')
        elif '--output' in argv_minus_csv:
            o_ind = argv_minus_csv.index('--output')
        if o_ind != -1:
            argv_minus_csv.pop(o_ind); argv_minus_csv.pop(o_ind)
        for row in reader(open(args.fastq_files[0])):
            if len(row) < 3:
                stderr.write("Invalid ViReflow sample CSV\n"); exit(1)
            parts = [s.replace('\ufeff','').strip() for s in row]
            clear_argv(keep_first_arg=False)
            sys.argv += argv_minus_csv + ['-o', '%s.rf' % parts[0], '-id'] + parts
            main()
        exit()

    # check input files (FASTQs, ref FASTA, ref GFF, and primer BED)
    input_files = [
        ('ref_fas',    '%s.reference.fas' % args.run_id, args.reference_fasta),
        ('ref_gff',    '%s.reference.gff' % args.run_id, args.reference_gff),
        ('primer_bed', '%s.primers.bed'   % args.run_id, args.primer_bed),
    ]
    for i, f in enumerate(args.fastq_files):
        if f.lower().endswith('.gz'):
            new_fn = '%s.fq%d.fastq.gz' % (args.run_id, i+1)
        else:
            new_fn = '%s.fq%d.fastq' % (args.run_id, i+1)
        input_files.append(('fq%d' % (i+1), new_fn, f))

    # check destination
    if not args.destination.lower().startswith('s3://'):
        stderr.write("Invalid output s3 directory: %s\n" % args.destination); exit(1)
    args.destination = args.destination.rstrip('/')
    final_output_s3 = "%s/%s.tar.gz" % (args.destination, args.run_id)

    # initialize output file
    if args.output == 'stdout':
        from sys import stdout as rf_file
    else:
        rf_file = open(args.output, 'w')
    rf_file.write('// Run ID: %s\n' % args.run_id)
    rf_file.write('// Created using ViReflow (%s)\n' % VERSION)
    rf_file.write('// ViReflow command: %s\n' % ' '.join(sys.argv))
    rf_file.write('@requires(disk := %s)\n' % INSTANCE_INFO['disk'])
    rf_file.write('val Main = {\n')
    rf_file.write('    files := make("$/files")\n\n')

    # handle input files that are S3 paths
    tmp = False # check if header has been written
    for rf_var, local_fn, p in input_files:
        if p.lower().startswith('s3://'):
            if not tmp:
                rf_file.write('    // Using the following S3 input files\n'); tmp = True
            rf_file.write('    %s := file("%s")\n' % (rf_var, p))
    if tmp:
        rf_file.write('\n')

    # begin main exec
    outdir = "%s_output" % args.run_id
    rf_file.write('    // Main ViReflow pipeline execution\n')
    rf_file.write('    out_file := exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (INSTANCE_INFO['docker_image'], INSTANCE_INFO['mem'], INSTANCE_INFO['cpu']))

    # copy input files locally
    local_fns = dict()
    local_fns['vireflow_log'] = "%s/%s.vireflow.log" % (outdir, args.run_id)
    rf_file.write('        # Copy input files locally\n')
    rf_file.write('        mkdir "%s"\n' % outdir)
    rf_file.write('        %s "Starting ViReflow %s pipeline" >> %s\n' % (DATE_COMMAND_BASH, VERSION, local_fns['vireflow_log']))
    rf_file.write('        %s "Copying input files locally" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
    for rf_var, local_fn, p in input_files:
        local_fns[rf_var] = "%s/%s" % (outdir, local_fn)
        rf_file.write('        ')
        if p.startswith('s3://'):
            rf_file.write('cp "{{%s}}" "%s"\n' % (rf_var, local_fns[rf_var]))
        else:
            if p.split('://')[0].lower() not in {'http', 'https', 'ftp'}:
                if rf_var == 'ref_fas': # assume GenBank accession number
                    p = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta' % p
                else:
                    stderr.write("Invalid input file: %s\n" % p)
                    if args.output != 'stdout':
                        rf_file.close(); remove(args.output)
                    exit(1)
            rf_file.write('wget -q -O "%s" "%s"\n' % (local_fns[rf_var], p))
    rf_file.write('\n')

    # handle primer FASTA (if needed)
    if args.read_trimmer in READ_TRIMMERS['primers']['fasta']:
        local_fns['primer_fas'] = "%s/%s.primers.fas" % (outdir, args.run_id)
        rf_file.write('        # Create primer FASTA file\n')
        rf_file.write('        %s "Creating primer FASTA file" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        bedtools getfasta -fi "%s" -bed "%s" -fo "%s" 1>&2\n' % (local_fns['ref_fas'], local_fns['primer_bed'], local_fns['primer_fas']))
        rf_file.write('\n')

    # handle primer TXT (if needed, just pTrimmer for now)
    if args.read_trimmer == 'ptrimmer':
        local_fns['primer_txt'] = "%s/%s.primers.txt" % (outdir, args.run_id)
        rf_file.write('        # Create primer TXT file\n')
        rf_file.write('        %s "Creating primer TXT file" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        faidx "%s" > /dev/null\n' % local_fns['ref_fas'])
        rf_file.write('        Bed2Amplicon.py "%s" "%s" "%s" 1>&2\n' % (local_fns['ref_fas'], local_fns['primer_bed'], local_fns['primer_txt']))
        rf_file.write('\n')

    # trim unmapped reads (if using FASTQ trimmer)
    if args.read_trimmer in READ_TRIMMERS['reads']['fastq']:
        rf_file.write('        # Trim reads using %s\n' % args.read_trimmer)
        rf_file.write('        %s "Trimming reads using %s" >> %s\n' % (DATE_COMMAND_BASH, args.read_trimmer, local_fns['vireflow_log']))
        for rf_var, fq_fn, p in input_files:
            local_fn = '%s/%s' % (outdir, fq_fn)
            if not rf_var.startswith('fq'):
                continue
            trimmed_rf_var = 'trimmed_%s' % rf_var
            local_fns[trimmed_rf_var] = "%s/%s.%s.trimmed.fastq" % (outdir, args.run_id, rf_var)
            rf_file.write('        ')
            if args.read_trimmer == 'fastp':
                rf_file.write('fastp --thread %d --adapter_fasta "%s" -i "%s" -o "%s" 1>&2\n' % (args.threads, local_fns['primer_fas'], local_fn, local_fns[trimmed_rf_var]))
            elif args.read_trimmer == 'prinseq':
                rf_file.write('prinseq-lite.pl -fastq "%s" -ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10 -out_format 3 -out_good stdout -out_bad null -min_len 0 > "%s"\n' % (local_fn, local_fns[trimmed_rf_var]))
            elif args.read_trimmer == 'ptrimmer':
                rf_file.write('pTrimmer -t single -a "%s" -f "%s" -d "%s" 1>&2\n' % (local_fns['primer_txt'], local_fn, local_fns[trimmed_rf_var]))
            else:
                stderr.write("Invalid read trimmer: %s\n" % args.read_trimmer)
                if args.output != 'stdout':
                    rf_file.close(); remove(args.output)
                exit(1)
        rf_file.write('\n')

    # organize FASTQ filenames to map
    if args.read_trimmer in READ_TRIMMERS['reads']['fastq']:
        local_fq_fns_to_map = [local_fns[k] for k in local_fns if k.startswith('trimmed_fq')]
    else:
        local_fq_fns_to_map = [local_fns[k] for k in local_fns if k.startswith('fq')]

    # map reads
    local_fns['trimmed_bam'] = "%s/%s.trimmed.bam" % (outdir, args.run_id)
    rf_file.write('        # Map reads using %s' % args.read_mapper)
    if args.mapped_read_cap is not None:
        rf_file.write(' (cap at %d successfully-mapped reads)' % args.mapped_read_cap)
    rf_file.write('\n')
    rf_file.write('        %s "Mapping reads using %s' % (DATE_COMMAND_BASH, args.read_mapper))
    if args.mapped_read_cap is not None:
        rf_file.write(' (cap at %d successfully-mapped reads)' % args.mapped_read_cap)
    rf_file.write('" >> %s\n' % local_fns['vireflow_log'])
    if args.read_trimmer in READ_TRIMMERS['reads']['fastq']:
        curr_bam_var = 'trimmed_bam'
    else:
        local_fns['untrimmed_bam'] = "%s/%s.untrimmed.bam" % (outdir, args.run_id)
        curr_bam_var = 'untrimmed_bam'
    if args.read_mapper == 'bowtie2':
        rf_file.write('        bowtie2-build --threads %d -f "%s" ref 1>&2\n' % (args.threads, local_fns['ref_fas']))
        rf_file.write('        bowtie2 --threads %d -x ref -U "%s"' % (args.threads, ','.join(local_fq_fns_to_map)))
    elif args.read_mapper == 'bwa':
        rf_file.write('        bwa index "%s" 1>&2\n' % local_fns['ref_fas'])
        rf_file.write('        bwa mem -p -t %d "%s" %s' % (args.threads, local_fns['ref_fas'], ' '.join('"%s"' % fn for fn in local_fq_fns_to_map)))
    elif args.read_mapper == 'hisat2':
        rf_file.write('        hisat2-build -p %d "%s" ref 1>&2\n' % (args.threads, local_fns['ref_fas']))
        rf_file.write('        hisat2 -p %d -x ref -q "%s"' % (args.threads, ','.join(local_fq_fns_to_map)))
    elif args.read_mapper == 'minimap2':
        rf_file.write('        minimap2 -t %d -a -x sr "%s" %s' % (args.threads, local_fns['ref_fas'], ' '.join('"%s"' % fn for fn in local_fq_fns_to_map)))
    else:
        stderr.write("Invalid read mapper: %s\n" % args.read_mapper)
        if args.output != 'stdout':
            rf_file.close(); remove(args.output)
        exit(1)
    if args.mapped_read_cap is not None:
        rf_file.write(' | samhead %d successful' % args.mapped_read_cap)
    rf_file.write(' | samtools view -bS - > "%s"' % local_fns[curr_bam_var])
    if args.mapped_read_cap is not None:
        rf_file.write(' || true') # samhead truncation seems to cause non-zero exit code in read mapper
    rf_file.write('\n')
    rf_file.write('\n')

    # sort mapped reads
    local_fns['trimmed_sorted_bam'] = "%s/%s.trimmed.sorted.bam" % (outdir, args.run_id)
    rf_file.write('        # Sort mapped reads\n')
    rf_file.write('        %s "Sorting mapped reads" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
    if args.read_trimmer in READ_TRIMMERS['reads']['fastq']:
        curr_sorted_bam_var = 'trimmed_sorted_bam'
    else:
        local_fns['untrimmed_sorted_bam'] = "%s/%s.untrimmed.sorted.bam" % (outdir, args.run_id)
        curr_sorted_bam_var = 'untrimmed_sorted_bam'
    rf_file.write('        samtools sort --threads %d -O bam -o "%s" "%s" 1>&2\n' % (args.threads, local_fns[curr_sorted_bam_var], local_fns[curr_bam_var]))
    rf_file.write('\n')

    # trim mapped reads (if using BAM trimmer)
    if args.read_trimmer in READ_TRIMMERS['reads']['bam']:
        rf_file.write('        # Trim mapped reads using %s\n' % args.read_trimmer)
        rf_file.write('        %s "Trimming reads using %s" >> %s\n' % (DATE_COMMAND_BASH, args.read_trimmer, local_fns['vireflow_log']))
        if args.read_trimmer == 'ivar':
            rf_file.write('        ivar trim -x 5 -e -i "%s" -b "%s" -p "%s" 1>&2\n' % (local_fns['untrimmed_sorted_bam'], local_fns['primer_bed'], local_fns['trimmed_bam'].rstrip('.bam')))
        else:
            stderr.write("Invalid read trimmer: %s\n" % args.read_trimmer)
            if args.output != 'stdout':
                rf_file.close(); remove(args.output)
            exit(1)
        rf_file.write('\n')

    # sort trimmed BAM (if using BAM trimmer)
    if args.read_trimmer in READ_TRIMMERS['reads']['bam']:
        rf_file.write('        # Sort trimmed mapped reads\n')
        rf_file.write('        %s "Sorting trimmed mapped reads" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        samtools sort --threads %d -O bam -o "%s" "%s" 1>&2\n' % (args.threads, local_fns['trimmed_sorted_bam'], local_fns['trimmed_bam']))
        rf_file.write('\n')

    # generate pile-up from sorted trimmed BAM
    local_fns['pileup_txt'] = "%s/%s.pileup.txt" % (outdir, args.run_id)
    rf_file.write('        # Generate pile-up from sorted trimmed BAM\n')
    rf_file.write('        %s "Generating pile-up from sorted trimmed BAM" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
    rf_file.write('        samtools mpileup -A -aa -d 0 -Q 0 --reference "%s" "%s" > "%s"\n' % (local_fns['ref_fas'], local_fns['trimmed_sorted_bam'], local_fns['pileup_txt']))
    rf_file.write('\n')

    # call variants
    local_fns['variants_vcf'] = "%s/%s.variants.vcf" % (outdir, args.run_id)
    rf_file.write('        # Call variants using %s"\n' % args.variant_caller)
    rf_file.write('        %s "Calling variants using %s" >> %s\n' % (DATE_COMMAND_BASH, args.variant_caller, local_fns['vireflow_log']))
    if args.variant_caller == 'freebayes':
        rf_file.write('        freebayes --min-alternate-fraction 0.001 --pooled-continuous --ploidy 1 -f "%s" "%s" > "%s"\n' % (local_fns['ref_fas'], local_fns['trimmed_sorted_bam'], local_fns['variants_vcf']))
    elif args.variant_caller == 'ivar':
        local_fns['variants_tsv'] = "%s/%s.variants.tsv" % (outdir, args.run_id)
        rf_file.write('        cat "%s" | ivar variants -r "%s" -g "%s" -p "%s" -m 10 1>&2\n' % (local_fns['pileup_txt'], local_fns['ref_fas'], local_fns['ref_gff'], local_fns['variants_tsv']))
        rf_file.write('        ivar_variants_to_vcf.py "%s" "%s" 1>&2\n' % (local_fns['variants_tsv'], local_fns['variants_vcf']))
    elif args.variant_caller == 'lofreq':
        rf_file.write('        lofreq call -f "%s" --call-indels "%s" > "%s"\n' % (local_fns['ref_fas'], local_fns['trimmed_sorted_bam'], local_fns['variants_vcf']))
    else:
        stderr.write("Invalid variant caller: %s\n" % args.variant_caller)
        if args.output != 'stdout':
            rf_file.close(); remove(args.output)
        exit(1)
    rf_file.write('\n')

    # call depth
    local_fns['depth_txt'] = "%s/%s.depth.txt" % (outdir, args.run_id)
    rf_file.write('        # Call depth from trimmed BAM\n')
    rf_file.write('        %s "Calling depth from trimmed BAM" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
    rf_file.write('        samtools depth -d 0 -Q 0 -q 0 -aa "%s" > "%s"\n' % (local_fns['trimmed_sorted_bam'], local_fns['depth_txt']))
    rf_file.write('\n')

    # find low-depth regions
    local_fns['low_depth_tsv'] = "%s/%s.low_depth.tsv" % (outdir, args.run_id)
    rf_file.write('        # Find low-depth regions\n')
    rf_file.write('        %s "Finding low-depth regions" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
    rf_file.write('        low_depth_regions "%s" "%s" 10 1>&2\n' % (local_fns['depth_txt'], local_fns['low_depth_tsv'])) # minimum depth of 10
    rf_file.write('\n')

    # generate consensus sequence
    local_fns['consensus_fas'] = "%s/%s.consensus.fas" % (outdir, args.run_id)
    rf_file.write('        # Generate consensus sequence\n')
    rf_file.write('        %s "Generating consensus sequence" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
    rf_file.write('        alt_vars.py -i "%s" -o tmp.vcf -v %s\n' % (local_fns['variants_vcf'], args.variant_caller))
    rf_file.write('        bgzip tmp.vcf\n')
    rf_file.write('        bcftools index tmp.vcf.gz\n')
    rf_file.write('        cat "%s" | bcftools consensus -m "%s" tmp.vcf.gz > "%s"\n' % (local_fns['ref_fas'], local_fns['low_depth_tsv'], local_fns['consensus_fas']))
    rf_file.write('\n')

    # optional: run pangolin
    if args.optional_pangolin:
        local_fns['pangolin_csv'] = "%s/%s.pangolin.lineage_report.csv" % (outdir, args.run_id)
        rf_file.write('        # Run pangolin (optional)\n')
        rf_file.write('        %s "Running pangolin (optional)" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        pangolin --threads %d --outfile "%s" "%s"\n' % (args.threads, local_fns['pangolin_csv'], local_fns['consensus_fas']))
        rf_file.write('\n')

    # optional: compute pi
    if args.optional_pi_metric:
        local_fns['pi_tsv'] = "%s/%s.pi.diversity.tsv" % (outdir, args.run_id)
        rf_file.write('        # Compute pi diversity metric (optional)\n')
        rf_file.write('        %s "Computing pi diversity metric (optional)" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        pi_from_pileup "%s" 10 > "%s"\n' % (local_fns['pileup_txt'], local_fns['pi_tsv']))
        rf_file.write('\n')

    # optional: convert trimmed BAM to FASTQ
    if args.optional_virstrain is not None \
        or args.optional_spades_coronaspades \
        or args.optional_spades_metaviralspades \
        or args.optional_spades_rnaviralspades \
        or args.optional_unicycler \
        or args.optional_megahit \
        or args.optional_minia:
        local_fns['trimmed_sorted_fastq'] = 'tmp.trimmed.sorted.fastq'
        rf_file.write('        # Convert trimmed BAM to FASTQ (optional)\n')
        rf_file.write('        %s "Converting trimmed sorted BAM to FASTQ (optional)" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        bedtools bamtofastq -i "%s" -fq "%s"\n' % (local_fns['trimmed_sorted_bam'], local_fns['trimmed_sorted_fastq']))
        rf_file.write('\n')

    # optional: run VirStrain
    if args.optional_virstrain is not None:
        local_fns['virstrain_targz'] = "%s/%s.virstrain.output.%s.tar.gz" % (outdir, args.run_id, args.optional_virstrain)
        rf_file.write('        # Run VirStrain (optional) (DB: %s)\n' % args.optional_virstrain)
        rf_file.write('        %s "Running VirStrain (optional) (DB: %s)" >> %s\n' % (DATE_COMMAND_BASH, args.optional_virstrain, local_fns['vireflow_log']))
        rf_file.write('        virstrain -i "%s" -d "/usr/local/bin/VirStrain_DB/%s" -o "virstrain_out"\n' % (local_fns['trimmed_sorted_fastq'], args.optional_virstrain))
        rf_file.write('        tar -cf - "virstrain_out" | pigz -%d -p %d > "%s"\n' % (args.compression_level, args.threads, local_fns['virstrain_targz']))
        rf_file.write('\n')

    # optional: run SPAdes (coronaSPAdes mode)
    if args.optional_spades_coronaspades:
        local_fns['coronaspades_targz'] = "%s/%s.spades.coronaspades.output.tar.gz" % (outdir, args.run_id)
        rf_file.write('        # Run coronaSPAdes (optional)\n')
        rf_file.write('        %s "Running coronaSPAdes (optional)" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        coronaspades.py --threads %d -s "%s" -o "coronaspades_out"\n' % (args.threads, local_fns['trimmed_sorted_fastq']))
        rf_file.write('        tar -cf - "coronaspades_out" | pigz -%d -p %d > "%s"\n' % (args.compression_level, args.threads, local_fns['coronaspades_targz']))
        rf_file.write('\n')

    # optional: run SPAdes (metaviralSPAdes mode)
    if args.optional_spades_metaviralspades:
        local_fns['metaviralspades_targz'] = "%s/%s.spades.metaviralspades.output.tar.gz" % (outdir, args.run_id)
        rf_file.write('        # Run metaviralSPAdes (optional)\n')
        rf_file.write('        %s "Running metaviralSPAdes (optional)" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        metaviralspades.py --threads %d -s "%s" -o "metaviralspades_out"\n' % (args.threads, local_fns['trimmed_sorted_fastq']))
        rf_file.write('        tar -cf - "metaviralspades_out" | pigz -%d -p %d > "%s"\n' % (args.compression_level, args.threads, local_fns['metaviralspades_targz']))
        rf_file.write('\n')

    # optional: run SPAdes (rnaviralSPAdes mode)
    if args.optional_spades_rnaviralspades:
        local_fns['rnaviralspades_targz'] = "%s/%s.spades.rnaviralspades.output.tar.gz" % (outdir, args.run_id)
        rf_file.write('        # Run rnaviralSPAdes (optional)\n')
        rf_file.write('        %s "Running rnaviralSPAdes (optional)" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        rnaviralspades.py --threads %d -s "%s" -o "rnaviralspades_out"\n' % (args.threads, local_fns['trimmed_sorted_fastq']))
        rf_file.write('        tar -cf - "rnaviralspades_out" | pigz -%d -p %d > "%s"\n' % (args.compression_level, args.threads, local_fns['rnaviralspades_targz']))
        rf_file.write('\n')

    # optional: run Unicycler
    if args.optional_unicycler:
        local_fns['unicycler_targz'] = "%s/%s.unicycler.output.tar.gz" % (outdir, args.run_id)
        rf_file.write('        # Run Unicycler (optional)\n')
        rf_file.write('        %s "Running Unicycler (optional)" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        unicycler -t %d -s "%s" -o "unicycler_out"\n' % (args.threads, local_fns['trimmed_sorted_fastq']))
        rf_file.write('        tar -cf - "unicycler_out" | pigz -%d -p %d > "%s"\n' % (args.compression_level, args.threads, local_fns['unicycler_targz']))
        rf_file.write('\n')

    # optional: run MEGAHIT
    if args.optional_megahit:
        local_fns['megahit_targz'] = "%s/%s.megahit.output.tar.gz" % (outdir, args.run_id)
        rf_file.write('        # Run MEGAHIT (optional)\n')
        rf_file.write('        %s "Running MEGAHIT (optional)" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        megahit --num-cpu-threads %d --read "%s" --out-dir "megahit_out"\n' % (args.threads, local_fns['trimmed_sorted_fastq']))
        rf_file.write('        tar -cf - "megahit_out" | pigz -%d -p %d > "%s"\n' % (args.compression_level, args.threads, local_fns['megahit_targz']))
        rf_file.write('\n')

    # optional: run minia
    if args.optional_minia:
        local_fns['minia_targz'] = "%s/%s.minia.output.tar.gz" % (outdir, args.run_id)
        rf_file.write('        # Run minia (optional)\n')
        rf_file.write('        %s "Running minia (optional)" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
        rf_file.write('        mkdir -p "minia_out"\n')
        rf_file.write('        cd "minia_out"\n')
        rf_file.write('        minia -nb-cores %d -in "../%s"\n' % (args.threads, local_fns['trimmed_sorted_fastq']))
        rf_file.write('        cd ..\n')
        rf_file.write('        tar -cf - "minia_out" | pigz -%d -p %d > "%s"\n' % (args.compression_level, args.threads, local_fns['minia_targz']))
        rf_file.write('\n')

    # remove redundant files before compressing output
    rf_file.write('        # Remove redundant output files before compressing\n')
    rf_file.write('        %s "Removing redundant output files before compressing" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
    rf_file.write('        rm */*trimmed.bam\n') # remove unsorted BAMs
    rf_file.write('\n')

    # archive + compress output files
    rf_file.write('        # Compress output files\n')
    rf_file.write('        %s "Compressing output files" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
    rf_file.write('        tar cf - "%s" | pigz -%d -p %d > "{{out}}"\n' % (outdir, args.compression_level, args.threads))

    # end main exec, copy output file, and end run file
    rf_file.write('    "}\n')
    rf_file.write('    cp_out_file := files.Copy(out_file, "%s")\n' % final_output_s3)
    rf_file.write('    (cp_out_file)\n')
    rf_file.write('}\n')

# GUI
def run_gui():
    try:
        # imports
        from tkinter import Button, Checkbutton, END, Entry, Frame, IntVar, Label, OptionMenu, StringVar, Tk
        from tkinter.filedialog import askdirectory, askopenfilename

        # helper function to make a popup
        def gui_popup(message, title=None):
            popup = Tk()
            if title:
                popup.wm_title(title)
            label = Label(popup, text=message)
            label.pack()
            button_close = Button(popup, text="Close", command=popup.destroy)
            button_close.pack(padx=3, pady=3)
            popup.mainloop()

        # create applet
        root = Tk()
        root.geometry("600x750")

        # set up main frame
        frame = Frame(root)
        frame.pack()

        # add header
        header = Label(frame, text="ViReflow (%s)" % VERSION, font=('Arial',24))
        header.pack()

        # handle input CSV selection
        button_csv_prefix = "Input Sample CSV:\n"
        button_csv_nofile = "<none selected>"
        def find_filename_csv():
            fn = askopenfilename(title="Select Input Samples (CSV format)",filetypes=(("CSV files","*.csv"),))
            if len(fn) == 0:
                button_csv.configure(text="%s%s" % (button_csv_prefix,button_csv_nofile))
            else:
                button_csv.configure(text="%s%s" % (button_csv_prefix,fn))
        button_csv = Button(frame, text="%s%s" % (button_csv_prefix,button_csv_nofile), command=find_filename_csv)
        button_csv.pack(padx=3, pady=3)

        # handle output folder selection
        button_out_prefix = "Output Directory:\n"
        button_out_nofolder = "<none selected>"
        def find_directory_out():
            dn = askdirectory(title="Select Output Directory")
            if len(dn) == 0:
                button_out.configure(text="%s%s" % (button_out_prefix,button_out_nofolder))
            else:
                button_out.configure(text="%s%s" % (button_out_prefix,dn))
        button_out = Button(frame, text="%s%s" % (button_out_prefix,button_out_nofolder), command=find_directory_out)
        button_out.pack(padx=3, pady=3)

        # handle S3 destination
        entry_dest_default = "Enter Destination (S3 path)"
        entry_dest = Entry(frame, width=60)
        entry_dest.insert(END, entry_dest_default)
        entry_dest.pack()

        # handle reference FASTA
        entry_ref_fas_default = "Enter Reference FASTA (S3/HTTP/HTTPS/FTP path)"
        entry_ref_fas = Entry(frame, width=60)
        entry_ref_fas.insert(END, entry_ref_fas_default)
        entry_ref_fas.pack()

        # handle reference GFF
        entry_ref_gff_default = "Enter Reference GFF (S3/HTTP/HTTPS/FTP path)"
        entry_ref_gff = Entry(frame, width=60)
        entry_ref_gff.insert(END, entry_ref_gff_default)
        entry_ref_gff.pack()

        # handle primer BED
        entry_bed_default = "Enter Primer BED (S3/HTTP/HTTPS/FTP path)"
        entry_bed = Entry(frame, width=60)
        entry_bed.insert(END, entry_bed_default)
        entry_bed.pack()

        # handle threads selection
        entry_threads_default = "Enter Number of Threads (recommended: %d)" % DEFAULT_THREADS
        entry_threads = Entry(frame, width=60)
        entry_threads.insert(END, entry_threads_default)
        entry_threads.pack()

        # handle compression level selection
        dropdown_compress_prefix = "Compression Level: "
        dropdown_compress_var = StringVar(frame)
        dropdown_compress_var.set("%s%s" % (dropdown_compress_prefix,DEFAULT_COMPRESS))
        dropdown_compress = OptionMenu(frame, dropdown_compress_var, *[("%s%d" % (dropdown_compress_prefix,i)) for i in range(1,10)])
        dropdown_compress.pack()

        # handle mapped read cap selection
        entry_cap_default = "Enter Mapped Read Cap (or None to not cap)"
        entry_cap = Entry(frame, width=60)
        entry_cap.insert(END, entry_cap_default)
        entry_cap.pack()

        # handle min alt freq selection
        entry_min_alt_default = "Enter Minimum Alternate Frequency (recommended: %s)" % DEFAULT_MIN_ALT_FREQ
        entry_min_alt = Entry(frame, width=60)
        entry_min_alt.insert(END, entry_min_alt_default)
        entry_min_alt.pack()

        # handle read mapper selection
        dropdown_mapper_prefix = "Read Mapper: "
        dropdown_mapper_var = StringVar(frame)
        dropdown_mapper_var.set("%s%s" % (dropdown_mapper_prefix,DEFAULT_READ_MAPPER))
        dropdown_mapper = OptionMenu(frame, dropdown_mapper_var, *sorted(("%s%s" % (dropdown_mapper_prefix,v)) for v in READ_MAPPERS))
        dropdown_mapper.pack()

        # handle read trimmer selection
        dropdown_trimmer_prefix = "Read Trimmer: "
        dropdown_trimmer_var = StringVar(frame)
        dropdown_trimmer_var.set("%s%s" % (dropdown_trimmer_prefix,DEFAULT_READ_TRIMMER))
        dropdown_trimmer = OptionMenu(frame, dropdown_trimmer_var, *sorted(("%s%s" % (dropdown_trimmer_prefix,v)) for v in READ_TRIMMERS_ALL))
        dropdown_trimmer.pack()

        # handle variant caller selection
        dropdown_variant_caller_prefix = "Variant Caller: "
        dropdown_variant_caller_var = StringVar(frame)
        dropdown_variant_caller_var.set("%s%s" % (dropdown_variant_caller_prefix,DEFAULT_VARIANT_CALLER))
        dropdown_variant_caller = OptionMenu(frame, dropdown_variant_caller_var, *sorted(("%s%s" % (dropdown_variant_caller_prefix,v)) for v in VARIANT_CALLERS))
        dropdown_variant_caller.pack()

        # optional header
        optional_header = Label(frame, text="\n=== Optional Analyses ===", font=('Arial',14))
        optional_header.pack()

        # handle pangolin toggle
        check_pangolin_var = IntVar(frame)
        check_pangolin = Checkbutton(frame, text=HELP_TEXT_PANGOLIN.split(" (optional)")[0].strip(), variable=check_pangolin_var, onvalue=1, offvalue=0)
        check_pangolin.pack()

        # handle virstrain selection
        dropdown_virstrain_prefix = "VirStrain DB: "
        dropdown_virstrain_var = StringVar(frame)
        dropdown_virstrain_var.set("%sNone" % dropdown_virstrain_prefix)
        dropdown_virstrain = OptionMenu(frame, dropdown_virstrain_var, *(["%sNone" % dropdown_virstrain_prefix] + sorted(("%s%s" % (dropdown_virstrain_prefix,db)) for db in VIRSTRAIN_DBS)))
        dropdown_virstrain.pack()

        # handle pi toggle
        check_pi_var = IntVar(frame)
        check_pi = Checkbutton(frame, text=HELP_TEXT_PI.split(" (optional)")[0].strip(), variable=check_pi_var, onvalue=1, offvalue=0)
        check_pi.pack()

        # handle coronaspades toggle
        check_coronaspades_var = IntVar(frame)
        check_coronaspades = Checkbutton(frame, text=HELP_TEXT_CORONASPADES.split(" (optional)")[0].strip(), variable=check_coronaspades_var, onvalue=1, offvalue=0)
        check_coronaspades.pack()

        # handle metaviralspades toggle
        check_metaviralspades_var = IntVar(frame)
        check_metaviralspades = Checkbutton(frame, text=HELP_TEXT_METAVIRALSPADES.split(" (optional)")[0].strip(), variable=check_metaviralspades_var, onvalue=1, offvalue=0)
        check_metaviralspades.pack()

        # handle rnaviralspades toggle
        check_rnaviralspades_var = IntVar(frame)
        check_rnaviralspades = Checkbutton(frame, text=HELP_TEXT_RNAVIRALSPADES.split(" (optional)")[0].strip(), variable=check_rnaviralspades_var, onvalue=1, offvalue=0)
        check_rnaviralspades.pack()

        # handle unicycler toggle
        check_unicycler_var = IntVar(frame)
        check_unicycler = Checkbutton(frame, text=HELP_TEXT_UNICYCLER.split(" (optional)")[0].strip(), variable=check_unicycler_var, onvalue=1, offvalue=0)
        check_unicycler.pack()

        # handle megahit toggle
        check_megahit_var = IntVar(frame)
        check_megahit = Checkbutton(frame, text=HELP_TEXT_MEGAHIT.split(" (optional)")[0].strip(), variable=check_megahit_var, onvalue=1, offvalue=0)
        check_megahit.pack()

        # handle minia toggle
        check_minia_var = IntVar(frame)
        check_minia = Checkbutton(frame, text=HELP_TEXT_MINIA.split(" (optional)")[0].strip(), variable=check_minia_var, onvalue=1, offvalue=0)
        check_minia.pack()

        # add generate button
        def finish_applet():
            valid = True
            # check CSV
            try:
                if button_csv['text'] == "%s%s" % (button_csv_prefix,button_csv_nofile):
                    gui_popup("ERROR: Input Sample CSV file not selected", title="ERROR"); valid = False
            except:
                pass
            # check output directory
            try:
                if button_out['text'] == "%s%s" % (button_out_prefix,button_out_nofolder):
                    gui_popup("ERROR: Output Directory not selected", title="ERROR"); valid = False
            except:
                pass
            # check destination
            try:
                if not entry_dest.get().strip().lower().startswith('s3://'):
                    gui_popup("ERROR: Invalid destination S3 path:\n%s" % entry_dest.get()); valid = False
            except:
                pass
            # check reference FASTA
            try:
                if entry_ref_fas.get().strip().lower().split('://')[0] not in {'s3', 'http', 'https', 'ftp'}:
                    gui_popup("ERROR: Invalid reference FASTA path:\n%s" % entry_ref_fas.get()); valid = False
            except:
                pass
            # check reference GFF
            try:
                if entry_ref_gff.get().strip().lower().split('://')[0] not in {'s3', 'http', 'https', 'ftp'}:
                    gui_popup("ERROR: Invalid reference GFF path:\n%s" % entry_ref_gff.get()); valid = False
            except:
                pass
            # check primer BED
            try:
                if entry_bed.get().strip().lower().split('://')[0] not in {'s3', 'http', 'https', 'ftp'}:
                    gui_popup("ERROR: Invalid primer BED path:\n%s" % entry_bed.get()); valid = False
            except:
                pass
            # check threads
            try:
                if not check_int(entry_threads.get()) or int(entry_threads.get()) < 1:
                    gui_popup("ERROR: Invalid number of threads:\n%s" % entry_threads.get()); valid = False
            except:
                pass
            # check mapped read cap
            try:
                if entry_cap.get() != 'None' and (not check_int(entry_cap.get()) or int(entry_cap.get()) < 1):
                    gui_popup("ERROR: Invalid mapped read cap:\n%s" % entry_cap.get()); valid = False
            except:
                pass
            # check min alt freq
            try:
                if not check_float(entry_min_alt.get()) or float(entry_min_alt.get()) < 0 or float(entry_min_alt.get()) > 1:
                    gui_popup("ERROR: Invalid minimum alternate allele frequency:\n%s" % entry_min_alt.get()); valid = False
            except:
                pass
            # close applet to run ViReflow
            if valid:
                sys.argv.append('-d'); sys.argv.append(entry_dest.get().strip())
                sys.argv.append('-rf'); sys.argv.append(entry_ref_fas.get().strip())
                sys.argv.append('-rg'); sys.argv.append(entry_ref_gff.get().strip())
                sys.argv.append('-p'); sys.argv.append(entry_bed.get().strip())
                sys.argv.append('-t'); sys.argv.append(entry_threads.get().strip())
                sys.argv.append('-cl'); sys.argv.append(dropdown_compress_var.get().split(':')[-1].strip())
                #try:
                #    read_cap = int(entry_cap.get())
                #    sys.argv.append('--mapped_read_cap'); sys.argv.append(str(read_cap))
                #except:
                #    pass
                sys.argv.append('--min_alt_freq'); sys.argv.append(entry_min_alt.get().strip())
                sys.argv.append('--read_mapper'); sys.argv.append(dropdown_mapper_var.get().split(':')[-1].strip())
                sys.argv.append('--read_trimmer'); sys.argv.append(dropdown_trimmer_var.get().split(':')[-1].strip())
                sys.argv.append('--variant_caller'); sys.argv.append(dropdown_trimmer_var.get().split(':')[-1].strip())
                if dropdown_virstrain_var.get().split(':')[-1].strip() != 'None':
                    sys.argv.append('--optional_virstrain'); sys.argv.append(dropdown_virstrain.get().split(':')[-1].strip())
                if check_pangolin_var.get() == 1:
                    sys.argv.append('--optional_pangolin')
                if check_pi_var.get() == 1:
                    sys.argv.append('--optional_pi_metric')
                if check_coronaspades_var.get() == 1:
                    sys.argv.append('--optional_spades_coronaspades')
                if check_metaviralspades_var.get() == 1:
                    sys.argv.append('--optional_spades_metaviralspades')
                if check_rnaviralspades_var.get() == 1:
                    sys.argv.append('--optional_spades_rnaviralspades')
                if check_unicycler_var.get() == 1:
                    sys.argv.append('--optional_unicycler')
                if check_megahit_var.get() == 1:
                    sys.argv.append('--optional_megahit')
                if check_minia_var.get() == 1:
                    sys.argv.append('--optional_minia')
                sys.argv.append(button_csv['text'].split(':')[-1].strip())
                chdir(button_out['text'].split(':')[-1].strip())
                try:
                    root.destroy()
                except:
                    pass
        button_generate = Button(frame, text="Generate Reflow run files", command=finish_applet)
        button_generate.pack(padx=3, pady=3)

        # add title and execute GUI
        root.title("ViReflow (%s)" % VERSION)
        root.mainloop()
    except:
        print("ERROR: Unable to import Tkinter", file=stderr); exit(1)
    if len(sys.argv) == 1:
        exit()

# when script is executed
if __name__ == "__main__":
    if len(sys.argv) == 1:
        run_gui()
    main()
