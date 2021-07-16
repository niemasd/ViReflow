#! /usr/bin/env python3
'''
ViReflow: An elastically-scaling parallelized AWS pipeline for viral consensus sequence generation
'''

# imports
from json import load as jload
from os import remove
from sys import argv, stderr
from urllib.request import urlopen
import argparse

# useful constants
VERSION = '1.0.10'
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
READ_MAPPERS = {'bowtie2', 'bwa', 'minimap2'}
VARIANT_CALLERS = {'freebayes', 'ivar', 'lofreq'}
INSTANCE_INFO = {
    'docker_image': 'niemasd/vireflow:latest', # TODO CHANGE TO VERSION ONCE STABLE
    'mem':          '1*GiB',                   # TODO FIND APPROPRIATE VALUE
    'disk':         '3*GiB',                   # TODO FIND APPROPRIATE VALUE
}
DATE_COMMAND_BASH = 'echo "[$(date +"%Y-%m-%d %T")]"'

# convert a ViReflow version string to a tuple of integers
def parse_version(s):
    return tuple(int(v) for v in s.split('.'))

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
    if '-u' in argv or '--update' in argv:
        update_vireflow()

    # use argparse to parse user arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-id', '--run_id', required=True, type=str, help="Unique Run Identifier (for output file naming)")
    parser.add_argument('-d', '--destination', required=True, type=str, help="Destination for Results (s3 folder)")
    parser.add_argument('-rf', '--reference_fasta', required=True, type=str, help="Reference Genome Sequence (s3/http/https/ftp to FASTA)")
    parser.add_argument('-rg', '--reference_gff', required=True, type=str, help="Reference Genome Annotation (s3/http/https/ftp to GFF3)")
    parser.add_argument('-p', '--primer_bed', required=True, type=str, help="Primer (s3/http/https/ftp to BED)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Reflow File (rf)")
    parser.add_argument('-t', '--threads', required=False, type=int, default=1, help="Number of Threads")
    parser.add_argument('-cl', '--compression_level', required=False, type=int, default=1, help="Compression Level (1 = fastest, 9 = best)")
    parser.add_argument('--mapped_read_cap', required=False, type=int, default=None, help="Successfully-Mapped Read Cap")
    parser.add_argument('--min_alt_freq', required=False, type=float, default=0.5, help="Minimum Alt Allele Frequency for consensus sequence")
    parser.add_argument('--read_mapper', required=False, type=str, default='minimap2', help="Read Mapper (options: %s)" % ', '.join(sorted(READ_MAPPERS)))
    parser.add_argument('--read_trimmer', required=False, type=str, default='ivar', help="Read Trimmer (options: %s)" % ', '.join(sorted(READ_TRIMMERS_ALL)))
    parser.add_argument('--variant_caller', required=False, type=str, default='lofreq', help="Variant Caller (options: %s)" % ', '.join(sorted(VARIANT_CALLERS)))
    parser.add_argument('--optional_pangolin', action="store_true", help="Run Pangolin (optional)")
    parser.add_argument('-u', '--update', action="store_true", help="Update ViReflow (current version: %s)" % VERSION)
    parser.add_argument('fastq_files', metavar='FQ', type=str, nargs='+', help="Input FASTQ Files (s3/http/https/ftp; single biological sample)")
    args = parser.parse_args()
    args.variant_caller = args.variant_caller.lower()

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

    # check 
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

    # user args are valid, so return
    return args

# main content
if __name__ == "__main__":
    # parse user args
    args = parse_args()
    INSTANCE_INFO['cpu'] = args.threads

    # check input files (FASTQs, ref FASTA, ref GFF, and primer BED)
    input_files = [
        ('ref_fas',    '%s.reference.fas' % args.run_id, args.reference_fasta),
        ('ref_gff',    '%s.reference.gff' % args.run_id, args.reference_gff),
        ('primer_bed', '%s.primers.bed'   % args.run_id, args.primer_bed),
    ]
    input_files += [('fq%d' % (i+1), '%s.fq%d.fastq' % (args.run_id, i+1), f) for i,f in enumerate(args.fastq_files)]

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
    rf_file.write('// Created using ViReflow %s\n' % VERSION)
    rf_file.write('// ViReflow command: %s\n' % ' '.join(argv))
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
        for rf_var, local_fn, p in input_files:
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
        rf_file.write('        bwa mem -t %d "%s" %s' % (args.threads, local_fns['ref_fas'], ' '.join('"%s"' % fn for fn in local_fq_fns_to_map)))
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

    # remove redundant files before compressing output
    rf_file.write('        # Remove redundant output files before compressing\n')
    rf_file.write('        %s "Removing redundant output files before compressing" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
    rf_file.write('        rm */*trimmed.bam\n') # remove unsorted BAMs
    rf_file.write('\n')

    # archive + compress output files
    rf_file.write('        # Compress output files\n')
    rf_file.write('        %s "Compressing output files" >> %s\n' % (DATE_COMMAND_BASH, local_fns['vireflow_log']))
    rf_file.write('        tar cvf - "%s" | pigz -%d -p %d > "{{out}}"\n' % (outdir, args.compression_level, args.threads))

    # end main exec, copy output file, and end run file
    rf_file.write('    "}\n')
    rf_file.write('    cp_out_file := files.Copy(out_file, "%s")\n' % final_output_s3)
    rf_file.write('    (cp_out_file)\n')
    rf_file.write('}\n')
