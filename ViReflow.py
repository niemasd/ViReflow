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
VERSION = '0.0.1'
RELEASES_URL = 'https://api.github.com/repos/niemasd/ViReflow/tags'
TOOL = {
    'base': {
        'docker_image':  'alpine:latest',                    # Base Docker image (Alpine)
        'cpu_wget':      1,                                  # Num CPUs for wget
        'mem_wget':      '50*MiB',                           # Memory for wget
        'disk':          '5*GiB',                            # Overall, shouldn't need more than 5 GB of disk
    },

    'fastqc': {
        'docker_image':  'niemasd/fastqc:latest',            # Docker image for fastqc
        'cpu': 1,                                            # Num CPUs (just going to run single-threaded for memory reasons)
        'mem': '2*GiB',                                      # Memory per FASTQ (Java allocates a ton of memory)
    },

    'ivar': {
        'docker_image':  'niemasd/ivar:latest',              # Docker image for iVar
        'cpu_trim':      1,                                  # Num CPUs for trimming (iVar is still single-threaded)
        'mem_trim':      '50*MiB',                           # Memory for trimming (takes 7 MB on demo)
        'cpu_variants':  1,                                  # Num CPUs for variant-calling (iVar is still single-threaded)
        'mem_variants':  '20*MiB',                           # Memory for variant-calling (takes 1 MB on demo)
        'cpu_consensus': 1,                                  # Num CPUs for consensus-calling (iVar is still single-threaded)
        'mem_consensus': '20*MiB',                           # Memory for consensus-calling (takes 1 MB on demo)
    },

    'minimap2': {
        'docker_image':  'niemasd/minimap2:latest',          # Docker image for Minimap2
        'cpu_index':     1,                                  # Num CPUs for indexing reference genome (Minimap2 indexing doesn't benefit from more than 1 CPU for viral genomes)
        'mem_index':     '50*MiB',                           # Memory for indexing reference genome (Minimap2 Peak RSS to index the SARS-CoV-2 reference was 0.003 GB)
        'cpu_map':       32,                                 # Num CPUs for mapping reads (can be increased/decreased by user as desired)
        'mem_map':       '100*MiB',                          # Memory for mapping reads (takes 30 MB on demo)
    },

    'minimap2_samtools': {
        'docker_image':  'niemasd/minimap2_samtools:latest', # Docker image for Minimap2 + samtools
    },

    'samtools': {
        'docker_image':  'niemasd/samtools:latest',          # Docker image for samtools
        'cpu_sort':      32,                                 # Num CPUs for sorting BAM
        'mem_sort':      '100*MiB',                          # Memory for sorting BAM (takes 20 MB on demo)
        'cpu_pileup':    1,                                  # Num CPUs for generating pileup
        'mem_pileup':    '50*MiB',                           # Memory for generating pileup (takes 5 MB on demo)
        'cpu_depth':     1,                                  # Num CPUs for computing depth
        'mem_depth':     '50*MiB',                           # Memory for computing depth (takes 4 MB on demo)
    },
}

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
    parser.add_argument('-d', '--destination', required=True, type=str, help="Destination for Results (s3 folder)")
    parser.add_argument('-rf', '--reference_fasta', required=True, type=str, help="Reference Genome Sequence (s3/http/https/ftp to FASTA)")
    parser.add_argument('-rm', '--reference_mmi', required=False, type=str, default=None, help="Reference Genome Minimap2 Index (s3 to MMI)")
    parser.add_argument('-rg', '--reference_gff', required=True, type=str, help="Reference Genome Annotation (s3/http/https/ftp to GFF3)")
    parser.add_argument('-p', '--primer_bed', required=True, type=str, help="Primer (s3/http/https/ftp to BED)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Reflow File (rf)")
    parser.add_argument('-mt', '--max_threads', required=False, type=int, default=TOOL['minimap2']['cpu_map'], help="Max Threads")
    parser.add_argument('--include_fastqc', action="store_true", help="Include FastQC")
    parser.add_argument('--include_depth', action="store_true", help="Include Depth Calling")
    parser.add_argument('-u', '--update', action="store_true", help="Update ViReflow (current version: %s)" % VERSION)
    parser.add_argument('fastq_files', metavar='FQ', type=str, nargs='+', help="Input FASTQ Files (s3 paths; single biological sample)")
    args = parser.parse_args()

    # user args are valid, so return
    return args

# main content
if __name__ == "__main__":
    # parse user args and prepare run
    args = parse_args()
    out_list = ['cp_sorted_untrimmed_bam', 'cp_sorted_trimmed_bam', 'cp_pileup', 'cp_variants', 'cp_consensus']

    # check max threads
    if args.max_threads < 1:
        stderr.write("Invalid maximum number of threads: %d\n" % args.max_threads); exit(1)
    for k1 in TOOL:
        for k2 in TOOL[k1]:
            if k2.startswith('cpu_') and TOOL[k1][k2] != 1: # some tools only run single-threaded
                TOOL[k1][k2] = args.max_threads

    # check destination folder
    if not args.destination.lower().startswith('s3://'):
        stderr.write("Invalid output s3 directory: %s\n" % args.destination); exit(1)
    args.destination = args.destination.rstrip('/')

    # handle output file
    if args.output == 'stdout':
        from sys import stdout as rf_file
    else:
        rf_file = open(args.output, 'w')
    rf_file.write('// Created using ViReflow %s\n' % VERSION)
    rf_file.write('@requires(disk := %s)\n' % TOOL['base']['disk'])
    rf_file.write('val Main = {\n')
    rf_file.write('    files := make("$/files")\n')

    # handle input FASTQs
    fqs = list() # (Reflow variable, s3 path) tuples
    for i,fq in enumerate(args.fastq_files):
        if not fq.lower().startswith('s3://'):
            stderr.write("Invalid s3 path to FASTQ file: %s" % fq)
            if args.output != 'stdout':
                rf_file.close(); remove(args.output)
            exit(1)
        fqs.append(('fq%d' % (i+1), fq))
    rf_file.write('    // Use the following FASTQ s3 path(s)\n')
    for tup in fqs:
        rf_file.write('    %s := file("%s")\n' % tup)
    rf_file.write('\n')

    # handle reference sequence
    ref_fas_lower = args.reference_fasta.lower()
    rf_file.write('    // Use reference FASTA file: %s\n' % args.reference_fasta)
    rf_file.write('    ref_fas := ')
    if ref_fas_lower.startswith('s3://'): # Amazon S3 path
        rf_file.write('file("%s")' % args.reference_fasta)
    else:
        if ref_fas_lower.startswith('http://') or ref_fas_lower.startswith('https://') or ref_fas_lower.startswith('ftp://'): # URL of FASTA
            ref_fasta_url = args.reference_fasta
        else: # assume GenBank accession number
            ref_fasta_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta' % args.reference_fasta
        try:
            urlopen(ref_fasta_url)
        except:
            stderr.write("Invalid reference genome FASTA: %s\n" % args.reference_fasta)
            if args.output != 'stdout':
                rf_file.close(); remove(args.output)
            exit(1)
        rf_file.write('exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['base']['docker_image'], TOOL['base']['mem_wget'], TOOL['base']['cpu_wget']))
        rf_file.write('        wget -O "{{out}}" "%s" 1>&2\n' % ref_fasta_url)
        rf_file.write('    "}\n')
        rf_file.write('    cp_ref_fas := files.Copy(ref_fas, "%s/reference.fas")' % args.destination)
        out_list.append('cp_ref_fas')
    rf_file.write('\n\n')

    # handle reference index
    if args.reference_mmi is None:
        rf_file.write('    // Create Minimap2 reference index\n')
        rf_file.write('    ref_mmi := exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['minimap2']['docker_image'], TOOL['minimap2']['mem_index'], TOOL['minimap2']['cpu_index']))
        rf_file.write('        minimap2 -t %d -d "{{out}}" "{{ref_fas}}" 1>&2\n' % TOOL['minimap2']['cpu_index'])
        rf_file.write('    "}\n')
        rf_file.write('    cp_ref_mmi := files.Copy(ref_mmi, "%s/reference.mmi")' % args.destination)
        out_list.append('cp_ref_mmi')
    else:
        ref_mmi_lower = args.reference_mmi.lower()
        rf_file.write('    // Use existing Minimap2 reference index: %s\n' % args.reference_mmi)
        rf_file.write('    ref_mmi := ')
        if ref_mmi_lower.startswith('s3://'): # Amazon S3 path
            rf_file.write('file("%s")' % args.reference_mmi)
        else:
            try:
                urlopen(args.reference_mmi)
            except:
                stderr.write("Invalid reference genome MMI: %s\n" % args.reference_mmi)
                if args.output != 'stdout':
                    rf_file.close(); remove(args.output)
                exit(1)
            rf_file.write('exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['base']['docker_image'], TOOL['base']['mem_wget'], TOOL['base']['cpu_wget']))
            rf_file.write('        wget -O "{{out}}" "%s" 1>&2\n' % args.reference_mmi)
            rf_file.write('    "}\n')
            rf_file.write('    cp_ref_mmi := files.Copy(ref_mmi, "%s/reference.mmi")' % args.destination)
            out_list.append('cp_ref_mmi')
    rf_file.write('\n\n')

    # handle reference annotation
    ref_gff_lower = args.reference_gff.lower()
    rf_file.write('    // Use reference GFF3 file: %s\n' % args.reference_gff)
    rf_file.write('    ref_gff := ')
    if ref_gff_lower.startswith('s3://'): # Amazon S3 path
        rf_file.write('file("%s")' % args.reference_gff)
    else:
        try:
            urlopen(args.reference_gff)
        except:
            stderr.write("Invalid reference genome FASTA: %s\n" % args.reference_gff)
            if args.output != 'stdout':
                rf_file.close(); remove(args.output)
            exit(1)
        rf_file.write('exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['base']['docker_image'], TOOL['base']['mem_wget'], TOOL['base']['cpu_wget']))
        rf_file.write('        wget -O "{{out}}" "%s" 1>&2\n' % args.reference_gff)
        rf_file.write('    "}\n')
        rf_file.write('    cp_ref_gff := files.Copy(ref_gff, "%s/reference.gff")' % args.destination)
        out_list.append('cp_ref_gff')
    rf_file.write('\n\n')

    # handle primer BED
    primer_bed_lower = args.primer_bed.lower()
    rf_file.write('    // Use primer BED file: %s\n' % args.primer_bed)
    rf_file.write('    primer_bed := ')
    if primer_bed_lower.startswith('s3://'): # Amazon S3 path
        rf_file.write('file("%s")' % args.primer_bed)
    else:
        try:
            urlopen(args.primer_bed)
        except:
            stderr.write("Invalid primer BED: %s" % args.primer_bed)
            if args.output != 'stdout':
                rf_file.close(); remove(args.output)
            exit(1)
        rf_file.write('exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['base']['docker_image'], TOOL['base']['mem_wget'], TOOL['base']['cpu_wget']))
        rf_file.write('        wget -O "{{out}}" "%s" 1>&2\n' % args.primer_bed)
        rf_file.write('    "}\n')
        rf_file.write('    cp_primer_bed := files.Copy(primer_bed, "%s/primers.bed")' % args.destination)
        out_list.append('cp_primer_bed')
    rf_file.write('\n\n')

    # map reads using Minimap2 and sort using samtools
    rf_file.write('    // Map reads using Minimap2 and sort using samtools\n')
    rf_file.write('    sorted_untrimmed_bam := exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['minimap2_samtools']['docker_image'], TOOL['minimap2']['mem_map'], TOOL['minimap2']['cpu_map']))
    rf_file.write('        minimap2 -t %d -a -x sr "{{ref_mmi}}" %s | samtools sort --threads %d -o "{{out}}" 1>&2\n' % (TOOL['minimap2']['cpu_map'], ' '.join('"{{%s}}"' % var for var,s3 in fqs), TOOL['samtools']['cpu_sort']))
    rf_file.write('    "}\n')
    rf_file.write('    cp_sorted_untrimmed_bam := files.Copy(sorted_untrimmed_bam, "%s/sorted.untrimmed.bam")\n' % args.destination)
    rf_file.write('\n')

    # trim reads using iVar
    rf_file.write('    // Trim reads using iVar\n')
    rf_file.write('    trimmed_bam := exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['ivar']['docker_image'], TOOL['ivar']['mem_trim'], TOOL['ivar']['cpu_trim']))
    rf_file.write('        ivar trim -x 5 -e -i "{{sorted_untrimmed_bam}}" -b "{{primer_bed}}" -p trimmed 1>&2 && mv trimmed.bam "{{out}}"\n')
    rf_file.write('    "}\n\n')

    # sort trimmed BAM
    rf_file.write('    // Sort trimmed BAM\n')
    rf_file.write('    sorted_trimmed_bam := exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['samtools']['docker_image'], TOOL['samtools']['mem_sort'], TOOL['samtools']['cpu_sort']))
    rf_file.write('        samtools sort --threads %d -o "{{out}}" "{{trimmed_bam}}" 1>&2\n' % TOOL['samtools']['cpu_sort'])
    rf_file.write('    "}\n')
    rf_file.write('    cp_sorted_trimmed_bam := files.Copy(sorted_trimmed_bam, "%s/sorted.trimmed.bam")\n' % args.destination)
    rf_file.write('\n')

    # generate pile-up from sorted trimmed BAM
    rf_file.write('    // Generate pile-up from sorted trimmed BAM\n')
    rf_file.write('    pileup := exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['samtools']['docker_image'], TOOL['samtools']['mem_pileup'], TOOL['samtools']['cpu_pileup']))
    rf_file.write('        samtools mpileup -A -aa -d 0 -Q 0 --reference "{{ref_fas}}" "{{sorted_trimmed_bam}}" > "{{out}}"\n')
    rf_file.write('    "}\n')
    rf_file.write('    cp_pileup := files.Copy(pileup, "%s/pileup.txt")\n' % args.destination)
    rf_file.write('\n')

    # call variants from pile-up
    rf_file.write('    // Call variants from pile-up\n')
    rf_file.write('    variants := exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['ivar']['docker_image'], TOOL['ivar']['mem_variants'], TOOL['ivar']['cpu_variants']))
    rf_file.write('        cp "{{ref_fas}}" ref.fas && cp "{{ref_gff}}" ref.gff\n')
    rf_file.write('        cat "{{pileup}}" | ivar variants -r ref.fas -g ref.gff -p tmp.tsv -m 10 1>&2\n')
    rf_file.write('        mv tmp.tsv "{{out}}"\n')
    rf_file.write('    "}\n')
    rf_file.write('    cp_variants := files.Copy(variants, "%s/variants.tsv")\n' % args.destination)
    rf_file.write('\n')

    # call consensus from pile-up
    rf_file.write('    // Call consensus from pile-up\n')
    rf_file.write('    consensus := exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['ivar']['docker_image'], TOOL['ivar']['mem_consensus'], TOOL['ivar']['cpu_consensus']))
    rf_file.write('        cat "{{pileup}}" | ivar consensus -p consensus -m 10 -n N -t 0.5 1>&2 && mv consensus.fa "{{out}}" 1>&2\n')
    rf_file.write('    "}\n')
    rf_file.write('    cp_consensus := files.Copy(consensus, "%s/consensus.fas")\n' % args.destination)
    rf_file.write('\n')

    # run FastQC (optional)
    if args.include_fastqc:
        out_list.append('cp_fastqc')
        rf_file.write('    // Run FastQC\n')
        rf_file.write('    fastqc := exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['fastqc']['docker_image'], TOOL['fastqc']['mem'], TOOL['fastqc']['cpu']))
        rf_file.write('        mkdir fastqc\n')
        rf_file.write('        fastqc -o fastqc %s\n' % ' '.join('"{{%s}}"' % var for var,s3 in fqs))
        rf_file.write('        zip -9 "{{out}}" fastqc/*\n')
        rf_file.write('    "}\n')
        rf_file.write('    cp_fastqc := files.Copy(fastqc, "%s/fastqc.zip")\n' % args.destination)
        rf_file.write('\n')

    # call depth (optional)
    if args.include_depth:
        out_list.append('cp_depth')
        rf_file.write('    // Call depth from trimmed BAM\n')
        rf_file.write('    depth := exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['samtools']['docker_image'], TOOL['samtools']['mem_depth'], TOOL['samtools']['cpu_depth']))
        rf_file.write('        samtools depth -d 0 -Q 0 -q 0 -aa "{{sorted_trimmed_bam}}" > "{{out}}"\n')
        rf_file.write('    "}\n')
        rf_file.write('    cp_depth := files.Copy(depth, "%s/depth.txt")\n' % args.destination)
        rf_file.write('\n')

    # finish Main
    rf_file.write('    // Finish workflow\n')
    rf_file.write('    (%s)\n' % ', '.join(out_list))
    rf_file.write('}\n')
