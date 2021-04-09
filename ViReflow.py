#! /usr/bin/env python3
'''
ViReflow: An elastically-scaling automated AWS pipeline for viral consensus sequence generation
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
    'minimap2': {
        'docker_image':  'niemasd/minimap2:latest', # Docker image for Minimap2
        'num_cpu_index': 1,                         # AWS instance num CPUs for indexing reference genome (Minimap2 indexing doesn't benefit from more than 1 CPU for viral genomes)
        'mem_index':     '50*MiB',                  # AWS instance memory for indexing reference genome (Minimap2 Peak RSS to index the SARS-CoV-2 reference was 0.003 GB)
    }
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
    parser.add_argument('-r', '--reference', required=True, type=str, help="Reference")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Reflow File")
    parser.add_argument('-u', '--update', action="store_true", help="Update ViReflow (current version: %s)" % VERSION)
    args = parser.parse_args()

    # user args are valid, so return
    return args

# main content
if __name__ == "__main__":
    # parse user args and prepare run
    args = parse_args()

    # handle output file
    if args.output == 'stdout':
        from sys import stdout as rf_file
    else:
        rf_file = open(args.output, 'w')
    rf_file.write('// Created using ViReflow %s\n\n' % VERSION)

    # handle reference
    if args.reference.lower().startswith('s3://'): # Amazon S3 path
        if not args.reference.lower().endswith('.mmi'):
            stderr.write("Invalid Minimap2 annotation file extension (should be .mmi): %s\n" % args.reference)
            if args.output != 'stdout':
                rf_file.close(); remove(args.output)
            exit(1)
        rf_file.write('// using "%s" as reference\n' % args.reference)
        rf_file.write('val ref_mmi = file("%s")\n' % args.reference)
        rf_file.flush()
    else: # assume GenBank accession number
        ref_fasta_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta' % args.reference
        try:
            urlopen(ref_fasta_url)
        except:
            stderr.write("Invalid reference genome GenBank accession: %s\n" % args.reference)
            if args.output != 'stdout':
                rf_file.close(); remove(args.output)
            exit(1)
        rf_file.write('// create Minimap2 reference from GenBank accession "%s"\n' % args.reference)
        rf_file.write('val ref_mmi = exec(image := "%s", mem := %s, cpu := %d) (out file) {"\n' % (TOOL['minimap2']['docker_image'], TOOL['minimap2']['mem_index'], TOOL['minimap2']['num_cpu_index'])) # TODO
        rf_file.write('    wget -O "%s.fas" "%s"\n' % (args.reference, ref_fasta_url))
        rf_file.write('    minimap2 -t %d -d {{out}} "%s.fas"\n' % (TOOL['minimap2']['num_cpu_index'], args.reference))
        rf_file.write('"}\n')
    rf_file.write('\n')
    exit(1) # TODO
