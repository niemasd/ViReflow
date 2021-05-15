#! /usr/bin/env python3
'''
Create a Reflow run file that batch-executes many other Reflow run files.
'''
import argparse
from os.path import abspath, isfile
from sys import stderr, stdout

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Batch Reflow Run File (RF)")
    parser.add_argument('rf_files', metavar='RF', type=str, nargs='+', help="Input Reflow Run Files (RF)")
    args = parser.parse_args()
    for rf in args.rf_files:
        if not isfile(rf):
            stderr.write("ERROR: File not found: %s\n" % args.output); exit(1)
    if args.output == 'stdout':
        args.output = stdout
    else:
        if isfile(args.output):
            stderr.write("ERROR: Output file exists: %s\n" % args.output); exit(1)
        args.output = open(args.output, 'w')
    return args

# main content
if __name__ == "__main__":
    args = parse_args()
    args.output.write("val Main = [\n")
    for rf in args.rf_files:
        args.output.write('    make("%s").Main,\n' % abspath(rf))
    args.output.write("]\n")
    args.output.close()
