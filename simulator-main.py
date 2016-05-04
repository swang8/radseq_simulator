from Util import *
import sys
import argparse

## parse command line parameters
parser = argparse.ArgumentParser(description='A read simulator for RADseq, introducing variations and sequencing errors.')
parser.add_argument('-ref', help='reference file in fasta format', type=str)
parser.add_argument('-depth', help='coverage depth to be simulated', type=int)
parser.add_argument('-len', help='the length of reads to be generated', type=int)
parser.add_argument('-type', help='single-end (S) or pair-end (P), default is single-end', type=str, default='S')
parser.add_argument('-e1', help='The 1st enzyme', type=str)
parser.add_argument('-e2', help='The 2nd enzyme, optional', type=str)

args = parser.parse_args()

if args.type not in ('s', 'S', 'p', 'P'):
    print 'Use S for single end; P for pair end'
    exit(1);

(ref_fasta, depth, read_len, read_type, enz_first, enz_second) = (args.ref, args.depth, args.len, args.type, args.e1, args.e2)

