import Util
import sys
import argparse
import datetime
from Bio.Seq import Seq
from Bio import SeqIO

## parse command line parameters
parser = argparse.ArgumentParser(description='A read simulator for RADseq, introducing variations and sequencing errors.')
parser.add_argument('-ref', help='reference file in fasta format', type=str)
parser.add_argument('-depth', help='coverage depth to be simulated', type=int)
parser.add_argument('-len', help='the length of reads to be generated', type=int)
parser.add_argument('-type', help='single-end (S) or pair-end (P), default is single-end', type=str, default='S')
parser.add_argument('-e1', help='The 1st enzyme', type=str)
parser.add_argument('-e2', help='The 2nd enzyme, optional', type=str)
parser.add_argument('-output', help='prefix of the output file in fastq format', type=str)

args = parser.parse_args()

if args.type not in ('s', 'S', 'p', 'P'):
    print 'Use S for single end; P for pair end'
    print 'Not recognizing type: ' + args.type
    print 'Will generate single end reads!'

(ref_fasta, depth, read_len, read_type, enz_first, enz_second, output_fasta) = \
    (args.ref, args.depth, args.len, args.type, args.e1, args.e2, args.output)

read_len += 10
## check the enzyme(s)

if not (enz_first and Util.enzyme.enzymeDict().checkEnz(enz_first)):
    print 'The first Enzyme not provided or not in the enzyme dictionary'
    print Util.enzyme.enzymeDict().ids()
    parser.print_help()
    exit(1)

if enz_second and not Util.enzyme.enzymeDict().checkEnz(enz_second):
    print 'Enzyme ' + enz_second + ' not in the enzyme dictionary'
    print Util.enzyme.enzymeDict().ids()
    parser.print_help()
    exit(1)

## get the fragments
digest = Util.digest.Digst([enz_first], ref_fasta)
if enz_second: digest.addEnzyme(enz_second)
## start digestion
print_current_time('Start digestion!')
digest.runDigest()
print_current_time('Done digestion!')

fragments = digest.getFragments()

## order the positions
fragments_ordered = {}
for ctg in fragments:
    if ctg not in fragments_ordered:
        fragments_ordered[ctg] = []
    for enz in fragments[ctg]:
        for i  in fragments[ctg][enz]:
            fragments_ordered[ctg].append([i, enz])
    fragments_ordered[ctg].sort()

## generate fasta file

out_fh1 = open (output_fasta + "_R1.fastq", 'w')
out_fh2 = ''
if read_type in ['p', 'P']: out_fh2 = open (output_fasta + "_R2.fastq", 'w')

enz_leftover =get_re_leftover([enz_first, enz_second])

for seq_record in SeqIO.parse(ref_fasta, "fasta"):
    #print(seq_record.id)
    #print(repr(seq_record.seq))
    #print(len(seq_record))
    if seq_record.id in fragments_ordered:
        (mut_seq, mut_pos) = Util.simulate.mutate_genome(seq_record.seq)
        for frag_index in range(len(fragments_ordered[seq_record.id])):
            frag = fragments_ordered[seq_record.id][frag_index]
            (pos, enz) = frag

            if type in ['s', 'S']: ## single end
                if enz == enz_first:
                    leftover_seq = enz_leftover[enz]
                    # forward
                    forward_read = leftover_seq + mut_seq[pos: read_len + pos - len(leftover_seq)]
                    out_fh1.write(">" + seq_record.id + ':' + pos + "\n")
                    out_fh1.write(forward_read + '\n')
                    # backward
                    back_pos = pos - len(Util.enzyme.enzymeDict().getProperty(enz)['site'])
                    back_start = back_pos - read_len
                    if back_start < 0:
                        back_start = 0
                    backward_read = leftover_seq + Seq(mut_seq[back_start:back_pos]).reverse_complment()
                    out_fh1.write('>' + seq_record.id + ':' + back_pos + '\n' + backward_read + '\n')
            else: ## pair end
                



out_fh1.close()
out_fh2.close()

## helpers
def get_re_leftover (enzs):
    leftover = {}
    for enz in enzs:
        if not enz:
            continue
        site = Util.enzyme.enzymeDict().getProperty(enz)['site']
        ovhg = Util.enzyme.enzymeDict().getProperty(enz)['ovhgseq']
        leftover[enz] = site[site.index(ovhg):]
    return leftover

def print_current_time(s):
    print datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), "\t", s