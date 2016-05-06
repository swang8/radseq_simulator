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

(ref_fasta, depth, read_len, read_type, enz_first, enz_second, output_prefix) = \
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

out_fh1 = open (output_prefix + "_R1.fastq", 'w')
out_fh2 = ''
if read_type in ['p', 'P']: out_fh2 = open (output_prefix + "_R2.fastq", 'w')

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

            # generate reads
            if type in ['s', 'S']:
                ## single end
                if enz == enz_first:
                    leftover_seq = enz_leftover[enz]
                    # forward
                    forward_read = leftover_seq + mut_seq[pos: read_len + pos - len(leftover_seq)]
                    read_qual = Util.simulate.generate_phred_scores(len(forward_read))
                    out_fh1.write("@" + seq_record.id + ':' + pos + " 1\n")
                    out_fh1.write(forward_read + '\n+\n' +read_qual + '\n')
                    # backward
                    back_pos = pos - len(Util.enzyme.enzymeDict().getProperty(enz)['site'])
                    back_start = back_pos - read_len
                    if back_start < 0:
                        back_start = 0
                    backward_read = leftover_seq + Seq(mut_seq[back_start:back_pos]).reverse_complment()
                    read_qual = Util.simulate.generate_phred_scores(len(backward_read))
                    out_fh1.write('@' + seq_record.id + ':' + back_pos + ' 1\n' + backward_read + '\n+\n', read_qual + '\n')

            if type in ['p', 'P']:
                ## pair end
                if enz == enz_first:
                    # plus strand
                    frag_next = frag_index + 1
                    (next_pos, next_enz) = ()
                    if frag_next < len(fragments_ordered[seq_record.id]):
                        (next_pos, next_enz) = fragments_ordered[seq_record.id][frag_next]
                    else:
                        next_pos = len(mut_seq) - 1
                    if next_enz == enz_second:
                        leftover_seq = enz_leftover[enz]
                        # forward
                        forward_read = leftover_seq + mut_seq[pos: read_len + pos - len(leftover_seq)]
                        read_qual = Util.simulate.generate_phred_scores(len(forward_read))
                        out_fh1.write("@" + seq_record.id + ':' + pos + " 1\n")
                        out_fh1.write(forward_read + '\n+\n' + read_qual + '\n')
                        # second read
                        sec_read = Util.seq.reverse_complement(mut_seq[next_pos - read_len: next_pos])
                        sec_read_qual = Util.simulate.generate_phred_scores(len(sec_read))
                        out_fh2.write("@" + seq_record.id + ':' + pos + " 2\n")
                        out_fh2.write(sec_read + '\n+\n' + sec_read_qual + '\n')

                    # minus strand
                    frag_prev = frag_index - 1
                    if frag_prev >= 0:
                        (prev_pos, prev_enz) = fragments_ordered[seq_record.id][frag_prev]
                    else:
                        prev_pos = 0
                    if prev_enz == enz_second:
                        # backward
                        back_pos = pos - len(Util.enzyme.enzymeDict().getProperty(enz)['site'])
                        back_start = back_pos - read_len
                        if back_start < 0: back_start = 0
                        backward_read = leftover_seq + Seq(mut_seq[back_start:back_pos]).reverse_complment()
                        read_qual = Util.simulate.generate_phred_scores(len(backward_read))
                        out_fh1.write('@' + seq_record.id + ':' + back_pos + ' 1\n' + backward_read + '\n+\n',
                                      read_qual + '\n')
                        # second read
                        sec_read = mut_seq[prev_pos: prev_pos + read_len]
                        sec_read_qual = Util.simulate.generate_phred_scores(len(sec_read))
                        out_fh2.write("@" + seq_record.id + ':' + pos + " 2\n")
                        out_fh2.write(sec_read + '\n+\n' + sec_read_qual + '\n')

out_fh1.close()
if out_fh2: out_fh2.close()

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