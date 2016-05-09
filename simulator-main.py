import Util.enzyme
import Util.seq
import Util.digest
import Util.fasta
import Util.simulate
import sys
import argparse
import datetime
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np

## parse command line parameters
parser = argparse.ArgumentParser(description='A read simulator for RADseq, introducing variations and sequencing errors.')
parser.add_argument('-ref', help='reference file in fasta format', type=str)
parser.add_argument('-depth', help='coverage depth to be simulated', type=int)
parser.add_argument('-len', help='the length of reads to be generated', type=int)
parser.add_argument('-type', help='single-end (S) or pair-end (P), default is single-end', type=str, default='S')
parser.add_argument('-e1', help='The 1st enzyme', type=str)
parser.add_argument('-e2', help='The 2nd enzyme, optional', type=str)
parser.add_argument('-output', help='prefix of the output file in fastq format', type=str)

if len(sys.argv) < 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

if args.type not in ('s', 'S', 'p', 'P'):
    print 'Use S for single end; P for pair end'
    print 'Not recognizing type: ' + args.type
    print 'Will generate single end reads!'

(ref_fasta, depth, read_len, read_type, enz_first, enz_second, output_prefix) = \
    (args.ref, args.depth, args.len, args.type, args.e1, args.e2, args.output)

##read_len += 10
## check the enzyme(s)

if not (enz_first and Util.enzyme.enzymeDict().checkEnz(enz_first)):
    print '\nThe first Enzyme not provided or not in the enzyme dictionary!!\n'
    #print Util.enzyme.enzymeDict().ids()
    parser.print_help()
    exit(1)

if enz_second and not Util.enzyme.enzymeDict().checkEnz(enz_second):
    print '\nEnzyme ' + enz_second + ' not in the enzyme dictionary!!\n'
    #print Util.enzyme.enzymeDict().ids()
    parser.print_help()
    exit(1)


## helpers
def get_re_leftover(enzs):
    leftover = {}
    for enz in enzs:
        if not enz:
            continue
        site = Util.enzyme.enzymeDict().getProperty(enz)['site']
        ovhg = Util.enzyme.enzymeDict().getProperty(enz)['ovhgseq']
        leftover[enz] = site[site.index(ovhg):]
    return leftover

def check_read_length(read, read_length):
    if len(read) < read_length:
        read = read + 'N' * (read_length - len(read))
    else:
        read = read[0:read_length]
    return read

def print_current_time(s):
    print datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), "\t", s

## get the fragments
digest = Util.digest.Digest([enz_first], ref_fasta)
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

## generate fastq file

out_fh1 = open (output_prefix + "_R1.fastq", 'w')
out_fh2 = ''
if read_type in ['p', 'P']: out_fh2 = open (output_prefix + "_R2.fastq", 'w')

enz_leftover =get_re_leftover([enz_first, enz_second])

for seq_record in SeqIO.parse(ref_fasta, "fasta"):
    #print(seq_record.id)
    #print(repr(seq_record.seq))
    #print(len(seq_record))
    if seq_record.id in fragments_ordered:
        (mut_seq, mut_pos) = Util.simulate.mutate_genome(seq_record.id, seq_record.seq, output_prefix)

        #poisson distribution of coverage
        dep_poisson = np.random.poisson(depth, len(fragments_ordered[seq_record.id]))

        for frag_index in range(len(fragments_ordered[seq_record.id])):
            frag = fragments_ordered[seq_record.id][frag_index]
            (pos, enz) = frag
            #print 'pos, enz, enz_first:', pos, enz, enz_first
            # generate reads
            if read_type in ['s', 'S']:
                ## single end
                if enz == enz_first:
                    print 'enz, enz_first:', enz, enz_first
                    leftover_seq = enz_leftover[enz]
                    for d in range(dep_poisson[frag_index]):
                        # forward
                        forward_read = leftover_seq + \
                            Util.simulate.simulate_error(mut_seq[pos: read_len + pos - len(leftover_seq)])
                        forward_read = check_read_length(forward_read, read_len)
                        read_qual = Util.simulate.generate_phred_scores(len(forward_read))
                        out_fh1.write("@" + seq_record.id + ':' + str(pos) + " 1\n")
                        out_fh1.write(forward_read + '\n+\n' +read_qual + '\n')
                        # backward
                        back_pos = pos - len(Util.enzyme.enzymeDict().getProperty(enz)['site'])
                        back_start = back_pos - read_len
                        if back_start < 0:
                            back_start = 0
                        backward_read = leftover_seq + \
                            Util.simulate.simulate_error(str(Seq(mut_seq[back_start:back_pos]).reverse_complement()))
                        backward_read = check_read_length(backward_read, read_len)
                        read_qual = Util.simulate.generate_phred_scores(len(backward_read))
                        read = '@' + seq_record.id + ':' + str(back_pos) + ' 1\n' + backward_read + '\n+\n' + read_qual + '\n'
                        print type(read)
                        out_fh1.write(str(read))

            if read_type in ['p', 'P']:
                ## pair end
                if enz == enz_first:
                    for d in range(dep_poisson[frag_index]):
                        # plus strand
                        frag_next = frag_index + 1
                        (next_pos, next_enz) = ('', '')
                        if frag_next < len(fragments_ordered[seq_record.id]):
                            (next_pos, next_enz) = fragments_ordered[seq_record.id][frag_next]
                        else:
                            next_pos = len(mut_seq) - 1
                        if next_enz == enz_second:
                            leftover_seq = enz_leftover[enz]
                            # forward
                            forward_read = leftover_seq + \
                                        Util.simulate.simulate_error(mut_seq[pos: read_len + pos - len(leftover_seq)])
                            forward_read = check_read_length(forward_read, read_len)
                            read_qual = Util.simulate.generate_phred_scores(len(forward_read))
                            out_fh1.write(str("@" + seq_record.id + ':' + str(pos) + " 1\n"))
                            out_fh1.write(forward_read + '\n+\n' + read_qual + '\n')
                            # second read
                            sec_read = Util.seq.reverse_complement(mut_seq[next_pos - read_len: next_pos])
                            sec_read = Util.simulate.simulate_error(sec_read)
                            sec_read = check_read_length(sec_read, read_len)
                            sec_read_qual = Util.simulate.generate_phred_scores(len(sec_read))
                            out_fh2.write(str("@" + seq_record.id + ':' + str(pos) + " 2\n"))
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
                            backward_read = leftover_seq + \
                                            Util.simulate.simulate_error(str(Seq(mut_seq[back_start:back_pos]).reverse_complement()))
                            backward_read = check_read_length(backward_read, read_len)
                            read_qual = Util.simulate.generate_phred_scores(len(backward_read))
                            out_fh1.write(str('@' + seq_record.id + ':' + str(back_pos) + ' 1\n' + backward_read + '\n+\n' +
                                          read_qual + '\n'))
                            # second read
                            sec_read = mut_seq[prev_pos: prev_pos + read_len]
                            sec_read = check_read_length(sec_read, read_len)
                            sec_read_qual = Util.simulate.generate_phred_scores(len(sec_read))
                            out_fh2.write(str("@" + seq_record.id + ':' + str(pos) + " 2\n"))
                            out_fh2.write(sec_read + '\n+\n' + sec_read_qual + '\n')

out_fh1.close()
if out_fh2: out_fh2.close()

