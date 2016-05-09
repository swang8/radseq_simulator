# radseq_simulator
#### Simulate RADseq/GBS reads allowing variations and sequencing errors.
##### Usage:
<pre>
python simulator-main.py -h
usage: simulator-main.py [-h] [-ref REF] [-depth DEPTH] [-len LEN]
                         [-type TYPE] [-e1 E1] [-e2 E2] [-output OUTPUT]

A read simulator for RADseq, introducing variations and sequencing errors.

optional arguments:
  -h, --help      show this help message and exit
  -ref REF        reference file in fasta format
  -depth DEPTH    coverage depth to be simulated
  -len LEN        the length of reads to be generated
  -type TYPE      single-end (S) or pair-end (P), default is single-end
  -e1 E1          The 1st enzyme
  -e2 E2          The 2nd enzyme, optional
  -output OUTPUT  prefix of the output file in fastq format

Example:
python simulator-main.py  -ref test/ecoli.contigs.fasta -depth 100 -len 50 -type P -e1 PstI -e2 MseI -output test/t1
</pre>
