#!/usr/bin/env python

from Bio import SeqIO
import re
from collections import defaultdict
import sys

if len(sys.argv) != 2:
    sys.exit()
ref_path = sys.argv[1]

fasta_file = ref_path + '/hg19-tRNAs.fa'
# remove duplicates and intron
out_seq = defaultdict(set)
out_id = set()
with open(fasta_file, 'r') as fasta:
    for record in SeqIO.parse(fasta, 'fasta'):
        sequence = str(record.seq)
        sequence = re.sub('[actgn]','', sequence)
        name = record.id
        print '>%s\n%sCCAA' %(name, sequence)
