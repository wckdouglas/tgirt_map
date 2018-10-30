#!/usr/bin/env python

import sys
from Bio import SeqIO 


amino_acid = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
          'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
          'SEC': 'U', 'SUP':'SUP', 'UNDET':'UND', 'IMET': 'IMET',
          'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

for record in SeqIO.parse(sys.stdin, 'fasta'):
    seq = str(record.seq).replace('U','T')
    name = record.id.replace('Homo_sapiens_tRNA-','')
    name = name.upper()
    fields = name.split('-')
    numbers = '-'.join(fields[2:])
    anticodon = fields[1]
    aa = amino_acid[fields[0]]
    name = 'TR{aa}-{anticodon}-{numbers}'\
        .format(aa = aa,
                anticodon = anticodon,
                numbers = numbers)
    print('>{name}\n{seq}CCAA'.format(name = name, seq = seq), file = sys.stdout)
