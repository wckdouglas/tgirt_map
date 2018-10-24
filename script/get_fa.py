#!/usr/bin/env python

import pysam
from sequencing_tools.fastq_tools import reverse_complement
from xopen import xopen
import re
import fileinput
import sys


def get_sequences(pysam_genome_fa, bed_line, out_bed_handle, out_fa_handle):
    fields = bed_line.strip().split('\t')
    strand = fields[5]
    chrom = fields[0]
    start, end = int(fields[1]), int(fields[2])
    gname = fields[3]
    gtype = fields[6]

    seq = pysam_genome_fa.fetch(chrom, start, end).upper()
    seq = seq if strand == '+' else reverse_complement(seq)
    gtype = 'Mt_protein' if 'protein_coding' in gtype and gname.startswith('Mt') else gtype
    
    print('>{name}\n{seq}'.format(name = gname,seq = seq), file = out_fa_handle)

    bed_line = '{name}\t0\t{gene_end}\t{name}\t0\t+\t{gene_type}\t{gene_id}'
    print(bed_line.format(name = gname,
                            gene_end = len(seq),
                            gene_type = gtype,
                            gene_id = fields[-1]),
            file = out_bed_handle)


def main():
    genome_fa = pysam.Fastafile(sys.argv[1])
    out_bed = open(sys.argv[2], 'w')
    out_fasta  = open(sys.argv[3], 'w')

    for line in sys.stdin:
        get_sequences(genome_fa, line, out_bed, out_fasta)

    
    out_bed.close()
    out_fasta.close()


if __name__=='__main__':
    if len(sys.argv) != 4:
        sys.exit('[usage] python %s <genome_fa> <out_bed_file> <out_fasta>' %sys.argv[0] )
    main()
