#!/usr/bin/env python

import pysam

infile1 = pysam.AlignmentFile("test_illumina_to_nanopore.bam",'rb')
infile2 = pysam.AlignmentFile("N2_OP50_bio1.scd-1.bam",'rb')
outfile = pysam.AlignmentFile("aligned_reads.bam",'wb',template=infile2)
aligned_reads = set()
for read in infile1.fetch():
    # print(read.query_name)
    if read.is_unmapped:
        continue
    aligned_reads.add(read.query_name)

for read in infile2.fetch():
    # print(read.query_name)
    if read.flag & 64:
        if read.query_name + '/1' in aligned_reads or read.query_name in aligned_reads:
            outfile.write(read)
    elif read.flag & 128:
        if read.query_name + '/2' in aligned_reads or read.query_name in aligned_reads:
            outfile.write(read)
outfile.close()
