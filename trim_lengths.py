#!/usr/bin/env python

infile = open("pass.fa")
outfile = open("pass_long_reads.fa",'w')

read_id = ""

for line in infile:
    if line[0] == '>':
        read_id = line
    else:
        if len(line.strip()) > 100:
            outfile.write(read_id)
            outfile.write(line)
outfile.close()

infile = open("consensus_reads.fa")
outfile = open("consensus_long_reads.fa",'w')

read_id = ""

for line in infile:
    if line[0] == '>':
        read_id = line
    else:
        if len(line.strip()) > 50:
            outfile.write(read_id)
            outfile.write(line)
outfile.close()