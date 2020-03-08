# import os
# import osproc
import parseopt

proc writeDefaultHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo "CONDUIT Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu, https://github.com/NatPRoach/ )"
  echo "Usage:"
  echo "\t./conduit <nano | hybrid>"

proc writeNanoHelp() =
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo "CONDUIT Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu, https://github.com/NatPRoach/ )"
  echo "Usage:"
  echo "  ./conduit nano [options] <clusters_directory>"
  echo "Options (defaults in parentheses):"
  echo "  Scaffold type:"
  echo "    --drna (default)"
  #echo "    --cdna" #Not supported yet

proc writeHybridHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo "CONDUIT Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu, https://github.com/NatPRoach/conduit/ )"
  echo "Usage:"
  echo "  ./conduit hybrid [options] <clusters_directory> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>}"
  echo "  <clusters_directory>   Directory containing the .fasta/.fa or .fastq/.fq files of reads separated by gene cluster"
  echo "  <m1>                   Files with #1 mates, paired with files in <m2>"
  echo "                         Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)."
  echo "  <m2>                   Files with #2 mates, paired with files in <m1>"
  echo "                         Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)."
  echo "  <r>                    Files with unpaired reads"
  echo "                         Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)."
  echo "  <i>                    File with interleaved paired-end FASTQ/FASTA reads"
  echo "                         Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)."
  echo "  <bam>                  Files are unaligned BAM sorted by read name."
  echo ""
  echo "  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be specified many times."
  echo "  E.g. '-U file1.fq,file2.fq -U file3.fq'."
  echo "Options (defaults in parentheses):"
  echo "  Scaffold type:"
  echo "    --drna (default)"
  echo "        Scaffold reads are stranded, and may contain U characters instead of Ts"
#  echo "    --cdna-fwd-stranded"
#  echo "        Scaffold reads are stranded, coding strand" #Is this possible for cDNA?
  echo "    --cdna-rev-stranded"
  echo "        Scaffold reads are stranded reverse complemented relative to coding strand"
  echo "    --cdna"
  echo "        Scaffold reads are NOT stranded"
  echo "    --sfq (default)"
  echo "        Scaffold reads are in FASTQ format"
  echo "    --sfa"
  echo "        Scaffold reads are in FASTA format"
  echo "  Illumina type:"
  echo "    -U, --unstranded "
  echo "        Illumina reads are unstranded"
  echo "    -f, --fwd_stranded"
  echo "        Illumina reads are stranded s.t. the first mate originates from the RNA strand"
  echo "        Ignored if scaffold reads are not stranded"
  echo "    -r, --rev_stranded"
  echo "        Illumina reads are stranded s.t. the first mate is the reverse complement of the RNA strand"
  echo "        Ignored if scaffold reads are not stranded"
  echo "    --ifq (default)"
  echo "        Illumina reads are in FASTQ format"
  echo "    --ifa"
  echo "        Illumina reads are in FASTA format"
  echo "  "



for kind, key, val in getopt():
  echo kind, " ", key, " ", val
writeDefaultHelp()
echo ""
writeNanoHelp()
echo ""
writeHybridHelp()

###
### Use walkFiles() or walkPattern from os to iterate through the .fasta / .fa files. 
###