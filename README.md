![image](CONDUIT.png)

# CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly
#### Builds a transcriptome independent of a reference genome using Oxford Nanopore Technologies and optionally Illumina RNA-seq data (recommended)

Usage statement for hybrid (ONT + Illumina) assembly:
```
CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:
CONDUIT Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu, https://github.com/NatPRoach/conduit/ )
Usage:
  ./conduit hybrid [options] <clusters_directory> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>}
  <clusters_directory>   Directory containing the .fasta/.fa or .fastq/.fq files of reads separated by gene cluster

  Illumina data is aligned with Bowtie2, therefore Illumina data is provided in the same format as Bowtie2, namely:

    <m1>                   Files with #1 mates, paired with files in <m2>
                           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
    <m2>                   Files with #2 mates, paired with files in <m1>
                           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
    <r>                    Files with unpaired reads
                           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
    <i>                    File with interleaved paired-end FASTQ/FASTA reads
                           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
    <bam>                  Files are unaligned BAM sorted by read name.

  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be specified many times.
  E.g. '-U file1.fq,file2.fq -U file3.fq'.

Options (defaults in parentheses):
  Scaffold Type:
    --drna (default)
        Scaffold reads are stranded forward relative to coding strand, and may contain U characters instead of Ts
    --cdna-rev-stranded
        Scaffold reads are stranded reverse complemented relative to coding strand
    --cdna
        Scaffold reads are NOT stranded
    --sfq (default)
        Scaffold reads are in FASTQ format
    --sfa
        Scaffold reads are in FASTA format
  Illumina Type:
    -u, --unstranded
        Illumina reads are unstranded
    -f, --fwd-stranded
        Illumina reads are stranded s.t. the first mate originates from the RNA strand
        Ignored if scaffold reads are not stranded
    -r, --rev-stranded (default)
        Illumina reads are stranded s.t. the first mate is the reverse complement of the RNA strand
        Ignored if scaffold reads are not stranded
    --ifq (default)
        Illumina reads are in FASTQ format; Mutually exclusive with --ifa
    --ifa
        Illumina reads are in FASTA format; Mutually exclusive with --ifq
  Consensus Collapsing:
    -d, --isoform-delta (35)
        Maximum indel size to be 'corrected', beyond this size a new isoform is declared. Must be between 0 and 255
    -e, --ends-delta (35)
        Maximum size at the ends of isoforms to 'correct' before splitting. Must be between 0 and 255
    -i, --max-iterations (5)
        Maximum number of iterations to align to and correct scaffolds. Does not include optional final polshing step
    --final-polish (default)
        Include a final correction of individual isoforms, not in a splice graph
    --no-final-polish
        Do not do a final correction of individual isoforms not in a splice graph
  Ouput:
    -o, --output-dir <path> (conduit/)
        <path> where corrected clusters will be written
        NOTE: THIS WILL OVERWRITE EXISTING FILES!
    -n, --no-intermediates (default)
        Does not save FASTA file generated for intermediate rounds of polishing
    -s, --save-intermediates
        Saves the FASTA file generated for intermediate rounds of polishing
  Bowtie2:
    --end-to-end (default)
        Align Illumina reads to ONT scaffolds in end-to-end alignment mode
    --local
        Align Illumina reads to ONT scaffolds in local alignment mode
  Miscellaneous:
    -h, --help
        Display this help message and exit
    -v, --version
        Display the installed version number of CONDUIT and exit
    --tmp-dir <path> (conduit-tmp/)
        <path> where temporary files will be created
    -t, --threads (4)
        Number of threads to run in parallel (used for both Bowtie2 and Partial Order Graph correction)
```