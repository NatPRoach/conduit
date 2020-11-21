Modified POA by Haibao
======================
New features:

- Allow piping through stdin and output stdout
- Disable debugging info with ``-silent``
- PIR/FASTA output gap symbols changed to ``-``
- Compile flags use ``-O2``


POA INSTALLATION NOTES
======================
September 2001, updated March 2004.

Chris Lee
Dept. of Chemistry & Biochemistry
UCLA


I. COMPILATION

To compile this program, simply type 'make poa'.
This produces an executable for sequence alignment (poa) and also a
linkable library liblpo.a.  The software has been compiled and tested
on LINUX and Mac OS X.


II. RUNNING POA

POA has a variety of command line options.  Running POA without any
arguments will print a list of the possible command line arguments.
POA may be used to construct a PO-MSA, or to analyze a PO-MSA.

A.  Constructing a PO-MSA
-------------------------

1.  Required Input:

i. An Alignment Score Matrix File:

A score matrix file is required, because POA uses it to get the
residue alphabet and indexing.  Even if POA is not being used to
perform multiple sequence alignment, this file must be provided.  Any
basic alignment matrix which may be used with BLAST may be used here.
This file must be the first command line argument without a flag in
order to be interpreted by POA as the score matrix file.  Two example
score matrix files, blosum80.mat and blosum80_trunc.mat, are provided
in this directory.  They includes scores for matching nucleotides, as
well as amino acids.  Header lines may be used to specify gap
parameters, as in the examples:

GAP-PENALTIES=A B C
GAP-TRUNCATION-LENGTH=T
GAP-DECAY-LENGTH=D

means that the gap opening penalty is A; the gap extension penalty is
B until the gap length reaches T; the gap extension penalty decreases
linearly from B to C for gap lengths between T and T+D; and the gap
extension penalty is C for all longer gaps.  An additional line,
"GAP-PENALTIES-X=Ax Bx Cx", can be inserted after the GAP-PENALTIES
line to specify the opening and extension penalties for gaps in the
first sequence relative to the second sequence in an alignment, i.e.,
for asymmetric gap scoring.  Use the -v flag to see what gap penalties
POA is using for a given run.

NOTE: POA is case-sensitive.
In order to distinguish amino acid residues from nucleic acid
residues, POA is case sensitive.  Residues that are uppercase are
interpreted as amino acids, while residues that are lowercase are
interpreted as nucleotides.  POA can handle mixed score matrices
containing both amino acid and nucleotide scores, as long as the
column and row labels of the matrix are case-sensitive.  The
blosum80.mat file is an example.

ii.  A FASTA file, or MSA Files in PO, CLUSTAL or PIR Format
A FASTA file is required only if POA is being used to construct a
new PO-MSA from a list of sequences, or to align a list of sequences to
an already existing PO-MSA (see Analyzing a PO-MSA below).  This FASTA file
should contain sequences to be aligned by POA.  The command line argument
to get POA to accept a FASTA file as input is '-read_fasta FILENAME'.  POA
will interpret FILENAME as the FASTA sequence file.  An example file,
multidom.seq, is provided in this directory.

POA is case-sensitive (see NOTE above).  All residues in the FASTA
file must be uppercase to be interpreted as amino acids by POA, or
lowercase to be interpreted as nucleotides.  To switch the case of all
of the letters in the FASTA file to uppercase, use the '-toupper'
command line argument.  To switch the case of all the letters in the
FASTA file to lowercase, use the '-tolower' command line argument.

POA will also read in a set of MSA files to be aligned (see below).


2.  MSA Construction Options:

i.  Global vs. Local Alignment

POA will build alignments using local or global alignment.  The
default is set to local aligment.  To call global alignment, use the
'-do_global' flag.

ii.  Iterative vs. Progressive Alignment

POA will build the alignment iteratively, aligning sequences and MSAs
in the order they are provided.  It will also align sequences and MSAs
in the order dictated by a guide tree built from a matrix of pairwise
similarity scores.  To call progressive alignment, use the
'-do_progressive' flag.

To provide POA with a set of pairwise similarity scores, use the
'-read_pairscores' flag followed by the name of the text file
containing the list of pairwise similarity scores.  This file should
be a tab-delimited file, where each row contains two sequence names
followed by the pairwise similarity score.  The included file
"multidom.pscore" shows an example.

Example row in pairscore file:
ABL1_HUMAN	MATK_HUMAN	260.0

To quickly compute a set of pairwise similarity scores, run BLAST, and
set the similarity scores to the set of BLAST bitscores.  A simple
BLAST driver/parser is provided as "make_pscores.pl"; you may need to
modify this script for your particular configuration.  If the
'-do_progressive' flag is specified without a corresponding pairscore
file, POA will compute pairwise similarity scores itself, by
performing all pairwise alignments.  This is very slow, and we do not
recommend it.

iii.  Aggressive Fusion:

By default, during the building up of a PO-MSA, if a node i with label
'A' is aligned to a node j with label 'B' that belongs to an align
ring containing a third node k with label 'A', POA simply adds node i
to the j-k align ring.  It is possible to force POA to do aggressive
fusion, so that node i is instead fused to node k.  Use the '-fuse_all'
flag to accomplish this.


3.  MSA Output Formats:

POA can output a PO-MSA in several formats simultaneously, including
CLUSTAL, PIR, and PO.  The PO format is the best format since it
contains all of the information in the PO-MSA.  The other formats
accurately represent the MSA, but since they are RC-MSA formats, they
may lose some of the information in the full PO-MSA.

i.  CLUSTAL format:
This format is the standard CLUSTAL format.  The command line argument
to get the MSA output in this format is '-clustal FILENAME'.

ii.  PIR format:
This format is the standard PIR format, which is like FASTA with a '.'
character representing gaps.  The command line argument to get the MSA
output in this format is '-pir FILENAME'.

iii.  PO format:
This format is the standard PO format.  It is described below in the
section PO format.  The command line argument to get the MSA output in
this format is '-po FILENAME'.

EXAMPLE: Constructing an MSA of Four Protein Sequences

Running POA with the following statement will take the FASTA-formatted
sequences in the multidom.seq file, construct a PO-MSA using the
scoring matrix in the file blosum80.mat, and then output the PO-MSA in
CLUSTAL format to the file multidom.aln.

poa -read_fasta multidom.seq -clustal multidom.aln blosum80.mat


4.  Other Output:

i.  Score Matrix

POA will also print to stdout the score matrix stored in the '.mat'
file.  The command line argument to get POA to do this is
'-printmatrix LETTERSET', where LETTERSET is a string of letters to be
printed with the score matrix.  For example, if the score matrix is
designed for protein alignment the letter set might be
'ARNDCQEGHILKMFPSTWYV'.

ii.  Verbose Mode

POA will run in verbose mode, printing additional information
generated during the run (such as the set of gap scores used) to
stdout.  The command line argument for verbose mode is '-v'.


B.  Analyzing a PO-MSA
-----------------------

POA can also take as input an MSA in PO, CLUSTAL or PIR file format
and rebuild the PO-MSA data structure.  Once this data structure has
been rebuilt, it may be analyzed for features.  In "liblpo.a", the
linkable POA library, we have included the functions necessary to do
heaviest bundling and thereby find consensus sequences in the PO-MSA
(the details of the heaviest bundling algorithm are described
elsewhere).  POA has been written so that users may create their own
functions for analyzing a PO-MSA.  We have not included in the
"liblpo.a" library the functions that we wrote to analyze PO-MSAs
constructed with ESTs and genome sequence to find snps and alternative
splice sites.  However, it is possible to design modular library
functions that will look for highly specific biological features in
any PO-MSA data structure.

1.  Required Input:

Before the PO-MSA data structure can be analyzed it must be built.  It can
be built iteratively or using a guide tree, or converted from another file
type.  POA can align a set of FASTA-formatted sequences to each other or
to an existing PO-MSA.  It can align two PO-MSAs.  It can also align an
arbitrary set of PO-MSAs.

Note:  POA Requires Either An MSA File or a FASTA File
If neither type of file is read in by POA it will terminate early, since it
has not received any sequence data.

i.  An MSA file in PO, CLUSTAL or PIR format:

POA will read in an MSA file in PO, CLUSTAL or PIR format.  The
command line argument to get poa to read in an MSA file and rebuild
the PO-MSA data structure is '-read_msa FILENAME'.  POA automatically
determines whether the MSA file is in PO, CLUSTAL, or PIR format.  POA
will read in a second MSA file when the '-read_msa2 FILENAME' flag is
used.  POA will read in a set of MSA files using the '-read_msa_list
FILENAME' flag.  The file should contain a list of names of MSA files.

It is possible to filter the PO-MSA data structure as it is being
rebuilt.  In order to filter the PO-MSA in the MSA file to include
only a subset of sequences use the command line argument '-subset
FILENAME', where the file named FILENAME contains the list of sequence
names to be included in the new PO-MSA. In order to filter the PO-MSA
in the MSA file to exclude a subset of sequences, use the command line
argument '-remove FILENAME', where the file named FILENAME contains
the list of sequences to be excluded from the new PO-MSA.  The names
of sequences to be included or excluded should be in the format
"SOURCENAME=*", as they are in a PO file.  Lists of sequence source
names can be created by using the unix grep utility on the PO file.
Each line in the list of sequences to be filtered should read,
"SOURCENAME=" followed by the name of the sequence,
e.g. "SOURCENAME=ABL1_HUMAN".  To filter the second PO-MSA read in
using the '-read_msa2 FILENAME', use the '-subset2' and '-remove2'
flags.

ii.  A FASTA File:

The FASTA file should contain sequences to be aligned by POA.  The
command line argument to get POA to accept a FASTA file as input is
'-read_fasta FILENAME'.  POA will interpret FILENAME as the FASTA
sequence file.  An example file, "multidom.seq", is provided in this
directory.  (See note above on case-sensitivity).

NOTE: POA Can Take Both An MSA File And A FASTA File As Input
If both the '-read_msa FILENAME' argument and the '-read_fasta FILENAME'
argument are given to POA on the command line, then POA will first rebuild the
PO-MSA in the MSA file, and then it will align the sequences in the FASTA file
to this PO-MSA.  Similarly, if both the '-read_msa_list FILENAME' flag
and the '-read_fasta FILENAME' flag are given to POA, then POA will
rebuild all of the PO-MSAs and will align them to each other and to the
sequences in the FASTA file.


2.  Additional PO Utilities:

i.  Consensus Generation Via Heaviest Bundling Algorithm:
The heaviest bundling algorithm finds consensus sequences in the
PO-MSA.  The command line argument for heaviest bundling is '-hb'.  This
function adds the new consensus sequences to the PO-MSA by storing new
consensus sequence indices on the in the PO-MSA nodes corresponding to
the consensus sequence paths.  The sequence source names for consensus
sequences generated by heaviest bundling are CONSENS'i' where 'i' is the
index of the bundle corresponding to the consensus sequence.

The heaviest bundling algorithm can also take as input a bundling
threshold value.  The command line argument for setting a bundling
threshold value for heaviest bundling is '-hbmin VALUE'.  This threshold
is used during the process of associating sequences with bundles.  If a
sequence has a percentage of nodes shared with bundle 'i' greater than this
threshold value, it is associated with bundle 'i'.  Iterative heaviest
bundling can also be affected by the bundling threshold.  A detailed
description of heaviest bundling and heaviest bundling thresholds is given
elsewhere.  The consensus sequences corresponding to bundles generated
by heaviest bundling are listed in the sequence source list.  Additionally,
in the SOURCEINFO line for each sequence the index of the bundle to which
that sequence belongs is give.  Finally, using the command line argument
'-best' restricts the MSA output to the consensus sequences generated
by heaviest bundling (NB: this applies to PIR output only).


III.  PO FILE FORMAT

****************************HEADER****************************************
VERSION= ~Current version of POA,e.g. LPO.1.0~
NAME=  ~Name of PO-MSA.  Defaults to name of 1st sequence in PO-MSA~
TITLE=  ~Title of PO-MSA.  Defaults to title of 1st sequence in PO-MSA~
LENGTH=  ~Number of nodes in PO-MSA~
SOURCECOUNT=  ~Number of sequences in PO-MSA~

*********************SEQUENCE SOURCE LIST*********************************

/* For each sequence in the PO-MSA: */
SOURCENAME= ~Name of sequence taken from FASTA sequence header~
SOURCEINFO= ~Number of nodes in sequence~
            ~Index of first node containing sequence~
            ~Sequence weight~
            ~Index of bundle containing sequence~
            ~Title of sequence taken from FASTA sequence header~

/* Example: */
SOURCENAME=GRB2_HUMAN
SOURCEINFO=217 10 0 3 GROWTH FACTOR RECEPTOR-BOUND PROTEIN 2 (GRB2 ADAPTOR PROTEIN)(SH2)

********************PO-MSA DATA STRUCTURE*********************************

/* For each node in the PO-MSA:  */
~Residue label~:~'L' delimited index list of other nodes with edges into node~
                ~'S' delimited index list of sequences stored in each node~
                ~'A' index of next node in same align ring~
                     NB: align ring indices must form a cycle.
                     e.g. if two nodes 121 and 122 are aligned, then
                     the line for node 121 indicates "A122", and
                     the line for node 122 indicates "A121".

/* Example: */
F:L156L155L22S2S3S7A158

********************END***************************************************


For more information, see http://www.bioinformatics.ucla.edu/poa.
