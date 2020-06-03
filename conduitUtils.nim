# import os
# import osproc
import parseopt
import strutils
import strformat
import tables
import sets
import poGraphUtils
import algorithm

type
  BLASTmatch* = object
    query_name : string
    query_len : uint
    match_names : seq[string]
    match_lens : seq[uint]
    scores : seq[seq[float]]
    evals : seq[seq[float]]
    identities_numerators : seq[seq[uint]]
    identities_denominators : seq[seq[uint]]
    identities_percentages : seq[seq[float]]
    positives_numerators : seq[seq[uint]]
    positives_denominators : seq[seq[uint]]
    positives_percentages : seq[seq[float]]
    gap_numerators : seq[seq[uint]]
    gap_denominators : seq[seq[uint]]
    gap_percentages : seq[seq[float]]
    query_seqs : seq[seq[string]]
    conss_seqs : seq[seq[string]]
    sbjct_seqs : seq[seq[string]]
  
  UtilOptions = object
    mode : string
    run_flag : bool
    infilepath : string
    outfilepath : string
    reference_infilepath : string
    min_length : uint64
    fastq : bool
    stranded : bool

  FastqRecord* = object
    read_id* : string
    sequence* : string
    qualities* : string



proc conduitUtilsVersion() : string =
  return "CONDUIT Utilities Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu, https://github.com/NatPRoach/conduit/ )"

proc writeDefaultHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "Usage:"
  echo "  ./conduitUtils <function>"
  echo "Where <function> is one of the following:"
  echo "  translate   - Translates FASTA/Q nucleotide sequences into protein based on their longest ORF"
  echo ""
  echo "  bed2gtf     - Converts BED12 files to well structured GTF file suitable for use in GFFcompare"
  echo ""
  echo "  parseBLASTP   - Parses BLASTP output and outputs closest match for each query transcript as determined by BLASTP"
  echo ""
  echo "  compareBLASTP - Compares BLASTP output and reference proteome to determine the # of TP, FP, and FN for a sample"
  echo ""
  echo "  compareFASTA  - Compares two FASTA files, an input and a reference, to determine the # of TP, FP, and FN for a sample"
  echo ""
  echo "  splitFASTA    - Splits CONDUIT produced FASTA file based on the number of reads supporting each isoform"
  echo ""
  echo "  filterFASTA   - Filters CONDUIT produced FASTA file based on number of reads supporting each isoform"
  echo ""
  echo "  extractIntrons   - Extracts out intronic sequences from BED12 formatted input and outputs as BED6"
  echo ""
  echo "  callNonCanonical - Reads in a FASTA file and reports the readIDs of sequences that dont begin with GT and end with AG"
  echo ""
  echo "  callNovelNonCanonical - Compares introns described by reference GTF file to introns described by a list of readIDs in the format produced by `bedtools getfasta -name` function, outputs the novel introns in BED format"
  echo ""
  echo "  callOverlapping - Compares two files of readIDs specifying introns in the format produced by `bedtools getfasta -name`, and reports the introns that are shared between the two files (not stranded)"
  echo ""


proc writeTranslateHelp() =
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "translate   - Translates FASTA/Q nucleotide sequences into protein based on their longest ORF"
  echo "Usage:"
  echo "  ./conduitUtils translate [options] -i <transcripts.fa> -o <predicted_protein.fa>"
  echo "  <transcripts.fa>         FASTA/Q infile containing putative transcripts to be translated"
  echo "  <predicted_protein.fa>   FASTA outfile containing in silico translated ORFs from transcripts.fa"
  echo ""
  echo "Options (defaults in parentheses):"
  echo "  Input Options:"
  echo "    -a, --fasta (default)"
  echo "        Input file is in FASTA format"
  echo "    -q, --fastq"
  echo "        Input file is in FASTQ format"
  echo "    -s, --stranded"
  echo "        Input reads are forward stranded"
  echo "  Filtering Options:"
  echo "    -l, --min-length (75)"
  echo "        Minimum length in Amino Acids necessary for a putative ORF to be reported"

proc writeStrandTranscriptsHelp() =
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "strandTranscripts  - Translates FASTA/Q nucleotide sequences into protein based on their longest ORF"
  echo "Usage:"
  echo "  ./conduitUtils strandTranscripts [options] -i <transcripts.fa> -o <predicted_protein.fa>"
  echo "  <transcripts.fa>            FASTA/Q infile containing putative transcripts to be translated"
  echo "  <stranded_trasncripts.fa>   FASTA outfile containing in silico translated ORFs from transcripts.fa"
  echo ""
  echo "Options (defaults in parentheses):"
  echo "  Input Options:"
  echo "    -a, --fasta (default)"
  echo "        Input file is in FASTA format"
  echo "    -q, --fastq"
  echo "        Input file is in FASTQ format"


proc writeBED2GTFHelp() =
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "bed2gtf     - Converts BED12 files to well structured GTF file suitable for use in GFFcompare"
  echo "Usage:"
  echo "  ./conduitUtils bed2gtf -i <infile.bed> -o <outfile.gtf>"
  echo "  <infile.bed>    BED12 infile to be converted in to GTF format"
  echo "  <outfile.gtf>   GTF outfile"

proc writeBLASTPHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "parseBLASTP   - Parses BLASTP output and outputs closest match for each query transcript as determined by BLASTP"
  echo "Usage:"
  echo "  ./conduitUtils parseBLASTP -i <inBLASTP.txt> -o <outPutativeOrthologs.tsv>"
  echo "  <inBLASTP.txt>              Default output of BLASTP search of translated protein products vs some reference proteome"
  echo "  <outPutativeOrthologs.tsv>  Tab separated file of putative ortholog matches"
  echo "                              In format: <Query ID>\t<Reference proteome top match ID>\t<E value>"

proc writeCompareBLASTPHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "compareBLASTP - Compares BLASTP output and reference proteome to determine the # of true positives, false positives, and false negatives for a sample"
  echo "Usage:"
  echo "  ./conduitUtils compareBLASTP -r <reference_proteome.fa> -i <inBLASTP.txt>"
  echo "  <reference_proteome.fa>     FASTA file describing the reference proteome used in the BLASTP search"
  echo "  <inBLASTP.txt>              Default output of BLASTP search of translated protein products vs some reference proteome"

proc writeCompareFASTAHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "compareFASTA  - Compares two FASTA files, an input and a reference, to determine the # of true positives, false positives, and false negatives for a sample"
  echo "Usage:"
  echo "  ./conduitUtils compareFASTA -r <reference.fa> -i <query.fa>"
  echo "  <reference.fa>              Reference FASTA file defining the truth set"
  echo "  <query.fa>                  Query FASTA files defining the query set"

proc writeSplitFASTAHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "splitFASTA    - Splits CONDUIT produced FASTA file based on the number of reads supporting each isoform"
  echo "Usage:"
  echo "  ./conduitUtils splitFASTA -i <conduit_output.fa> -o <outprefix>"
  echo "  <conduit_output.fa>         CONDUIT produced FASTA file to be split based on number of reads supporting each isoform"
  echo "  <outprefix>                 Prefix for the fasta files to be output, suffix will describe the bin being reported"

proc writeFilterFASTAHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "filterFASTA   - Filters CONDUIT produced FASTA file based on number of reads supporting each isoform"
  echo "Usage:"
  echo "  ./conduitUtils filterFASTA -i <inBLASTP.txt> -o <outPutativeOrthologs.tsv>"
  echo "  <conduit_output.fa>         CONDUIT produced FASTA file to be filtered based on number of reads supporting each isoform"
  echo "  <filtered.fa>               Output FASTA file for filtered reads"
  echo "Options: (defaults in parentheses)"
  echo "  Filtering options:"
  echo "     -n (5)"
  echo "        Minimum number of reads that must support an isoform for it to be reported in the filtered FASTA"

proc writeExtractIntronsHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "extractIntrons   - Extracts out intronic sequences from BED12 formatted input and outputs as BED6"
  echo "Usage:"
  echo "  ./conduitUtils extractIntrons -i <transcripts.bed12> -o <introns.bed>"
  echo "  <transcripts.bed12>         Transcripts in BED12 format to extract introns from"
  echo "  <introns.bed>               BED6 output of extracted introns"

proc writeCallNonCanonicalHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "callNonCanonical - Reads in a FASTA file and reports the readIDs of sequences that dont begin with GT and end with AG"
  echo "Usage:"
  echo "  ./conduitUtils callNonCanonical -i <introns.fa> -o <noncanonical.txt>"
  echo "  <introns.fa>              FASTA describing the stranded sequence of introns extracted from `extractIntrons`"
  echo "                            Introns sequences can be obtained using `bedtools getfasta -name -s`"
  echo "  <noncanonical.txt>        Read IDs of the sequences that didn't begin with GT and end with AG"

proc writeCallNovelNonCanonicalHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "callNovelNonCanonical - Compares introns described by reference GTF file to introns described by a list of readIDs in the format produced by `bedtools getfasta -name` function, outputs the novel introns in BED format"
  echo "Usage:"
  echo "  ./conduitUtils callNovelNonCanonical -r <reference.gtf> -i <noncanonical.txt> -o <novel.bed>"
  echo "  <reference.gtf>              Reference GTF file specifying the introns to compare against"
  echo "  <noncanonical.txt>           Read IDs specifying intron structure in the format produced by `bedtools getfasta -name`"
  echo "  <novel.bed>                  Output of introns found in the noncanonical.txt file but not found in the reference, in BED6 format"

proc writeCallOverlappingHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "callOverlapping - Compares two files of readIDs specifying introns in the format produced by `bedtools getfasta -name`, and reports the introns that are shared between the two files (not stranded)"
  echo "Usage:"
  echo "  ./conduitUtils callOverlapping -r <introns1.txt> -i <introns2.txt> -o <shared_introns.txt>"
  echo "  <introns1.txt>              Read IDs specifying introns in the format produced by `bedtools getfasta -name`"
  echo "  <introns2.txt>              Read IDs specifying introns in the format produced by `bedtools getfasta -name`"
  echo "  <shared_introns.txt>        The introns in common between the two files"


proc parseFASTQ(infile : File) : seq[FastqRecord] = 
  while true:
    try:
      let line1 = infile.readLine()
      try:
        assert line1[0] == '@'
      except AssertionError:
        echo "File not in FASTQ format"
        raise
      let read_id = line1[1..^1]
      let sequence = infile.readLine().replace(sub='U',by='T')
      discard infile.readLine()
      let quals = infile.readLine()
      result.add(FastqRecord( read_id   : read_id,
                              sequence  : sequence,
                              qualities : quals))
    except EOFError:
      break

proc convertFASTQtoFASTA(record : FastqRecord) : FastaRecord = 
  result.read_id = record.read_id
  result.sequence = record.sequence

proc convertFASTQtoFASTA(records : openArray[FastqRecord]) : seq[FastaRecord] =
  for record in records:
    result.add(convertFASTQtoFASTA(record))

proc parseBLASTPoutput*(infilepath : string) : seq[BLASTmatch] = 
  ##
  ## Takes in BLASTP default output file and outputs seq of matching proteins
  ##
  
  var infile : File
  discard open(infile,infilepath,fmRead)
  # Discard the header / citation info
  var query_line : string
  var first_iter1 = true
  while true:
    query_line = infile.readLine()
    if query_line.len >= 7:
      if query_line[0..6] == "Query= ":
        break
  while true:
    try:
      if not first_iter1:
        query_line = infile.readLine() # Query= cluster_10_0 <unknown description>
      else:
        first_iter1 = false
      # echo query_line
      if query_line[0..11] == "  Database: ":
        break
      assert query_line[0..6] == "Query= "
      var query_name = query_line[7..^1]
      while true:
        let next = infile.readLine() # ""
        if next == "":
          break
        else:
          query_name = query_name & next

      var match_names : seq[string]
      var match_lens  : seq[uint]
      var evals : seq[seq[float]]
      var scores : seq[seq[float]]
      var identities_numerators : seq[seq[uint]]
      var identities_denominators : seq[seq[uint]]
      var identities_percentages : seq[seq[float]]
      var positives_numerators : seq[seq[uint]]
      var positives_denominators : seq[seq[uint]]
      var positives_percentages : seq[seq[float]]
      var gap_numerators : seq[seq[uint]]
      var gap_denominators : seq[seq[uint]]
      var gap_percentages : seq[seq[float]]
      var query_seqs : seq[seq[string]]
      var conss_seqs : seq[seq[string]]
      var sbjct_seqs : seq[seq[string]]
      

      let query_length_line = infile.readLine() # "Length=229"
      # echo query_length_line
      assert query_length_line[0..6] == "Length="
      let query_length = parseUInt(query_length_line[7..^1])
      # echo query_length
      let splitline = infile.readLine()
      if splitline == "":
        # No significant alignments here
        discard infile.readLine() # ""
        assert infile.readLine() == "***** No hits found *****"
        discard infile.readLine() # ""
        discard infile.readLine() # ""
        discard infile.readLine() # ""
      elif splitline == "                                                                      Score        E":
        assert infile.readLine() == "Sequences producing significant alignments:                          (Bits)     Value"
        discard infile.readLine() # ""
        var nextline : string
        var num_seqs = 0
        while true:
          nextline =  infile.readLine()
          if nextline.len == 0:
            continue
          if nextline[0] == '>':
            break
          num_seqs += 1
        # echo num_seqs
        for i in 0..<num_seqs:
          # echo "seq num: ", i
          # echo nextline
          var match_name = nextline[2..^1]
          var first_iter = true
          while true:
            if nextline.len == 0:
              nextline = infile.readLine()
              continue
            # echo "nextline: ", nextline
            if nextline.len >= 7:
              if nextline[0..6] == "Length=":
                break
            if first_iter:
              first_iter = false
            else:
              match_name = match_name & nextline
            nextline = infile.readLine()
          match_names.add(match_name)
          # echo "match name - ", match_name
          let match_length = parseUInt(nextline[7..^1])
          # echo "match len - ", match_length
          match_lens.add(match_length)
          var scores_for_match : seq[float]
          var evals_for_match : seq[float]
          var ids_numerators_for_match : seq[uint]
          var ids_denominators_for_match : seq[uint]
          var ids_percentages_for_match : seq[float]
          var pos_numerators_for_match : seq[uint]
          var pos_denominators_for_match : seq[uint]
          var pos_percentages_for_match : seq[float]
          var gap_numerators_for_match : seq[uint]
          var gap_denominators_for_match : seq[uint]
          var gap_percentages_for_match : seq[float]
          var query_seqs_for_match : seq[string]
          var conss_seqs_for_match : seq[string]
          var sbjct_seqs_for_match : seq[string]
          discard infile.readLine()
          while true:
            nextline = infile.readLine().strip()
            if nextline.len == 0:
              break
            elif nextline[0..7] == "Score = ":
              var fields0 = nextline.split(seps={','})
              var fields1 = fields0[0].splitWhitespace()
              let score = parseFloat(fields1[2])
              scores_for_match.add(score)
              fields1 = fields0[1].strip().splitWhitespace()
              let evalue = parseFloat(fields1[2])
              evals_for_match.add(evalue)
              fields1 = fields0[2].strip().split(seps={':'})
              let blast_method = fields1[1].strip()
              nextline = infile.readLine().strip()
              fields0 = nextline.split(seps={','})
              fields1 = fields0[0].strip().splitWhitespace()
              var fields2 = fields1[2].strip().split(seps={'/'})
              let ids_numerator = parseUInt(fields2[0])
              let ids_denominator = parseUInt(fields2[1])
              let ids_percentage = 100.0 * float(ids_numerator) / float(ids_denominator)
              
              ids_numerators_for_match.add(ids_numerator)
              ids_denominators_for_match.add(ids_denominator)
              ids_percentages_for_match.add(ids_percentage)

              fields1 = fields0[1].strip().splitWhitespace()
              fields2 = fields1[2].strip().split(seps={'/'})
              let pos_numerator = parseUInt(fields2[0])
              let pos_denominator = parseUInt(fields2[1])
              let pos_percentage = 100.0 * float(pos_numerator) / float(pos_denominator)

              pos_numerators_for_match.add(pos_numerator)
              pos_denominators_for_match.add(pos_denominator)
              pos_percentages_for_match.add(pos_percentage)

              fields1 = fields0[2].strip().splitWhitespace()
              fields2 = fields1[2].strip().split(seps={'/'})
              let gap_numerator = parseUInt(fields2[0])
              let gap_denominator = parseUInt(fields2[1])
              let gap_percentage = 100.0 * float(gap_numerator) / float(gap_denominator)

              gap_numerators_for_match.add(gap_numerator)
              gap_denominators_for_match.add(gap_denominator)
              gap_percentages_for_match.add(gap_percentage)

              discard infile.readLine() # ""
              var query_seq : string
              var sbjct_seq : string
              var conss_seq : string
              var query_start_idx : uint
              var sbjct_start_idx : uint
              var first_iter = true
              while true:
                nextline = infile.readLine()
                if nextline == "":
                  break
                let queryline = nextline
                let consensusline = infile.readLine()
                let sbjctline = infile.readLine()
                discard infile.readLine()
                let queryfields = queryline.strip().splitWhitespace()
                let conss_subseq = consensusline[12..^1] # Have to do it this way because sometimes there'll be a space at the beginning
                let sbjctfields = sbjctline.strip().splitWhitespace()
                # echo queryfields
                if first_iter:
                  first_iter = false
                  query_start_idx = parseUInt(queryfields[1])
                  sbjct_start_idx = parseUInt(sbjctfields[1])
                let query_subseq = queryfields[2]
                let sbjct_subseq = sbjctfields[2]
                query_seq = query_seq & query_subseq
                conss_seq = conss_seq & conss_subseq
                sbjct_seq = sbjct_seq & sbjct_subseq
              # echo query_seq
              # echo conss_seq
              # echo sbjct_seq
              query_seqs_for_match.add(query_seq)
              conss_seqs_for_match.add(conss_seq)
              sbjct_seqs_for_match.add(sbjct_seq)
              # nextline = infile.readLine()
            else:
              break
          
          evals.add(evals_for_match)
          scores.add(scores_for_match)

          identities_numerators.add(ids_numerators_for_match)
          identities_denominators.add(ids_denominators_for_match)
          identities_percentages.add(ids_percentages_for_match)

          positives_numerators.add(pos_numerators_for_match)
          positives_denominators.add(pos_denominators_for_match)
          positives_percentages.add(pos_percentages_for_match)

          gap_numerators.add(gap_numerators_for_match)
          gap_denominators.add(gap_denominators_for_match)
          gap_percentages.add(gap_percentages_for_match)

          query_seqs.add(query_seqs_for_match)
          conss_seqs.add(conss_seqs_for_match)
          sbjct_seqs.add(sbjct_seqs_for_match)
        # assert nextline == ""
      else:
        echo "ERROR PARSING BLASTP OUTPUT"
        break
      # discard infile.readLine() # ""
      # let test = infile.readLine()
      # echo test
      # assert test == "Lambda      K        H        a         alpha"
      assert infile.readLine() == "Lambda      K        H        a         alpha"
      discard infile.readLine() #    0.308    0.126    0.365    0.792     4.96
      discard infile.readLine() # ""
      assert infile.readLine() == "Gapped"
      assert infile.readLine() == "Lambda      K        H        a         alpha    sigma"
      discard infile.readLine() #    0.267   0.0410    0.140     1.90     42.6     43.6
      discard infile.readLine() # ""
      discard infile.readLine() # "Effective search space used: 318188442"
      discard infile.readLine() # ""
      discard infile.readLine() # ""
      result.add(BLASTmatch( query_name : query_name,
                             query_len : query_length,
                             match_names : match_names,
                             match_lens : match_lens,
                             evals : evals,
                             scores : scores,
                             query_seqs : query_seqs,
                             conss_seqs : conss_seqs,
                             sbjct_seqs : sbjct_seqs,
                             identities_numerators : identities_numerators,
                             identities_denominators : identities_denominators,
                             identities_percentages : identities_percentages,
                             positives_numerators : positives_numerators,
                             positives_denominators : positives_denominators,
                             positives_percentages : positives_percentages,
                             gap_numerators : gap_numerators,
                             gap_denominators : gap_denominators,
                             gap_percentages : gap_percentages ))
    except EOFError:
      break
  infile.close()
  
proc translateORF*(nts : string,to_stop = true) : string =
  let translation_table = { "TTT" : 'F',
                            "TTC" : 'F',
                            "TTA" : 'L',
                            "TTG" : 'L',
                            "CTT" : 'L',
                            "CTC" : 'L',
                            "CTA" : 'L',
                            "CTG" : 'L',
                            "ATT" : 'I',
                            "ATC" : 'I',
                            "ATA" : 'I',
                            "ATG" : 'M',
                            "GTT" : 'V',
                            "GTC" : 'V',
                            "GTA" : 'V',
                            "GTG" : 'V',
                            "TCT" : 'S',
                            "TCC" : 'S',
                            "TCA" : 'S',
                            "TCG" : 'S',
                            "CCT" : 'P',
                            "CCC" : 'P',
                            "CCA" : 'P',
                            "CCG" : 'P',
                            "ACT" : 'T',
                            "ACC" : 'T',
                            "ACA" : 'T',
                            "ACG" : 'T',
                            "GCT" : 'A',
                            "GCC" : 'A',
                            "GCA" : 'A',
                            "GCG" : 'A',
                            "TAT" : 'Y',
                            "TAC" : 'Y',
                            "TAA" : '*',
                            "TAG" : '*',
                            "CAT" : 'H',
                            "CAC" : 'H',
                            "CAA" : 'Q',
                            "CAG" : 'Q',
                            "AAT" : 'N',
                            "AAC" : 'N',
                            "AAA" : 'K',
                            "AAG" : 'K',
                            "GAT" : 'D',
                            "GAC" : 'D',
                            "GAA" : 'E',
                            "GAG" : 'E',
                            "TGT" : 'C',
                            "TGC" : 'C',
                            "TGA" : '*',
                            "TGG" : 'W',
                            "CGT" : 'R',
                            "CGC" : 'R',
                            "CGA" : 'R',
                            "CGG" : 'R',
                            "AGT" : 'S',
                            "AGC" : 'S',
                            "AGA" : 'R',
                            "AGG" : 'R',
                            "GGT" : 'G',
                            "GGC" : 'G',
                            "GGA" : 'G',
                            "GGG" : 'G'}.toTable()
  var translation : seq[char]
  # assert nts.len mod 3 == 0
  for i in 0..<(nts.len div 3):
    let aa = translation_table[nts[i*3..((i+1)*3) - 1]]
    if aa == '*':
      result = translation.join("")
      break
    translation.add(aa)

proc findAll*(s,sub : string) : seq[int] = 
  #Very stupid means of finding ALL the matches of sub in s
  var start = 0
  while true:
    let idx = s.find(sub,start = start)
    if idx == -1:
      break
    result.add(idx)
    start = idx + 1 

proc revComp*(nts : string) : string =
  let wc_pairs = {'A' : 'T',
                  'C' : 'G',
                  'T' : 'A',
                  'G' : 'C'}.toTable()
  var revcomp : seq[char]
  for i in 1..nts.len:
    revcomp.add(wc_pairs[nts[^i]])
  result = revcomp.join("")

proc translateTranscript*(nts : string) : string =
  let start_codon_indices = findAll(nts,"ATG")
  result = ""
  for start_codon_index in start_codon_indices:
    let translation = translateORF(nts[start_codon_index..^1])
    if translation.len > result.len:
      result = translation

proc translateTranscripts*(transcripts : openArray[FastaRecord],outfilepath : string , threshold : int = 75,wrap_len : int = 60, stranded : bool = false) =
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  for transcript in transcripts:
    let upper = transcript.sequence.toUpperAscii
    if upper.find('N') != -1:
      continue
    var translation : string
    if stranded:
      translation = translateTranscript(upper)
    else:
      let translation1 = translateTranscript(upper)
      let translation2 = translateTranscript(revComp(upper))
      if translation1.len > translation2.len:
        translation = translation1
      else:
        translation = translation2
    if translation.len >= threshold:
      outfile.write(&">{transcript.read_id}\n")
      for i in 0..<(translation.len div wrap_len):
        outfile.write(&"{translation[i*wrap_len..(i+1)*wrap_len - 1]}\n")
      if translation.len mod wrap_len != 0 :
        outfile.write(&"{translation[wrap_len*(translation.len div wrap_len)..^1]}\n")
  outfile.close()

proc translateTranscripts*(transcripts : openArray[FastqRecord],outfilepath : string , threshold : int = 75,wrap_len : int = 60,stranded = false) =
  translateTranscripts(convertFASTQtoFASTA(transcripts),outfilepath,threshold,wrap_len,stranded)

proc convertBED12toGTF*(infilepath : string, outfilepath : string ,stranded : bool = false ) =
  ## Converts BED12 formatted file to well-formed GTF file suitable for evaluation with GFFcompare
  var infile,outfile : File
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  try:
    while true:
      let bedline = infile.readLine()
      let bedfields = bedline.split(sep='\t')
      let chr = bedfields[0]
      let start_idx = parseUInt(bedfields[1]) + 1'u
      let end_idx = parseUInt(bedfields[2])
      let txid = bedfields[3]
      let strand = bedfields[5]
      if stranded:
        outfile.write(&"{chr}\tBLANK\ttranscript\t{start_idx}\t{end_idx}\t.\t{strand}\t.\ttranscript_id \"{txid}\";\n")
      else:
        outfile.write(&"{chr}\tBLANK\ttranscript\t{start_idx}\t{end_idx}\t.\t.\t.\ttranscript_id \"{txid}\";\n")
      let block_sizes = bedfields[10].split(sep=',')
      let block_starts= bedfields[11].split(sep=',')

      # var exon_count = 1
      # var reverse_exon_count = block_starts.len
      for i in 0..<block_starts.len:
        var new_start_idx, new_end_idx, exon_number : uint
        if strand == "+" or not stranded:
          new_start_idx = start_idx + parseUInt(block_starts[i])
          new_end_idx = start_idx + parseUInt(block_starts[i]) + parseUInt(block_sizes[i]) - 1
          exon_number = uint(i + 1)
        elif strand == "-":
          new_start_idx = start_idx + parseUInt(block_starts[i])
          new_end_idx = start_idx + parseUInt(block_starts[i]) + parseUInt(block_sizes[i]) - 1
          exon_number = uint(block_starts.len - i)
        if stranded:
          outfile.write(&"{chr}\tBLANK\texon\t{new_start_idx}\t{new_end_idx}\t.\t{strand}\t.\ttranscript_id \"{txid}\"; exon_number \"{exon_number}\"\n")
        else:
          outfile.write(&"{chr}\tBLANK\texon\t{new_start_idx}\t{new_end_idx}\t.\t.\t.\ttranscript_id \"{txid}\"; exon_number \"{exon_number}\"\n")
  except EOFError:
    discard
  infile.close()
  outfile.close()

proc strandTranscripts*(transcripts : openArray[FastaRecord],outfilepath : string, wrap_len : int = 60) =
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  for transcript in transcripts:
    let upper = transcript.sequence.toUpperAscii
    if upper.find('N') != -1:
      continue
    let rev = revComp(upper)
    var stranded : string
    let translation1 = translateTranscript(upper)
    let translation2 = translateTranscript(rev)
    if translation1.len > translation2.len:
      stranded = upper
    else:
      stranded = rev
    outfile.write(&">{transcript.read_id}\n")
    for i in 0..<(stranded.len div wrap_len):
      outfile.write(&"{stranded[i*wrap_len..(i+1)*wrap_len - 1]}\n")
    if stranded.len mod wrap_len != 0 :
      outfile.write(&"{stranded[wrap_len*(stranded.len div wrap_len)..^1]}\n")
  outfile.close()


proc strandTranscripts*(transcripts : openArray[FastqRecord],outfilepath : string, wrap_len : int = 60) =
  translateTranscripts(convertFASTQtoFASTA(transcripts),outfilepath,wrap_len)

proc compareExactTranslations*(reference_infilepath : string, translation_infilepath : string) =
  var r_infile, t_infile : File
  discard open(r_infile,reference_infilepath,fmRead)
  discard open(t_infile,translation_infilepath,fmRead)
  let r_records = poGraphUtils.parseFasta(r_infile)
  r_infile.close()
  let t_records = poGraphUtils.parseFasta(t_infile)
  t_infile.close()
  var r_proteins,t_proteins : HashSet[string]
  for record in r_records:
    r_proteins.incl(record.sequence)
  for record in t_records:
    t_proteins.incl(record.sequence)
  let tp = intersection(r_proteins,t_proteins).len
  let fp = difference(t_proteins,r_proteins).len
  let fn = difference(r_proteins,t_proteins).len
  echo "TP: ", tp
  echo "FP: ", fp
  echo "FN: ", fn
  echo ""
  echo &"Precision: {float(tp) / float(tp+fp)}"
  echo &"Recall:    {float(tp) / float(tp+fn)}"

proc compareBLASTPTranslations*(reference_infilepath : string, blastp_infilepath : string,) =
  var fp = 0
  var reference_id_set : HashSet[string]
  var match_set : HashSet[string]
  var ref_infile : File
  discard open(ref_infile,reference_infilepath,fmRead)
  let reference_records = poGraphUtils.parseFasta(ref_infile)
  ref_infile.close()
  for record in reference_records:
    reference_id_set.incl(record.read_id)
  
  let blast_records = parseBLASTPoutput(blastp_infilepath)
  for record in blast_records:
    if record.match_names.len > 0:
      match_set.incl(record.match_names[0])
      # tp += 1
    else:
      fp += 1

  let tp = match_set.len
  let fn = difference(reference_id_set,match_set).len
  echo "TP: ", tp
  echo "FP: ", fp
  echo "FN: ", fn
  echo ""
  echo &"Precision: {float(tp) / float(tp+fp)}"
  echo &"Recall:    {float(tp) / float(tp+fn)}"

proc extractIntronsFromBED12*(infilepath : string, outfilepath : string) = 
  ## Grabs Introns from BED12 file and reports one per line in BED format, with ID = ID from the BED12 line
  var infile,outfile : File
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  try:
    while true:
      let bedline = infile.readLine()
      let bedfields = bedline.split(sep='\t')
      let chr = bedfields[0]
      let start_idx = parseUInt(bedfields[1])
      let end_idx = parseUInt(bedfields[2])
      let txid = bedfields[3]
      let strand = bedfields[5]
      let block_sizes = bedfields[10].split(sep=',')
      let block_starts= bedfields[11].split(sep=',')

      # var exon_count = 1
      # var reverse_exon_count = block_starts.len
      for i in 1..<block_starts.len:
        let intron_start_idx = start_idx + parseUInt(block_starts[i-1]) + parseUInt(block_sizes[i-1])
        let intron_end_idx = start_idx + parseUInt(block_starts[i])
        outfile.write(&"{chr}\t{intron_start_idx}\t{intron_end_idx}\t{txid}\t.\t{strand}\n")
  except EOFError:
    discard
  infile.close()
  outfile.close()

proc testOverlap(read : (uint64,uint64),single_exon_gene_list : seq[(uint64,uint64,string)]) : string = 
  result = ""
  for gene in single_exon_gene_list: #gene list must be sorted
    if read[0] >= gene[0] and read[0] <= gene[1]:
      if result != "":
        if result != gene[2]:
          break
      else:
          result = gene[2]
    elif read[1] >= gene[0] and read[1] <= gene[1]:
      if result != "":
        if result != gene[2]:
          break
      else:
          result = gene[2]
    elif read[0] <= gene[0] and read[1] >= gene[1]:
      if result != "":
        if result != gene[2]:
          break
      else:
          result = gene[2]
    elif read[1] <= gene[0]:
      break

proc callNonCanonicalSplicingFromFASTA*(infilepath : string,outfilepath : string) =
  var infile,outfile : File
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  let records = parseFasta(infile)
  infile.close()
  for record in records:
    if record.sequence.len < 2:
      echo "WARNING - Very short intron detected"
    else:
      if record.sequence[0..1].toUpperAscii() != "GT" or record.sequence[^2..^1].toUpperAscii() != "AG":
        outfile.writeLine(record.read_id)
  outfile.close()

proc splitFASTAByReadCounts*(infilepath : string, outfile_prefix : string, bins : openArray[uint64] = [1'u64,2'u64,5'u64,10'u64,20'u64,40'u64,80'u64,160'u64,320'u64,640'u64]) = 
  var infile : File
  discard open(infile,infilepath, fmRead)
  let fasta_records = parseFasta(infile)
  infile.close()
  var split_records : seq[seq[FastaRecord]]
  for i in 0..<bins.len:
    split_records.add(@[])
  for record in fasta_records:
    let num_reads = uint64(parseUInt(record.read_id.split('_')[^1]))
    var bin = -1
    for i in 0..<(bins.len - 1):
      if num_reads >= bins[i] and num_reads < bins[i+1]:
        bin = i
    if num_reads >= bins[^1]:
      bin = bins.len - 1
    if bin != -1:
      split_records[bin].add(record)
  for i,split in split_records:
    var bin_id = ""
    if i == bins.len - 1:
      bin_id = &"{bins[^1]}+"
    else:
      bin_id = &"{bins[i]}-{bins[i+1] - 1}"
    var outfile : File
    discard open(outfile,&"{outfile_prefix}_{bin_id}.fa",fmWrite)
    for record in split:
      outfile.write(&">{record.read_id}\n")
      outfile.writeLine(record.sequence)
    outfile.close()

proc filterFASTAByReadCounts*(infilepath,outfilepath : string, filter : uint64 = 5'u64) =
  var infile,outfile : File
  discard open(infile,infilepath, fmRead)
  let fasta_records = parseFasta(infile)
  infile.close()
  discard open(outfile,outfilepath,fmWrite)
  for record in fasta_records:
    let num_reads = uint64(parseUInt(record.read_id.split('_')[^1]))
    if num_reads >= filter:
      outfile.write(&">{record.read_id}\n")
      outfile.writeLine(record.sequence)
  outfile.close()

proc parseAttributes(s : string) : Table[string,string] = 
  let fields = s.split(';')[0..^2]
  for field in fields:
    let fields1 = field.strip(leading=true,trailing=true,chars={' '}).split(' ')
    let key = fields1[0]
    let val = fields1[1].strip(leading=true,trailing=true,chars = {'"'})
    result[key] = val

proc callNovelNonCanonical(reference_infilepath, infilepath,outfilepath : string,threshold : uint = 5) =
  var reference_infile,infile,outfile : File
  discard open(reference_infile,reference_infilepath,fmRead)
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  var reference_introns : HashSet[(string,uint64,uint64)]
  var last_exons : Table[string,(string,uint64,uint64)]
  try:
    while true:
      let line = reference_infile.readLine()
      if line.len == 0:
        continue
      elif line[0] == '#':
        continue
      let fields0 = line.split('\t')
      # echo fields0
      # if fields0.len == 0:
      #   echo line
      if fields0[2] != "exon":
        continue
      let attributes = parseAttributes(fields0[8])
      if "transcript_id" in attributes:
        let new_chr = fields0[0]
        let new_start_idx =  uint64(parseUInt(fields0[3])) - 1'u32
        let new_end_idx = uint64(parseUInt(fields0[4]))
        if attributes["transcript_id"] in last_exons:
          let (chr, start_idx,end_idx)= last_exons[attributes["transcript_id"]]
          try:
            assert chr == new_chr
            reference_introns.incl((chr, end_idx, new_start_idx))
          except AssertionError:
            echo "ERROR - same transcript on different chromosomes:"
            let tx_id = attributes["transcript_id"]
            echo &"{chr}:{start_idx}-{end_idx} {tx_id}"
            echo &"{new_chr}:{new_start_idx}-{new_end_idx} {tx_id}"
          # outfilepath.writeLine(chr,"\t",end_idx,"\t",new_start_idx,"\t","reference_intron")
        last_exons[attributes["transcript_id"]] = (new_chr, new_start_idx, new_end_idx)
      else:
        echo "ERROR - no field transcript_id"
  except EOFError:
    discard
  reference_infile.close()
  var novel_counter = 0
  var total_counter = 0
  try:
    while true:
      let line = infile.readLine()
      if line.len == 0:
        continue
      let fields0 = line.split(':')
      let cluster_id = fields0[0]
      let read_support = parseUInt(cluster_id.split('_')[^1])
      let chr = fields0[2]
      let indices = fields0[3].strip(chars={'(',')','+','-'})
      let split_indices = indices.split('-')
      let start_idx = uint64(parseUInt(split_indices[0])) 
      let end_idx = uint64(parseUInt(split_indices[1]))
      if read_support >= threshold:
        total_counter += 1
        if (chr,start_idx,end_idx) notin  reference_introns:
          outfile.writeLine(&"{chr}\t{start_idx}\t{end_idx}\tquery_intron")
          novel_counter += 1
  except EOFError:
    discard
  infile.close()
  outfile.close()
  echo &"Total introns above threshold - {total_counter}"
  echo &"Novel introns - {novel_counter}"

proc callOverlappingNonCanonical(reference_infilepath, infilepath,outfilepath : string) =
  var reference_infile,infile,outfile : File
  discard open(reference_infile,reference_infilepath,fmRead)
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  var reference_introns : HashSet[(string,uint64,uint64)]
  # var last_exons : Table[string,(string,uint64,uint64)]
  try:
    while true:
      let line = reference_infile.readLine()
      if line.len == 0:
        continue
      let fields0 = line.split(':')
      let chr = fields0[2]
      let indices = fields0[3].strip(chars={'(',')','+','-'})
      let split_indices = indices.split('-')
      let start_idx = uint64(parseUInt(split_indices[0])) 
      let end_idx = uint64(parseUInt(split_indices[1]))
      reference_introns.incl((chr,start_idx,end_idx))
  except EOFError:
    discard
  reference_infile.close()
  var novel_counter = 0
  var overlapping_counter = 0
  try:
    while true:
      let line = infile.readLine()
      if line.len == 0:
        continue
      let fields0 = line.split(':')
      let chr = fields0[2]
      let indices = fields0[3].strip(chars={'(',')','+','-'})
      let split_indices = indices.split('-')
      let start_idx = uint64(parseUInt(split_indices[0])) 
      let end_idx = uint64(parseUInt(split_indices[1]))
      if (chr,start_idx,end_idx) notin  reference_introns:
        novel_counter += 1
      else:
        outfile.writeLine(line)
        overlapping_counter += 1
  except EOFError:
    discard
  infile.close()
  outfile.close()
  echo &"Non-overlapping introns - {novel_counter}"
  echo &"Overlapping introns - {overlapping_counter}"

proc assignTxIDs(reference_infilepath,infilepath,outfilepath : string) = 
  var reference_infile,infile,outfile : File
  discard open(reference_infile,reference_infilepath,fmRead)
  var tx_exons : Table[(string,string),seq[(uint64,uint64)]]
  try:
    while true:
      let line = reference_infile.readLine()
      if line.len == 0:
        continue
      elif line[0] == '#':
        continue
      let fields0 = line.split('\t')
      # echo fields0
      # if fields0.len == 0:
      #   echo line
      if fields0[2] != "exon":
        continue
      let attributes = parseAttributes(fields0[8])
      if "transcript_id" in attributes:
        let new_chr = fields0[0]
        let new_start_idx =  uint64(parseUInt(fields0[3])) - 1'u32
        let new_end_idx = uint64(parseUInt(fields0[4]))
        if (attributes["transcript_id"],new_chr) in tx_exons:
          tx_exons[(attributes["transcript_id"],new_chr)].add((new_start_idx,new_end_idx))
        else:
          tx_exons[(attributes["transcript_id"],new_chr)] = @[(new_start_idx,new_end_idx)]
      else:
        echo "ERROR - no field transcript_id"
  except EOFError:
    discard
  reference_infile.close()
  var tx_introns : Table[(string,seq[(uint64,uint64)]), string]
  var single_exon_genes : Table[string,seq[(uint64,uint64,string)]]
  for (tx_id,chr) in tx_exons.keys:
    # echo chr
    let exon_chain = tx_exons[(tx_id,chr)]
    if exon_chain.len > 1:
      var intron_chain : seq[(uint64,uint64)]
      for i in 1..<exon_chain.len:
        intron_chain.add((exon_chain[i-1][1],exon_chain[i][0]))
      tx_introns[(chr,intron_chain)] = tx_id
    else:
      if chr in single_exon_genes:
        single_exon_genes[chr].add((exon_chain[0][0],exon_chain[0][1],tx_id))
      else:
        single_exon_genes[chr] = @[(exon_chain[0][0],exon_chain[0][1],tx_id)]
  for chr in single_exon_genes.keys:
    single_exon_genes[chr].sort
  
  #Read in BED12 file and get intron chains
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  try:
    while true:
      let bedline = infile.readLine()
      let bedfields = bedline.split(sep='\t')
      let chr = bedfields[0]
      let start_idx = parseUInt(bedfields[1])
      let end_idx = parseUInt(bedfields[2])
      let txid = bedfields[3]
      let strand = bedfields[5]
      let block_sizes = bedfields[10].split(sep=',')
      let block_starts= bedfields[11].split(sep=',')
      var reference_id = ""
      if block_starts.len > 1:
        var intron_chain : seq[(uint64,uint64)]
        for i in 1..<block_starts.len:
          let intron_start_idx = uint64(start_idx + parseUInt(block_starts[i-1]) + parseUInt(block_sizes[i-1]))
          let intron_end_idx = uint64(start_idx + parseUInt(block_starts[i]))
          intron_chain.add((intron_start_idx,intron_end_idx))
        if (chr,intron_chain) in tx_introns:
          reference_id = tx_introns[(chr,intron_chain)]
      else:
        reference_id = testOverlap((uint64(start_idx),uint64(end_idx)),single_exon_genes[chr])
      if reference_id == "":
        outfile.write(&"{txid}\t.\n")
      else:
        outfile.write(&"{txid}\t{reference_id}\n")
  except EOFError:
    discard
  infile.close()
  outfile.close()

proc parseOptions() : UtilOptions = 
  
  var i = 0
  var mode = ""
  var last = ""

  var help_flag = true
  var run_flag = true

  var fastq = false
  var fastq_flag = false

  var min_length = 75'u64
  var min_length_flag = false

  var infilepath = ""
  var infilepath_flag = false

  var stranded = false

  var reference_infilepath = ""
  var reference_infilepath_flag = false
  
  var outfilepath = ""
  var outfilepath_flag = false

  for kind, key, val in getopt():
    # echo kind," ", key," ", val
    if i == 0:
      mode = key
      i += 1
      continue
    if i == 1:
      help_flag = false
    i += 1
    case mode:
      of "translate", "strandTranscripts", "bed2gtf", "parseBLASTP","compareBLASTP","compareFASTA","splitFASTA","filterFASTA","extractIntrons","callNonCanonical","callNovelNonCanonical","callOverlapping","assignIDs":
        case kind:
          of cmdEnd:
            break
          of cmdShortOption, cmdLongOption:
            if last != "":
              echo &"ERROR - Option {last} provided without an argument"
              run_flag = false
              break
            case key:
              of "l", "min-length":
                if not min_length_flag:
                  min_length_flag = true
                  if val != "":
                    min_length = parseUInt(val)
                  else:
                    last = "min-length"
                else:
                  echo "ERROR - Multiple min-length provided"
                  run_flag = false
              of "i":
                if not infilepath_flag:
                  infilepath_flag = true
                  if val != "":
                    infilepath = val
                  else:
                    last = "infile"
                else:
                  echo "ERROR - Multiple infiles provided"
                  run_flag = false
                  break
              of "q", "fastq":
                if not fastq_flag:
                  fastq = true
                  fastq_flag = true
                elif not fastq:
                  echo "ERROR - Conflicting flags -q / --fastq and -a / --fasta"
                  run_flag = false
                  break
              of "a", "fasta":
                if not fastq_flag:
                  fastq_flag = true
                elif fastq:
                  echo "ERROR - Conflicting flags -q / --fastq and -a / --fasta"
                  run_flag = false
                  break
              of "o":
                if not outfilepath_flag:
                  outfilepath_flag = true
                  if val != "":
                    infilepath = val
                  else:
                    last = "outfile"
                else:
                  echo "ERROR - Multiple outfiles provided"
                  run_flag = false
                  break
              of "r":
                if not reference_infilepath_flag:
                  reference_infilepath_flag = true
                  if val != "":
                    reference_infilepath = val
                  else:
                    last = "reference"
                else:
                  echo "ERROR Multiple references provided"
                  run_flag = false
                  break
              of "s", "stranded":
                stranded = true
              of "h", "help":
                help_flag = true
                run_flag = false
                break
          of cmdArgument:
            case last:
              of "min-length":
                min_length = parseUInt(key)
              of "infile":
                infilepath = key
              of "outfile":
                outfilepath = key
              of "reference":
                reference_infilepath = key
              else:
                echo &"ERROR - unknown option {last} provided"
                run_flag = false
                break
            last = ""
      else:
        help_flag = true
        break
  if infilepath == "" and not help_flag:
    echo "ERROR - infilepath must be specified"
    run_flag = false
    help_flag = true
  var require_outfile,require_reference = false
  case mode:
    of "translate","bed2gtf","parseBLASTP","splitFASTA","filterFASTA","extractIntrons":
      require_outfile = true
    of "compareBLASTP","compareFASTA":
      require_reference = true
    of "callNonCanonical","callNovelNonCanonical","callOverlapping":
      require_reference = true
      require_outfile = true
  
  if require_outfile and outfilepath == "" and not help_flag:
    echo "ERROR - outfilepath must be specified"
    run_flag = false
    help_flag = true
  if require_reference and reference_infilepath == "" and not help_flag:
    echo "ERROR - referencefilepath must be specified"
    run_flag = false
    help_flag = true
  if help_flag:
    case mode:
      of "translate":
        writeTranslateHelp()
      of "strandTranscripts":
        writeStrandTranscriptsHelp()
      of "bed2gtf":
        writeBED2GTFHelp()
      of "parseBLASTP":
        writeBLASTPHelp()
      of "compareBLASTP":
        writeCompareBLASTPHelp()
      of "compareFASTA":
        writeCompareFASTAHelp()
      of "splitFASTA":
        writeSplitFASTAHelp()
      of "filterFASTA":
        writeFilterFASTAHelp()
      of "extractIntrons":
        writeExtractIntronsHelp()
      of "callNonCanonical":
        writeCallNonCanonicalHelp()
      of "callNovelNonCanonical":
        writeCallNovelNonCanonicalHelp()
      of "callOverlapping":
        writeCallOverlappingHelp()
      else:
        echo "ERROR - first argument must specify utility function"
        writeDefaultHelp()
  return UtilOptions(mode : mode,
                     run_flag : run_flag,
                     infilepath : infilepath,
                     outfilepath : outfilepath,
                     reference_infilepath : reference_infilepath,
                     min_length : min_length,
                     fastq : fastq,
                     stranded : stranded
                     )

proc main() =
  let opt = parseOptions()
  if opt.run_flag:
    case opt.mode:
      of "translate":
        var infile : File
        discard open(infile,opt.infilepath,fmRead)
        if opt.fastq:
          let records = parseFASTQ(infile)
          infile.close()
          translateTranscripts(records,opt.outfilepath,threshold = int(opt.min_length),stranded = opt.stranded)
        else:
          let records = poGraphUtils.parseFasta(infile)
          infile.close()
          translateTranscripts(records,opt.outfilepath,threshold = int(opt.min_length),stranded = opt.stranded)
      of "strandTranscripts":
        var infile : File
        discard open(infile,opt.infilepath,fmRead)
        if opt.fastq:
          let records = parseFASTQ(infile)
          infile.close()
          strandTranscripts(records,opt.outfilepath)
        else:
          let records = poGraphUtils.parseFasta(infile)
          infile.close()
          strandTranscripts(records,opt.outfilepath)
      of "bed2gtf":
        convertBED12toGTF(opt.infilepath,opt.outfilepath,opt.stranded)
      of "parseBLASTP":
        let blast_output = parseBLASTPoutput(opt.infilepath)
        var outfile : File
        discard open(outfile,opt.outfilepath,fmWrite)
        for blast_match in blast_output:
          if blast_match.identities_percentages.len != 0:
            outfile.write(&"{blast_match.query_name}\t{blast_match.match_names[0]}\t{blast_match.evals[0][0]}\n")
          else:
            outfile.write(&"{blast_match.query_name}\t----\t----\n")
        outfile.close()
      of "compareBLASTP":
        compareBLASTPTranslations(opt.reference_infilepath,opt.infilepath)
      of "compareFASTA":
        compareExactTranslations(opt.reference_infilepath,opt.infilepath)
      of "splitFASTA":
        splitFASTAByReadCounts(opt.infilepath,opt.outfilepath)
      of "filterFASTA":
        filterFASTAByReadCounts(opt.infilepath,opt.outfilepath)
      of "extractIntrons":
        extractIntronsFromBED12(opt.infilepath,opt.outfilepath)
      of "callNonCanonical":
        callNonCanonicalSplicingFromFASTA(opt.infilepath,opt.outfilepath)
      of "callNovelNonCanonical":
        callNovelNonCanonical(opt.reference_infilepath,opt.infilepath,opt.outfilepath)
      of "callOverlapping":
        callOverlappingNonCanonical(opt.reference_infilepath,opt.infilepath,opt.outfilepath)
      of "assignIDs":
        assignTxIDs(opt.reference_infilepath,opt.infilepath,opt.outfilepath)



main()