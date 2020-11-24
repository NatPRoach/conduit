# import os
# import osproc
import parseopt
import strutils
import strformat
import tables
import sets
import algorithm
import fasta
import fastq
import na
import version

type
  BLASTmatch* = object
    queryName* : string
    queryLen* : uint
    matchNames* : seq[string]
    matchLens* : seq[uint]
    scores* : seq[seq[float]]
    evals* : seq[seq[float]]
    identitiesNumerators* : seq[seq[uint]]
    identitiesDenominators* : seq[seq[uint]]
    identitiesPercentages* : seq[seq[float]]
    positivesNumerators* : seq[seq[uint]]
    positivesDenominators* : seq[seq[uint]]
    positivesPercentages* : seq[seq[float]]
    gapNumerators* : seq[seq[uint]]
    gapDenominators* : seq[seq[uint]]
    gapPercentages* : seq[seq[float]]
    querySeqs* : seq[seq[string]]
    conssSeqs* : seq[seq[string]]
    sbjctSeqs* : seq[seq[string]]
  
  UtilOptions = object
    mode : string
    runFlag : bool
    infilepath : string
    outfilepath : string
    referenceInfilepath : string
    referenceInfilepath2 : string
    minLength : uint64
    fastq : bool
    stranded : bool
  
  GTFTranscript* = object
    chr* : string
    strand* : char
    startIdx* : uint64
    endIdx* : uint64
    introns* : seq[(uint64,uint64)]


# proc intersect(a,b : HashSet[GTFTranscript]) =
#   var a_,b_ = HashSet[(string,char,seq[(uint64,uint64)])]
#   for tx in a.items:
#     a_.incl((tx.chr,tx.strand,tx.introns))
#   for tx in b.items:
#     b_.incl(())


proc conduitUtilsVersion() : string =
  return &"CONDUIT Utilities Version {version.ConduitVersion} by Nathan Roach" &
       "\n( nroach2@jhu.edu, https://github.com/NatPRoach/conduit/ )"


proc writeDefaultHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "Usage:"
  echo "  ./conduitUtils <function>"
  echo "Where <function> is one of the following:"
  echo "  translate     - Translates FASTA/Q nucleotide sequences into protein"
  echo "                  based on their longest ORF"
  echo ""
  echo "  bed2gtf       - Converts BED12 files to well structured GTF file"
  echo "                  suitable for use in GFFcompare"
  echo ""
  echo "  parseBLASTP   - Parses BLASTP output and outputs closest match for"
  echo "                  each query transcript as determined by BLASTP"
  echo ""
  echo "  compareBLASTP - Compares BLASTP output and reference proteome to "
  echo "                  determine the # of TP, FP, and FN for a sample"
  echo ""
  echo "  compareFASTA  - Compares two FASTA files, an input and a reference, "
  echo "                  to determine the # of TP, FP, and FN for a sample"
  echo ""
  echo "  splitFASTA    - Splits CONDUIT produced FASTA file based on the "
  echo "                  number of reads supporting each isoform"
  echo ""
  echo "  filterFASTA   - Filters CONDUIT produced FASTA file based on number "
  echo "                  of reads supporting each isoform"
  echo ""
  echo "  extractIntrons   - Extracts out intronic sequences from BED12"
  echo "                     formatted input and outputs as BED6"
  echo ""
  echo "  callOverlapping  - Compares two files of readIDs specifying introns "
  echo "                     in the format produced by:"
  echo "                       `bedtools getfasta -name`"
  echo "                     and reports the introns that are shared between"
  echo "                     the two files (not stranded)"
  echo ""
  echo "  callNonCanonical - Reads in a FASTA file and reports the readIDs of "
  echo "                     sequences that dont begin with GT and end with AG"
  echo ""
  echo "  callNovelNonCanonical - Compares introns described by reference GTF "
  echo "                          file to introns described by a list of "
  echo "                          readIDs in the format produced by:"
  echo "                            `bedtools getfasta -name`"
  echo "                          function, outputs the novel introns in BED"
  echo "                          format"
  echo ""


proc writeTranslateHelp() =
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "translate   - Translates FASTA/Q nucleotide sequences into protein "
  echo "              based on their longest ORF"
  echo "Usage:"
  echo "  ./conduitUtils translate [options] -i <transcripts.fa>"
  echo "                                     -o <predicted_protein.fa>"
  echo "  <transcripts.fa>         FASTA/Q infile containing putative"
  echo "                           transcripts to be translated"
  echo "  <predicted_protein.fa>   FASTA outfile containing in silico"
  echo "                           translated ORFs from transcripts.fa"
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
  echo "        Minimum length in Amino Acids necessary for a putative ORF to "
  echo "        be reported"


proc writeStrandTranscriptsHelp() =
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "strandTranscripts - Translates FASTA/Q nucleotide sequences into"
  echo "                    protein based on their longest ORF"
  echo "Usage:"
  echo "  ./conduitUtils strandTranscripts [options] -i <transcripts.fa> \\"
  echo "                                             -o <predicted_protein.fa>"
  echo "  <transcripts.fa>            FASTA/Q infile containing putative "
  echo "                              transcripts to be translated"
  echo "  <stranded_trasncripts.fa>   FASTA outfile containing in silico "
  echo "                              translated ORFs from transcripts.fa"
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
  echo "bed2gtf - Converts BED12 files to well structured GTF file suitable for"
  echo "          use in GFFcompare"
  echo "Usage:"
  echo "  ./conduitUtils bed2gtf -i <infile.bed> -o <outfile.gtf>"
  echo "  <infile.bed>    BED12 infile to be converted in to GTF format"
  echo "  <outfile.gtf>   GTF outfile"
  echo ""
  echo "Options (defaults in parentheses):"
  echo "  Input Options:"
  echo "    -s, --stranded"
  echo "        Report gtf fields with strand information"


proc writeBLASTPHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "parseBLASTP - Parses BLASTP output and outputs closest match for each"
  echo "              query transcript as determined by BLASTP"
  echo "Usage:"
  echo "  ./conduitUtils parseBLASTP -i <inBLASTP.txt> \\"
  echo "                             -o <outPutativeOrthologs.tsv>"
  echo "  <inBLASTP.txt>              Default output of BLASTP search of"
  echo "                              translated protein products vs some"
  echo "                              reference proteome"
  echo "  <outPutativeOrthologs.tsv>  Tab separated file of putative ortholog"
  echo "                              matches"
  echo "    Output will be in format:"
  echo "      <Query ID>\\t<Reference proteome top match ID>\\t<E value>"


proc writeCompareBLASTPHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "compareBLASTP - Compares BLASTP output and reference proteome to"
  echo "                determine the # of true positives, false positives, and"
  echo "                false negatives for a sample"
  echo "Usage:"
  echo "  ./conduitUtils compareBLASTP -r <reference_proteome.fa> \\"
  echo "                               -i <inBLASTP.txt>"
  echo "  <reference_proteome.fa>     FASTA file describing the reference"
  echo "                              proteome used in the BLASTP search"
  echo "  <inBLASTP.txt>              Default output of BLASTP search of"
  echo "                              translated protein products vs some"
  echo "                              reference proteome"


proc writeCompareFASTAHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "compareFASTA  - Compares two FASTA files, an input and a reference, to"
  echo "                determine the # of true positives, false positives, and"
  echo "                false negatives for a sample"
  echo "Usage:"
  echo "  ./conduitUtils compareFASTA -r <reference.fa> -i <query.fa>"
  echo "  <reference.fa>  Reference FASTA file defining the truth set"
  echo "  <query.fa>      Query FASTA files defining the query set"


proc writeSplitFASTAHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "splitFASTA - Splits CONDUIT produced FASTA file based on the number of"
  echo "             reads supporting each isoform"
  echo "Usage:"
  echo "  ./conduitUtils splitFASTA -i <conduit_output.fa> -o <outprefix>"
  echo "  <conduit_output.fa> CONDUIT produced FASTA file to be split based on"
  echo "                      number of reads supporting each isoform"
  echo "  <outprefix>         Prefix for the fasta files to be output, suffix "
  echo "                      will describe the bin being reported"


proc writeFilterFASTAHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "filterFASTA - Filters CONDUIT produced FASTA file based on number of"
  echo "              reads supporting each isoform"
  echo "Usage:"
  echo "  ./conduitUtils filterFASTA -i <inBLASTP.txt> \\"
  echo "                             -o <outPutativeOrthologs.tsv>"
  echo "  <conduit_output.fa> CONDUIT produced FASTA file to be filtered based"
  echo "                      on number of reads supporting each isoform"
  echo "  <filtered.fa>       Output FASTA file for filtered reads"
  echo "Options: (defaults in parentheses)"
  echo "  Filtering options:"
  echo "     -n (5)"
  echo "        Minimum number of reads that must support an isoform for it to"
  echo "        be reported in the filtered FASTA"


proc writeExtractIntronsHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "extractIntrons - Extracts out intronic sequences from BED12 formatted "
  echo "                 input and outputs as BED6"
  echo "Usage:"
  echo "  ./conduitUtils extractIntrons -i <transcripts.bed12> -o <introns.bed>"
  echo "  <transcripts.bed12>         Transcripts in BED12 format to extract "
  echo "                              introns from"
  echo "  <introns.bed>               BED6 output of extracted introns"


proc writeCallNonCanonicalHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "callNonCanonical - Reads in a FASTA file and reports the readIDs of "
  echo "                   sequences that dont begin with GT and end with AG"
  echo "Usage:"
  echo "  ./conduitUtils callNonCanonical -i <introns.fa> -o <noncanonical.txt>"
  echo "  <introns.fa>       FASTA describing the stranded sequence of introns"
  echo "                     extracted from `extractIntrons` Introns sequences"
  echo "                     can be obtained using `bedtools getfasta -name -s`"
  echo "  <noncanonical.txt> Read IDs of the sequences that didn't begin with"
  echo "                     GT and end with AG"


proc writeCallNovelNonCanonicalHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "callNovelNonCanonical - Compares introns described by reference GTF"
  echo "                        file to introns described by a list of readIDs"
  echo "                        in the format produced by"
  echo "                          `bedtools getfasta -name`"
  echo "                        function, outputs the novel introns in BED"
  echo "                        format"
  echo "Usage:"
  echo "  ./conduitUtils callNovelNonCanonical -r <reference.gtf> \\"
  echo "                                       -i <noncanonical.txt> \\"
  echo "                                       -o <novel.bed>"
  echo "  <reference.gtf>              Reference GTF file specifying the"
  echo "                               introns to compare against"
  echo "  <noncanonical.txt>           Read IDs specifying intron structure in"
  echo "                               the format produced by"
  echo "                               `bedtools getfasta -name`"
  echo "  <novel.bed>                  Output of introns found in the"
  echo "                               noncanonical.txt file but not found in"
  echo "                               the reference, in BED6 format"


proc writeCallOverlappingHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "callOverlapping - Compares two files of readIDs specifying introns in"
  echo "                  the format produced by `bedtools getfasta -name`, and"
  echo "                  reports the introns that are shared between the two"
  echo "                  files (not stranded)"
  echo "Usage:"
  echo "  ./conduitUtils callOverlapping -r <introns1.txt> \\"
  echo "                                 -i <introns2.txt> \\"
  echo "                                 -o <shared_introns.txt>"
  echo "  <introns1.txt>              Read IDs specifying introns in the format"
  echo "                              produced by `bedtools getfasta -name`"
  echo "  <introns2.txt>              Read IDs specifying introns in the format"
  echo "                              produced by `bedtools getfasta -name`"
  echo "  <shared_introns.txt>        The introns in common between the two"
  echo "                              files"


proc parseBLASTPoutput*(infilepath : string) : seq[BLASTmatch] = 
  ##
  ## Takes in BLASTP default output file and outputs seq of matching proteins
  ##
  
  var infile : File
  discard open(infile,infilepath,fmRead)
  # Discard the header / citation info
  var queryLine : string
  var firstIter1 = true
  while true:
    queryLine = infile.readLine()
    if queryLine.len >= 7:
      if queryLine[0..6] == "Query= ":
        break
  while true:
    try:
      if not firstIter1:
        # Query= cluster_10_0 <unknown description>
        queryLine = infile.readLine()
      else:
        firstIter1 = false
      # echo queryLine
      if queryLine[0..11] == "  Database: ":
        break
      assert queryLine[0..6] == "Query= "
      var queryName = queryLine[7..^1]
      while true:
        let next = infile.readLine() # ""
        if next == "":
          break
        else:
          queryName = queryName & next

      var matchNames : seq[string]
      var matchLens  : seq[uint]
      var evals : seq[seq[float]]
      var scores : seq[seq[float]]
      var identitiesNumerators : seq[seq[uint]]
      var identitiesDenominators : seq[seq[uint]]
      var identitiesPercentages : seq[seq[float]]
      var positivesNumerators : seq[seq[uint]]
      var positivesDenominators : seq[seq[uint]]
      var positivesPercentages : seq[seq[float]]
      var gapNumerators : seq[seq[uint]]
      var gapDenominators : seq[seq[uint]]
      var gapPercentages : seq[seq[float]]
      var querySeqs : seq[seq[string]]
      var conssSeqs : seq[seq[string]]
      var sbjctSeqs : seq[seq[string]]
      

      let queryLengthLine = infile.readLine() # "Length=229"
      # echo queryLengthLine
      assert queryLengthLine[0..6] == "Length="
      let queryLength = parseUInt(queryLengthLine[7..^1])
      # echo queryLength
      let splitline = infile.readLine()
      if splitline == "":
        # No significant alignments here
        discard infile.readLine()  # ""
        assert infile.readLine() == "***** No hits found *****"
        discard infile.readLine()  # ""
        discard infile.readLine()  # ""
        discard infile.readLine()  # ""
      elif splitline == "                                   " &
                        "                                   Score        E":
        assert infile.readLine() == "Sequences producing significant " &
          "alignments:                          (Bits)     Value"
        discard infile.readLine()  # ""
        var nextline : string
        var numSeqs = 0
        while true:
          nextline =  infile.readLine()
          if nextline.len == 0:
            continue
          if nextline[0] == '>':
            break
          numSeqs += 1
        # echo numSeqs
        for i in 0..<numSeqs:
          # echo "seq num: ", i
          # echo nextline
          var matchName = nextline[2..^1]
          var firstIter = true
          while true:
            if nextline.len == 0:
              nextline = infile.readLine()
              continue
            # echo "nextline: ", nextline
            if nextline.len >= 7:
              if nextline[0..6] == "Length=":
                break
            if firstIter:
              firstIter = false
            else:
              matchName = matchName & nextline
            nextline = infile.readLine()
          matchNames.add(matchName)
          # echo "match name - ", matchName
          let matchLength = parseUInt(nextline[7..^1])
          # echo "match len - ", matchLength
          matchLens.add(matchLength)
          var scoresForMatch : seq[float]
          var evalsForMatch : seq[float]
          var idsNumeratorsForMatch : seq[uint]
          var idsDenominatorsForMatch : seq[uint]
          var idsPercentagesForMatch : seq[float]
          var posNumeratorsForMatch : seq[uint]
          var posDenominatorsForMatch : seq[uint]
          var posPercentagesForMatch : seq[float]
          var gapNumeratorsForMatch : seq[uint]
          var gapDenominatorsForMatch : seq[uint]
          var gapPercentagesForMatch : seq[float]
          var querySeqsForMatch : seq[string]
          var conssSeqsForMatch : seq[string]
          var sbjctSeqsForMatch : seq[string]
          discard infile.readLine()
          while true:
            nextline = infile.readLine().strip()
            if nextline.len == 0:
              break
            elif nextline[0..7] == "Score = ":
              var fields0 = nextline.split(seps={','})
              var fields1 = fields0[0].splitWhitespace()
              let score = parseFloat(fields1[2])
              scoresForMatch.add(score)
              fields1 = fields0[1].strip().splitWhitespace()
              let evalue = parseFloat(fields1[2])
              evalsForMatch.add(evalue)
              fields1 = fields0[2].strip().split(seps={':'})
              # let blastMethod = fields1[1].strip()
              nextline = infile.readLine().strip()
              fields0 = nextline.split(seps={','})
              fields1 = fields0[0].strip().splitWhitespace()
              var fields2 = fields1[2].strip().split(seps={'/'})
              let idsNumerator = parseUInt(fields2[0])
              let idsDenominator = parseUInt(fields2[1])
              let idsPercentage = 100.0 * float(idsNumerator) /
                float(idsDenominator)
              
              idsNumeratorsForMatch.add(idsNumerator)
              idsDenominatorsForMatch.add(idsDenominator)
              idsPercentagesForMatch.add(idsPercentage)

              fields1 = fields0[1].strip().splitWhitespace()
              fields2 = fields1[2].strip().split(seps={'/'})
              let posNumerator = parseUInt(fields2[0])
              let posDenominator = parseUInt(fields2[1])
              let posPercentage = 100.0 * float(posNumerator) /
                float(posDenominator)

              posNumeratorsForMatch.add(posNumerator)
              posDenominatorsForMatch.add(posDenominator)
              posPercentagesForMatch.add(posPercentage)

              fields1 = fields0[2].strip().splitWhitespace()
              fields2 = fields1[2].strip().split(seps={'/'})
              let gapNumerator = parseUInt(fields2[0])
              let gapDenominator = parseUInt(fields2[1])
              let gapPercentage = 100.0 * float(gapNumerator) /
                float(gapDenominator)

              gapNumeratorsForMatch.add(gapNumerator)
              gapDenominatorsForMatch.add(gapDenominator)
              gapPercentagesForMatch.add(gapPercentage)

              discard infile.readLine() # ""
              var querySeq : string
              var sbjctSeq : string
              var conssSeq : string
              var queryStartIdx : uint
              var sbjctStartIdx : uint
              var firstIter = true
              while true:
                nextline = infile.readLine()
                if nextline == "":
                  break
                let queryline = nextline
                let consensusline = infile.readLine()
                let sbjctline = infile.readLine()
                discard infile.readLine()
                let queryfields = queryline.strip().splitWhitespace()
                let conssSubseq = consensusline[12..^1] # Have to do it this way
                #  because sometimes there'll be a space at the beginning
                let sbjctfields = sbjctline.strip().splitWhitespace()
                # echo queryfields
                if firstIter:
                  firstIter = false
                  queryStartIdx = parseUInt(queryfields[1])
                  sbjctStartIdx = parseUInt(sbjctfields[1])
                let querySubseq = queryfields[2]
                let sbjctSubseq = sbjctfields[2]
                querySeq = querySeq & querySubseq
                conssSeq = conssSeq & conssSubseq
                sbjctSeq = sbjctSeq & sbjctSubseq
              # echo querySeq
              # echo conssSeq
              # echo sbjctSeq
              querySeqsForMatch.add(querySeq)
              conssSeqsForMatch.add(conssSeq)
              sbjctSeqsForMatch.add(sbjctSeq)
              # nextline = infile.readLine()
            else:
              break
          
          evals.add(evalsForMatch)
          scores.add(scoresForMatch)

          identitiesNumerators.add(idsNumeratorsForMatch)
          identitiesDenominators.add(idsDenominatorsForMatch)
          identitiesPercentages.add(idsPercentagesForMatch)

          positivesNumerators.add(posNumeratorsForMatch)
          positivesDenominators.add(posDenominatorsForMatch)
          positivesPercentages.add(posPercentagesForMatch)

          gapNumerators.add(gapNumeratorsForMatch)
          gapDenominators.add(gapDenominatorsForMatch)
          gapPercentages.add(gapPercentagesForMatch)

          querySeqs.add(querySeqsForMatch)
          conssSeqs.add(conssSeqsForMatch)
          sbjctSeqs.add(sbjctSeqsForMatch)
        # assert nextline == ""
      else:
        echo "ERROR PARSING BLASTP OUTPUT"
        break
      # discard infile.readLine() # ""
      # let test = infile.readLine()
      # echo test
      # assert test == "Lambda      K        H        a         alpha"
      assert infile.readLine() == "Lambda      K        H        " &
        "a         alpha"
      discard infile.readLine() #    0.308    0.126    0.365    0.792     4.96
      discard infile.readLine() # ""
      assert infile.readLine() == "Gapped"
      assert infile.readLine() == "Lambda      K        H        " &
        "a         alpha    sigma"
      #    0.267   0.0410    0.140     1.90     42.6     43.6
      discard infile.readLine() 
      discard infile.readLine() # ""
      discard infile.readLine() # "Effective search space used: 318188442"
      discard infile.readLine() # ""
      discard infile.readLine() # ""
      result.add(BLASTmatch( queryName : queryName,
                             queryLen : queryLength,
                             matchNames : matchNames,
                             matchLens : matchLens,
                             evals : evals,
                             scores : scores,
                             querySeqs : querySeqs,
                             conssSeqs : conssSeqs,
                             sbjctSeqs : sbjctSeqs,
                             identitiesNumerators : identitiesNumerators,
                             identitiesDenominators : identitiesDenominators,
                             identitiesPercentages : identitiesPercentages,
                             positivesNumerators : positivesNumerators,
                             positivesDenominators : positivesDenominators,
                             positivesPercentages : positivesPercentages,
                             gapNumerators : gapNumerators,
                             gapDenominators : gapDenominators,
                             gapPercentages : gapPercentages ))
    except EOFError:
      break
  infile.close()


proc translateORF*(nts : string,to_stop = true) : string =
  let translationTable = { "TTT" : 'F',
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
    let aa = translationTable[nts[i*3..((i+1)*3) - 1]]
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


proc translateTranscript*(nts : string) : string =
  let startCodonIndices = findAll(nts,"ATG")
  result = ""
  for start_codon_index in startCodonIndices:
    let translation = translateORF(nts[start_codon_index..^1])
    if translation.len > result.len:
      result = translation


proc translateTranscripts*(transcripts : openArray[FastaRecord],
                           outfilepath : string ,
                           threshold : int = 75,
                           wrap_len : int = 60,
                           stranded : bool = false) =
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
      outfile.write(&">{transcript.readId}\n")
      for i in 0..<(translation.len div wrap_len):
        outfile.write(&"{translation[i*wrap_len..(i+1)*wrap_len - 1]}\n")
      if translation.len mod wrap_len != 0 :
        let line = translation[wrap_len * (translation.len div wrap_len)..^1]
        outfile.write(&"{line}\n")
  outfile.close()


proc translateTranscripts*(transcripts : openArray[FastqRecord],
                           outfilepath : string ,
                           threshold : int = 75,
                           wrap_len : int = 60,
                           stranded = false) =
  translateTranscripts(convertFASTQtoFASTA(transcripts),
                       outfilepath,
                       threshold,
                       wrap_len,
                       stranded)


proc convertBED12toGTF*(infilepath : string,
                        outfilepath : string,
                        stranded : bool = false ) =
  ## Converts BED12 formatted file to well-formed GTF file
  ## suitable for evaluation with GFFcompare
  var infile,outfile : File
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  try:
    while true:
      let bedline = infile.readLine()
      let bedfields = bedline.split(sep='\t')
      let chr = bedfields[0]
      let startIdx = parseUInt(bedfields[1]) + 1'u
      let endIdx = parseUInt(bedfields[2])
      let txid = bedfields[3]
      let strand = bedfields[5]
      if stranded:
        outfile.write(&"{chr}\tBLANK\ttranscript\t{startIdx}\t{endIdx}\t.\t" &
          &"{strand}\t.\ttranscript_id \"{txid}\";\n")
      else:
        outfile.write(&"{chr}\tBLANK\ttranscript\t{startIdx}\t{endIdx}\t.\t" &
          &".\t.\ttranscript_id \"{txid}\";\n")
      let blockSizes = bedfields[10].split(sep=',')
      let blockStarts= bedfields[11].split(sep=',')

      # var exon_count = 1
      # var reverse_exon_count = blockStarts.len
      for i in 0..<blockStarts.len:
        var newStartIdx, newEndIdx, exonNumber : uint
        if strand == "+" or not stranded:
          newStartIdx = startIdx + parseUInt(blockStarts[i])
          newEndIdx = startIdx +
                      parseUInt(blockStarts[i]) +
                      parseUInt(blockSizes[i]) - 1
          exonNumber = uint(i + 1)
        elif strand == "-":
          newStartIdx = startIdx +
                        parseUInt(blockStarts[i])
          newEndIdx = startIdx +
                      parseUInt(blockStarts[i]) +
                      parseUInt(blockSizes[i]) - 1
          exonNumber = uint(blockStarts.len - i)
        if stranded:
          outfile.write(&"{chr}\tBLANK\texon\t{newStartIdx}\t{newEndIdx}\t." &
            &"\t{strand}\t.\ttranscript_id \"{txid}\"; " &
            &"exonNumber \"{exonNumber}\"\n")
        else:
          outfile.write(&"{chr}\tBLANK\texon\t{newStartIdx}\t{newEndIdx}\t." &
            &"\t.\t.\ttranscript_id \"{txid}\"; " &
            &"exonNumber \"{exonNumber}\"\n")
  except EOFError:
    discard
  infile.close()
  outfile.close()


proc strandTranscripts*(transcripts : openArray[FastaRecord],
                        outfilepath : string,
                        wrap_len : int = 60) =
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
    outfile.write(&">{transcript.readId}\n")
    for i in 0..<(stranded.len div wrap_len):
      outfile.write(&"{stranded[i*wrap_len..(i+1)*wrap_len - 1]}\n")
    if stranded.len mod wrap_len != 0 :
      outfile.write(&"{stranded[wrap_len*(stranded.len div wrap_len)..^1]}\n")
  outfile.close()


proc strandTranscripts*(transcripts : openArray[FastqRecord],
                        outfilepath : string,
                        wrap_len : int = 60) =
  translateTranscripts(convertFASTQtoFASTA(transcripts),
                       outfilepath,
                       wrap_len)


proc compareExactTranslations*(referenceInfilepath : string,
                               translation_infilepath : string) =
  var rInfile, tInfile : File
  discard open(rInfile,referenceInfilepath,fmRead)
  discard open(tInfile,translation_infilepath,fmRead)
  let rRecords = parseFasta(rInfile)
  rInfile.close()
  let tRecords = parseFasta(tInfile)
  tInfile.close()
  var rProteins,tProteins : HashSet[string]
  for record in rRecords:
    rProteins.incl(record.sequence)
  for record in tRecords:
    tProteins.incl(record.sequence)
  let tp = intersection(rProteins,tProteins).len
  let fp = difference(tProteins,rProteins).len
  let fn = difference(rProteins,tProteins).len
  echo "TP: ", tp
  echo "FP: ", fp
  echo "FN: ", fn
  echo ""
  echo &"Precision: {float(tp) / float(tp+fp)}"
  echo &"Recall:    {float(tp) / float(tp+fn)}"


proc compareBLASTPTranslations*(referenceInfilepath : string,
                                blastpInfilepath : string) =
  var fp,tp1 = 0
  var referenceIdSet : HashSet[string]
  var matchSet : HashSet[string]
  var refInfile : File
  discard open(refInfile,referenceInfilepath,fmRead)
  let referenceRecords = parseFasta(refInfile)
  refInfile.close()
  for record in referenceRecords:
    referenceIdSet.incl(record.readId)
  
  let blastRecords = parseBLASTPoutput(blastpInfilepath)
  for record in blastRecords:
    if record.matchNames.len > 0:
      matchSet.incl(record.matchNames[0])
      tp1 += 1
    else:
      fp += 1

  let tp2 = matchSet.len
  let fn = difference(referenceIdSet,matchSet).len

  echo "TP1 - Total matches: ", tp1
  echo "TP2 - Match set len: ", tp2
  echo "FP  - Non matched query: ", fp
  echo "FN  - Non matched reference: ", fn
  echo ""
  echo &"Precision (uses TP1): {float(tp1) / float(tp1+fp)}"
  echo &"Recall (uses TP2):    {float(tp2) / float(tp2+fn)}"


proc extractIntronsFromBED12*(infilepath : string,
                              outfilepath : string) = 
  ## Grabs Introns from BED12 file and reports one per line in BED format
  ## with ID = ID from the BED12 line
  var infile,outfile : File
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  try:
    while true:
      let bedline = infile.readLine()
      let bedfields = bedline.split(sep='\t')
      let chr = bedfields[0]
      let startIdx = parseUInt(bedfields[1])
      # let endIdx = parseUInt(bedfields[2])
      let txid = bedfields[3]
      let strand = bedfields[5]
      let blockSizes = bedfields[10].split(sep=',')
      let blockStarts= bedfields[11].split(sep=',')

      # var exon_count = 1
      # var reverse_exon_count = blockStarts.len
      for i in 1..<blockStarts.len:
        let intronStartIdx = startIdx +
                             parseUInt(blockStarts[i-1]) +
                             parseUInt(blockSizes[i-1])
        let intronEndIdx = startIdx +
                           parseUInt(blockStarts[i])
        outfile.write(&"{chr}\t{intronStartIdx}\t{intronEndIdx}\t{txid}\t." &
          &"\t{strand}\n")
  except EOFError:
    discard
  infile.close()
  outfile.close()


proc testOverlap(read : (uint64,uint64),
                 single_exon_gene_list : seq[(uint64,uint64,string)]) : string =
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


proc callNonCanonicalSplicingFromFASTA*(infilepath : string,
                                        outfilepath : string) =
  var infile,outfile : File
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  let records = parseFasta(infile)
  infile.close()
  for record in records:
    if record.sequence.len < 2:
      echo "WARNING - Very short intron detected"
    else:
      if record.sequence[0..1].toUpperAscii() != "GT" or
         record.sequence[^2..^1].toUpperAscii() != "AG":
        outfile.writeLine(record.readId)
  outfile.close()


proc splitFASTAByReadCounts*(infilepath : string,
                             outfile_prefix : string,
                             bins : openArray[uint64] = [1'u64,
                                                         2'u64,
                                                         5'u64,
                                                         10'u64,
                                                         20'u64,
                                                         40'u64,
                                                         80'u64,
                                                         160'u64,
                                                         320'u64,
                                                         640'u64]) = 
  var infile : File
  discard open(infile,infilepath, fmRead)
  let fastaRecords = parseFasta(infile)
  infile.close()
  var splitRecords : seq[seq[FastaRecord]]
  for i in 0..<bins.len:
    splitRecords.add(@[])
  for record in fastaRecords:
    let numReads = uint64(parseUInt(record.readId.split('_')[^1]))
    var bin = -1
    for i in 0..<(bins.len - 1):
      if numReads >= bins[i] and numReads < bins[i+1]:
        bin = i
    if numReads >= bins[^1]:
      bin = bins.len - 1
    if bin != -1:
      splitRecords[bin].add(record)
  for i,split in splitRecords:
    var binId = ""
    if i == bins.len - 1:
      binId = &"{bins[^1]}+"
    else:
      binId = &"{bins[i]}-{bins[i+1] - 1}"
    var outfile : File
    discard open(outfile,&"{outfile_prefix}_{binId}.fa",fmWrite)
    for record in split:
      outfile.write(&">{record.readId}\n")
      outfile.writeLine(record.sequence)
    outfile.close()


proc filterFASTAByReadCounts*(infilepath,outfilepath : string,
                              filter : uint64 = 5'u64) =
  var infile,outfile : File
  discard open(infile,infilepath, fmRead)
  let fastaRecords = parseFasta(infile)
  infile.close()
  discard open(outfile,outfilepath,fmWrite)
  for record in fastaRecords:
    let numReads = uint64(parseUInt(record.readId.split('_')[^1]))
    if numReads >= filter:
      outfile.write(&">{record.readId}\n")
      outfile.writeLine(record.sequence)
  outfile.close()


proc parseAttributes(s : string) : Table[string,string] = 
  let fields = s.split(';')[0..^2]
  for field in fields:
    let fields1 = field.strip(leading=true,trailing=true,chars={' '}).split(' ')
    let key = fields1[0]
    let val = fields1[1].strip(leading=true,trailing=true,chars = {'"'})
    result[key] = val


proc callNovelNonCanonical(referenceInfilepath,
                           infilepath,
                           outfilepath : string,
                           threshold : uint = 5) =
  var referenceInfile,infile,outfile : File
  discard open(referenceInfile,referenceInfilepath,fmRead)
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  var referenceIntrons : HashSet[(string,uint64,uint64)]
  var lastExons : Table[string,(string,uint64,uint64)]
  try:
    while true:
      let line = referenceInfile.readLine()
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
        let newChr = fields0[0]
        let newStartIdx =  uint64(parseUInt(fields0[3])) - 1'u32
        let newEndIdx = uint64(parseUInt(fields0[4]))
        if attributes["transcript_id"] in lastExons:
          let (chr, startIdx,endIdx)= lastExons[attributes["transcript_id"]]
          try:
            assert chr == newChr
            referenceIntrons.incl((chr, endIdx, newStartIdx))
          except AssertionDefect:
            echo "ERROR - same transcript on different chromosomes:"
            let txId = attributes["transcript_id"]
            echo &"{chr}:{startIdx}-{endIdx} {txId}"
            echo &"{newChr}:{newStartIdx}-{newEndIdx} {txId}"
          # outfilepath.writeLine(chr,"\t",endIdx,"\t"," & 
          #   &"newStartIdx,"\t","reference_intron")
        lastExons[attributes["transcript_id"]] = (newChr,
                                                  newStartIdx,
                                                  newEndIdx)
      else:
        echo "ERROR - no field transcript_id"
  except EOFError:
    discard
  referenceInfile.close()
  var novelCounter = 0
  var totalCounter = 0
  try:
    while true:
      let line = infile.readLine()
      if line.len == 0:
        continue
      let fields0 = line.split(':')
      let clusterId = fields0[0]
      let readSupport = parseUInt(clusterId.split('_')[^1])
      let chr = fields0[2]
      let indices = fields0[3].strip(chars={'(',')','+','-'})
      let splitIndices = indices.split('-')
      let startIdx = uint64(parseUInt(splitIndices[0])) 
      let endIdx = uint64(parseUInt(splitIndices[1]))
      if readSupport >= threshold:
        totalCounter += 1
        if (chr,startIdx,endIdx) notin  referenceIntrons:
          outfile.writeLine(&"{chr}\t{startIdx}\t{endIdx}\tquery_intron")
          novelCounter += 1
  except EOFError:
    discard
  infile.close()
  outfile.close()
  echo &"Total introns above threshold - {totalCounter}"
  echo &"Novel introns - {novelCounter}"


proc callOverlappingNonCanonical(referenceInfilepath,
                                 infilepath,
                                 outfilepath : string) =
  var referenceInfile,infile,outfile : File
  discard open(referenceInfile,referenceInfilepath,fmRead)
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  var referenceIntrons : HashSet[(string,uint64,uint64)]
  # var lastExons : Table[string,(string,uint64,uint64)]
  try:
    while true:
      let line = referenceInfile.readLine()
      if line.len == 0:
        continue
      let fields0 = line.split(':')
      let chr = fields0[2]
      let indices = fields0[3].strip(chars={'(',')','+','-'})
      let splitIndices = indices.split('-')
      let startIdx = uint64(parseUInt(splitIndices[0])) 
      let endIdx = uint64(parseUInt(splitIndices[1]))
      referenceIntrons.incl((chr,startIdx,endIdx))
  except EOFError:
    discard
  referenceInfile.close()
  var novelCounter = 0
  var overlappingCounter = 0
  try:
    while true:
      let line = infile.readLine()
      if line.len == 0:
        continue
      let fields0 = line.split(':')
      let chr = fields0[2]
      let indices = fields0[3].strip(chars={'(',')','+','-'})
      let splitIndices = indices.split('-')
      let startIdx = uint64(parseUInt(splitIndices[0])) 
      let endIdx = uint64(parseUInt(splitIndices[1]))
      if (chr,startIdx,endIdx) notin  referenceIntrons:
        novelCounter += 1
      else:
        outfile.writeLine(line)
        overlappingCounter += 1
  except EOFError:
    discard
  infile.close()
  outfile.close()
  echo &"Non-overlapping introns - {novelCounter}"
  echo &"Overlapping introns - {overlappingCounter}"


proc assignTxIDs(referenceInfilepath,
                 infilepath,
                 outfilepath : string) = 
  var referenceInfile,infile,outfile : File
  discard open(referenceInfile,referenceInfilepath,fmRead)
  var txExons : Table[(string,string),seq[(uint64,uint64)]]
  try:
    while true:
      let line = referenceInfile.readLine()
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
        let newChr = fields0[0]
        let newStartIdx =  uint64(parseUInt(fields0[3])) - 1'u32
        let newEndIdx = uint64(parseUInt(fields0[4]))
        if (attributes["transcript_id"],newChr) in txExons:
          txExons[(attributes["transcript_id"],newChr)].add((newStartIdx,
                                                             newEndIdx))
        else:
          txExons[(attributes["transcript_id"],newChr)] = @[(newStartIdx,
                                                             newEndIdx)]
      else:
        echo "ERROR - no field transcript_id"
  except EOFError:
    discard
  referenceInfile.close()
  var txIntrons : Table[(string,seq[(uint64,uint64)]), string]
  var singleExonGenes : Table[string,seq[(uint64,uint64,string)]]
  for (txId,chr) in txExons.keys:
    # echo chr
    let exonChain = txExons[(txId,chr)]
    if exonChain.len > 1:
      var intronChain : seq[(uint64,uint64)]
      for i in 1..<exonChain.len:
        intronChain.add((exonChain[i-1][1],exonChain[i][0]))
      txIntrons[(chr,intronChain)] = txId
    else:
      if chr in singleExonGenes:
        singleExonGenes[chr].add((exonChain[0][0],exonChain[0][1],txId))
      else:
        singleExonGenes[chr] = @[(exonChain[0][0],exonChain[0][1],txId)]
  for chr in singleExonGenes.keys:
    singleExonGenes[chr].sort
  
  #Read in BED12 file and get intron chains
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  try:
    while true:
      let bedline = infile.readLine()
      let bedfields = bedline.split(sep='\t')
      let chr = bedfields[0]
      let startIdx = parseUInt(bedfields[1])
      let endIdx = parseUInt(bedfields[2])
      let txid = bedfields[3]
      # let strand = bedfields[5]
      let blockSizes = bedfields[10].split(sep=',')
      let blockStarts= bedfields[11].split(sep=',')
      var referenceId = ""
      if blockStarts.len > 1:
        var intronChain : seq[(uint64,uint64)]
        for i in 1..<blockStarts.len:
          let intronStartIdx = uint64(startIdx +
                                      parseUInt(blockStarts[i-1]) +
                                      parseUInt(blockSizes[i-1]))
          let intronEndIdx = uint64(startIdx +
                                    parseUInt(blockStarts[i]))
          intronChain.add((intronStartIdx,intronEndIdx))
        if (chr,intronChain) in txIntrons:
          referenceId = txIntrons[(chr,intronChain)]
      else:
        referenceId = testOverlap((uint64(startIdx),uint64(endIdx)),
                                  singleExonGenes[chr])
      if referenceId == "":
        outfile.write(&"{txid}\t.\n")
      else:
        outfile.write(&"{txid}\t{referenceId}\n")
  except EOFError:
    discard
  infile.close()
  outfile.close()


proc parseGTF*(infile : File) : HashSet[(string,char,seq[(uint64,uint64)])] = 
  var exons : Table[string,seq[(string,char,uint64,uint64)]]
  try:
    while true:
      let line = infile.readLine()
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
        let newChr = fields0[0]
        let newStrand = fields0[6][0]
        let newStartIdx =  uint64(parseUInt(fields0[3])) - 1'u32
        let newEndIdx = uint64(parseUInt(fields0[4]))
        if attributes["transcript_id"] in exons:
          let (chr, strand, startIdx,endIdx) =
            exons[attributes["transcript_id"]][^1]
          try:
            assert chr == newChr
            assert strand == newStrand
            # referenceIntrons.incl((chr, endIdx, newStartIdx))
          except AssertionDefect:
            echo "ERROR - same transcript on different strand or chromosomes:"
            let txId = attributes["transcript_id"]
            echo &"{chr}({strand}):{startIdx}-{endIdx} {txId}"
            echo &"{newChr}({newStrand}):{newStartIdx}-{newEndIdx} {txId}"
          exons[attributes["transcript_id"]].add((newChr,
                                                  newStrand,
                                                  newStartIdx,
                                                  newEndIdx))
        else:
          exons[attributes["transcript_id"]] = @[(newChr,
                                                  newStrand,
                                                  newStartIdx,
                                                  newEndIdx)]
      else:
        echo "ERROR - no field transcript_id"
  except EOFError:
    discard

  for txId in exons.keys:
    let totalChr = exons[txId][0][0]
    let totalStrand = exons[txId][0][1]
    # let totalStart = exons[txId][0][2]
    # let totalEnd = exons[txId][^1][3]
    var introns : seq[(uint64,uint64)]
    for i in 1..<exons[txId].len:
      let (_,_,_,endIdx1) = exons[txId][i-1]
      let (_,_,startIdx2,_) = exons[txId][i]
      introns.add((endIdx1,startIdx2))
    # result.incl((totalChr,totalStart,totalEnd,totalStrand,introns))
    result.incl((totalChr,totalStrand,introns))

#TODO - Add to help
#TODO - (?) add outfile writing (?)
proc idNovelIsoforms*(infilepath,
                      reference1_infilepath,
                      reference2_infilepath,
                      outfilepath : string) =
  # var infile, ref1file, ref2file, outfile : File
  var infile, ref1file, ref2file : File
  discard open(ref1file,reference1_infilepath,fmRead)
  let ref1set = parseGTF(ref1file)
  ref1file.close
  
  discard open(ref2file,reference2_infilepath,fmRead)
  let ref2set = parseGTF(ref2file)
  ref2file.close
  
  discard open(infile,infilepath,fmRead)
  let inset = parseGTF(infile)
  infile.close

  let novel = (inset - ref1set) - ref2set
  echo "Total in infile: ", inset.len
  echo "Total in ref1:   ", ref1set.len
  echo "Total in ref2:   ", ref2set.len
  echo "infile - ref1:   ", (inset - ref1set).len
  echo "infile - ref2:   ", (inset - ref2set).len
  echo "ref1 - ref2:     ", (ref1set - ref2set).len
  echo "Novel in infile: ", novel.len

  # discard open(outfile,outfilepath,fmWrite)
  # outfile.close


proc getTxId*(s : string) : string =
  result = s.split(':')[1].split('|')[0]


proc getNovelLociFASTA*(infilepath,
                        gffcompare_infilepath,
                        outfilepath : string,
                        field = 1) = 
  var novelLoci : HashSet[string]
  var infile,gfffile,outfile : File
  discard open(gfffile,gffcompare_infilepath,fmRead)
  try:
    while true:
      let line = gfffile.readLine()
      if line.len == 0:
        continue
      elif line[0] == '#':
        continue
      let fields0 = line.split('\t')
      if fields0[3] != "u" and
         fields0[3] != "p" and
         fields0[3] != "y" and
         fields0[3] != "i" and
         fields0[3] != "x" and
         fields0[3] != "s":
        continue
      if fields0[3+field] != "-":
        echo fields0[1]
        novelLoci.incl(getTxId(fields0[3+field]))
        # for txId in fields0[2+field].split(','):
          # novelLoci.incl(txId)
  except EOFError:
    discard
  gfffile.close
  
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  let records = parseFasta(infile)
  var newRecords : seq[FastaRecord]
  for record in records:
    if record.readId in novelLoci:
      newRecords.add(record)
  infile.close
  # writeCorrectedReads(newRecords,outfile)
  writeFASTArecordsToFile(outfile,newRecords)
  outfile.close


proc parseOptions() : UtilOptions = 
  var i = 0
  var mode = ""
  var last = ""

  var helpFlag = true
  var runFlag = true

  var fastq = false
  var fastqFlag = false

  var minLength = 75'u64
  var minLengthFlag = false

  var infilepath = ""
  var infilepathFlag = false

  var stranded = false

  var referenceInfilepath = ""
  var referenceInfilepathFlag = false
  
  var referenceInfilepath2 = ""
  var referenceInfilepath2Flag = false

  var outfilepath = ""
  var outfilepathFlag = false

  for kind, key, val in getopt():
    # echo kind," ", key," ", val
    if i == 0:
      if kind == cmdArgument:
        mode = key
        i += 1
        continue
      else:
        if key == "h" or key == "help":
          mode = "help"
        helpFlag = true
        runFlag = false
        break

    if i == 1:
      helpFlag = false
    i += 1
    case mode:
      of "translate",
          "strandTranscripts",
          "bed2gtf",
          "parseBLASTP",
          "compareBLASTP",
          "compareFASTA",
          "splitFASTA",
          "filterFASTA",
          "extractIntrons",
          "callNonCanonical",
          "callNovelNonCanonical",
          "callOverlapping",
          "assignIDs",
          "idNovelIsoforms",
          "getNovelLociFASTA":
        case kind:
          of cmdEnd:
            break
          of cmdShortOption, cmdLongOption:
            if last != "":
              echo &"ERROR - Option {last} provided without an argument"
              runFlag = false
              break
            case key:
              of "l", "min-length":
                if not minLengthFlag:
                  minLengthFlag = true
                  if val != "":
                    minLength = parseUInt(val)
                  else:
                    last = "min-length"
                else:
                  echo "ERROR - Multiple min-length provided"
                  runFlag = false
              of "i":
                if not infilepathFlag:
                  infilepathFlag = true
                  if val != "":
                    infilepath = val
                  else:
                    last = "infile"
                else:
                  echo "ERROR - Multiple infiles provided"
                  runFlag = false
                  break
              of "q", "fastq":
                if not fastqFlag:
                  fastq = true
                  fastqFlag = true
                elif not fastq:
                  echo "ERROR - Conflicting flags -q / --fastq and -a / --fasta"
                  runFlag = false
                  break
              of "a", "fasta":
                if not fastqFlag:
                  fastqFlag = true
                elif fastq:
                  echo "ERROR - Conflicting flags -q / --fastq and -a / --fasta"
                  runFlag = false
                  break
              of "o":
                if not outfilepathFlag:
                  outfilepathFlag = true
                  if val != "":
                    infilepath = val
                  else:
                    last = "outfile"
                else:
                  echo "ERROR - Multiple outfiles provided"
                  runFlag = false
                  break
              of "r":
                if not referenceInfilepathFlag:
                  referenceInfilepathFlag = true
                  if val != "":
                    referenceInfilepath = val
                  else:
                    last = "reference"
                else:
                  echo "ERROR Multiple references provided"
                  runFlag = false
                  break
              of "r2":
                if not referenceInfilepath2Flag:
                  referenceInfilepath2Flag = true
                  if val != "":
                    referenceInfilepath2 = val
                  else:
                    last = "reference2"
                else:
                  echo "ERROR Multiple references provided"
                  runFlag = false
                  break
              of "s", "stranded":
                stranded = true
              of "h", "help":
                helpFlag = true
                runFlag = false
                break
          of cmdArgument:
            case last:
              of "min-length":
                minLength = parseUInt(key)
              of "infile":
                infilepath = key
              of "outfile":
                outfilepath = key
              of "reference":
                referenceInfilepath = key
              of "reference2":
                referenceInfilepath2 = key
              else:
                echo &"ERROR - unknown option {last} provided"
                runFlag = false
                break
            last = ""
      else:
        helpFlag = true
        break
  if infilepath == "" and not helpFlag:
    echo "ERROR - infilepath must be specified"
    runFlag = false
    helpFlag = true
  var requireOutfile,requireReference = false
  case mode:
    of "translate",
       "bed2gtf",
       "parseBLASTP",
       "splitFASTA",
       "filterFASTA",
       "extractIntrons":
      requireOutfile = true
    of "compareBLASTP",
       "compareFASTA":
      requireReference = true
    of "callNonCanonical",
       "callNovelNonCanonical",
       "callOverlapping":
      requireReference = true
      requireOutfile = true
  
  if requireOutfile and outfilepath == "" and not helpFlag:
    echo "ERROR - outfilepath must be specified"
    runFlag = false
    helpFlag = true
  if requireReference and referenceInfilepath == "" and not helpFlag:
    echo "ERROR - referencefilepath must be specified"
    runFlag = false
    helpFlag = true
  if helpFlag:
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
      of "help":
        writeDefaultHelp()
        quit(QuitSuccess)
      else:
        echo "ERROR - first argument must specify utility function"
        writeDefaultHelp()
  return UtilOptions(mode : mode,
                     runFlag : runFlag,
                     infilepath : infilepath,
                     outfilepath : outfilepath,
                     referenceInfilepath : referenceInfilepath,
                     referenceInfilepath2 : referenceInfilepath2,
                     minLength : minLength,
                     fastq : fastq,
                     stranded : stranded)


proc main() =
  let opt = parseOptions()
  if opt.runFlag:
    case opt.mode:
      of "translate": #TODO - Wrap this in a proc
        var infile : File
        discard open(infile,opt.infilepath,fmRead)
        if opt.fastq:
          let records = parseFASTQ(infile)
          infile.close()
          translateTranscripts(records,
                               opt.outfilepath,
                               threshold = int(opt.minLength),
                               stranded = opt.stranded)
        else:
          let records = parseFasta(infile)
          infile.close()
          translateTranscripts(records,
                               opt.outfilepath,
                               threshold = int(opt.minLength),
                               stranded = opt.stranded)
      of "strandTranscripts": #TODO - Wrap this in a proc
        var infile : File
        discard open(infile,opt.infilepath,fmRead)
        if opt.fastq:
          let records = parseFASTQ(infile)
          infile.close()
          strandTranscripts(records,opt.outfilepath)
        else:
          let records = parseFasta(infile)
          infile.close()
          strandTranscripts(records,opt.outfilepath)
      of "bed2gtf":
        convertBED12toGTF(opt.infilepath,opt.outfilepath,opt.stranded)
      of "parseBLASTP": # TODO - Wrap this into a proc
        let blastOutput = parseBLASTPoutput(opt.infilepath)
        var outfile : File
        discard open(outfile,opt.outfilepath,fmWrite)
        for blast_match in blastOutput:
          if blast_match.identitiesPercentages.len != 0:
            outfile.write(&"{blast_match.queryName}\t" &
              &"{blast_match.matchNames[0]}\t{blast_match.evals[0][0]}\n")
          else:
            outfile.write(&"{blast_match.queryName}\t----\t----\n")
        outfile.close()
      of "compareBLASTP":
        compareBLASTPTranslations(opt.referenceInfilepath,
                                  opt.infilepath)
      of "compareFASTA":
        compareExactTranslations(opt.referenceInfilepath,
                                 opt.infilepath)
      of "splitFASTA":
        splitFASTAByReadCounts(opt.infilepath,
                               opt.outfilepath)
      of "filterFASTA":
        filterFASTAByReadCounts(opt.infilepath,
                                opt.outfilepath)
      of "extractIntrons":
        extractIntronsFromBED12(opt.infilepath,
                                opt.outfilepath)
      of "callNonCanonical":
        callNonCanonicalSplicingFromFASTA(opt.infilepath,
                                          opt.outfilepath)
      of "callNovelNonCanonical":
        callNovelNonCanonical(opt.referenceInfilepath,
                              opt.infilepath,
                              opt.outfilepath)
      of "callOverlapping":
        callOverlappingNonCanonical(opt.referenceInfilepath,
                                    opt.infilepath,
                                    opt.outfilepath)
      of "assignIDs":
        assignTxIDs(opt.referenceInfilepath,
                    opt.infilepath,
                    opt.outfilepath)
      of "idNovelIsoforms":
        idNovelIsoforms(opt.infilepath,
                        opt.referenceInfilepath,
                        opt.referenceInfilepath2,
                        opt.outfilepath)
      of "getNovelLociFASTA":
        getNovelLociFASTA(opt.infilepath,
                          opt.referenceInfilepath,
                          opt.outfilepath)



main()