import os
import osproc
import strutils
import strformat
import tables
import poGraphUtils

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
    min_length : uint

#TODO: Make these functions callable from the command line, such that other people can use these utils.


proc conduitUtilsVersion() : string =
  return "CONDUIT Utilities Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu, https://github.com/NatPRoach/conduit/ )"

proc writeDefaultHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "Usage:"
  echo "  ./conduitUtils <translate | bed2gtf | parseBLASTP>"

proc writeTranslateHelp() =
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "Usage:"
  echo "  ./conduitUtils translate [options] -i <transcripts.fa> -o <predicted_protein.fa>"
  echo "  <transcripts.fa>         FASTA infile containing putative transcripts to be translated"
  echo "  <predicted_protein.fa>   FASTA outfile containing in silico translated ORFs from transcripts.fa"
  echo ""
  echo "Options (defaults in parentheses):"
  echo "  Filtering Options:"
  echo "    -l, --min-length (75)"
  echo "        Minimum length in Amino Acids necessary for a putative ORF to be reported"

proc writeBED2GTFHelp() =
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "Usage:"
  echo "  ./conduitUtils bed2gtf -i <infile.bed> -o <outfile.gtf>"
  echo "  <infile.bed>    BED12 infile to be converted in to GTF format"
  echo "  <outfile.gtf>   GTF outfile"

proc writeBLASTPHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitUtilsVersion()
  echo "Usage:"
  echo "  ./conduitUtils parseBLASTP -i <inBLASTP.txt> -o <outPutativeOrthologs.tsv>"
  echo "  <inBLASTP.txt>              Default output of BLASTP search of translated protein products vs some reference proteome"
  echo "  <outPutativeOrthologs.tsv>  Tab separated file of putative ortholog matches"
  echo "                              In format: <Query ID>\t<Reference proteome top match ID>\t<E value>"

proc parseBLASTPoutput*(infilepath : string) : seq[BLASTmatch] = 
  ##
  ## Takes in BLASTP default output file and outputs seq of matching proteins
  ##
  # Discard the header / citation info
  var infile : File
  discard open(infile,infilepath,fmRead)

  for _ in 0..<23:
    discard infile.readLine()
  while true:
    try:
      let query_line = infile.readLine() # Query= cluster_10_0 <unknown description>
      # echo query_line
      if query_line[0..11] == "  Database: ":
        break
      assert query_line[0..6] == "Query= "
      let query_name = query_line[7..^1]

      discard infile.readLine() # ""

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
                            "AGT" : 'T',
                            "AGC" : 'T',
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

proc translateTranscript*(nts : string) : string =
  let start_codon_indices = findAll(nts,"ATG")
  result = ""
  for start_codon_index in start_codon_indices:
    let translation = translateORF(nts[start_codon_index..^1])
    if translation.len > result.len:
      result = translation

proc translateTranscripts*(transcripts : openArray[FastaRecord],outfilepath : string , threshold : int = 75,wrap_len : int = 60) =
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  for transcript in transcripts:
    let translation = translateTranscript(transcript.sequence)
    if translation.len >= threshold:
      outfile.write(&">{transcript.read_id}\n")
      for i in 0..<(translation.len div wrap_len):
        outfile.write(&"{translation[i*wrap_len..(i+1)*wrap_len - 1]}\n")
      if translation.len mod wrap_len != 0 :
        outfile.write(&"{translation[wrap_len*(translation.len div wrap_len)..^1]}\n")
  outfile.close()

proc convertBED12toGTF*(infilepath : string, outfilepath : string ) =
  ## Converts BED12 formatted file to well-formed GTF file suitable for evaluation with GFFcompare
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
      outfile.write(&"{chr}\tBLANK\ttranscript\t{start_idx}\t{end_idx}\t.\t{strand}\t.\ttranscript_id \"{txid}\";\n")
      let block_sizes = bedfields[10].split(sep=',')
      let block_starts= bedfields[11].split(sep=',')

      # var exon_count = 1
      # var reverse_exon_count = block_starts.len
      for i in 0..<block_starts.len:
        var new_start_idx, new_end_idx, exon_number : uint
        if strand == "+":
          new_start_idx = start_idx + parseUInt(block_starts[i])
          new_end_idx = start_idx + parseUInt(block_starts[i]) + parseUInt(block_sizes[i]) - 1
          exon_number = uint(i + 1)
        elif strand == "-":
          new_start_idx = start_idx + parseUInt(block_starts[i])
          new_end_idx = start_idx + parseUInt(block_starts[i]) + parseUInt(block_sizes[i]) - 1
          exon_number = uint(block_starts.len - i)
        outfile.write(&"{chr}\tBLANK\ttranscript\t{new_start_idx}\t{new_end_idx}\t.\t{strand}\t.\ttranscript_id \"{txid}\"; exon_number \"{exon_number}\"\n")
  except EOFError:
    discard
  infile.close()
  outfile.close()


proc parseOptions() : UtilOptions = 
  
  var i = 0
  var mode = ""
  var last = ""

  var min_length = 75'u64
  var min_length_flag = false

  var infilepath = ""
  var infilepath_flag = false
  
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
      of "translate", "bed2gtf", "parseBLASTP":
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
              of "-i":
                if not infilepath_flag:
                  infilepath_flag = true
                  if val != "":
                    infilepath = val
                  else:
                    last = "infile"
                else:
                  echo "ERROR - Multiple infiles provided"
                  run_flag = false
              of "-o":
                if not outfilepath_flag:
                  outfilepath_flag = true
                  if val != "":
                    infilepath = val
                  else:
                    last = "outfile"
                else:
                  echo "ERROR - Multiple outfiles provided"
                  run_flag = false
          of cmdArgument:
            case last:
              of "min-length":
                min_length = parseUInt(key)
              of "infile":
                infilepath = key
              of "outfile":
                outfilepath = key
              else:
                echo &"ERROR - unknown option {last} provided"
                run_flag = false
                break
            last = ""
      else:
        help_flag = true
        break
  if infilepath == "":
    echo "ERROR - infilepath must be specified"
    run_flag = false
    help_flag = true
  if outfilepath == "":
    echo "ERROR - outfilepath must be specified"
    run_flag = false
    help_flag = true
  if help_flag:
    case mode:
      of "translate":
        writeTranslateHelp()
      of "bed2gtf":
        writeBED2GTFHelp()
      of "parseBLASTP":
        writeBLASTPHelp()
      else:
        echo "ERROR - first argument must specify utility function: \"translate\", \"bed2gtf\", or \"parseBLASTP\""
        writeDefaultHelp()
  return UtilOptions(mode : mode,
                     run_flag : run_flag
                     infilepath : infilepath,
                     outfilepath : outfilepath,
                     min_length : min_length
                     )
proc main() =
  let opt = parseOptions()
  if opt.run_flag:
    case opt.mode:
      of "translate":
        var infile : File
        discard open(infile,infilepath,fmRead)
        records = poGraphUtils.parseFasta(infile)
        infile.close()
        translateTranscripts(records,opt.outfilepath,threshold = int(opt.min_length))
      of "bed2gtf":
        convertBED12toGTF(opt.infilepath,opt.outfilepath)
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



