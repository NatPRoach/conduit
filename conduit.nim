import os
import osproc
import parseopt
import strutils
import strformat
import poParser
import tables
import threadpool_simple as tps
import hts
import poaV2/header
import poaV2/poa
{.experimental.}

type
  ConduitOptions = object
    run_flag : bool
    final_polish : bool
    intermediates : bool
    mode : string
    nanopore_format : string
    bowtie_strand_constraint : string
    bowtie_alignment_mode : string
    bowtie_read_inputs : string
    score_matrix_path : string
    output_dir : string
    tmp_dir : string
    files : seq[string]
    trims : seq[string]
    isoform_delta : uint64
    ends_delta : uint64
    illumina_weight : uint64
    thread_num : uint64
    max_iterations : uint64


proc conduitVersion() : string =
  return "CONDUIT Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu, https://github.com/NatPRoach/conduit/ )"

proc writeDefaultHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitVersion()
  echo "Usage:"
  echo "  ./conduit <nano | hybrid>"

proc writeNanoHelp() =
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitVersion()
  echo "Usage:"
  echo "  ./conduit nano [options] <clusters_directory>"
  echo "  <clusters_directory>   Directory containing the .fasta/.fa or .fastq/.fq files of reads separated by gene cluster"
  echo ""
  echo "Options (defaults in parentheses):"
  echo "  Scaffold type:"
  echo "    --drna (default)"
  echo "        Scaffold reads are stranded forward relative to coding strand, and may contain U characters instead of Ts"
  echo "    --cdna-rev-stranded"
  echo "        Scaffold reads are stranded reverse complemented relative to coding strand"
  echo "    --cdna"
  echo "        Scaffold reads are NOT stranded"
  echo "    --sfq (default)"
  echo "        Scaffold reads are in FASTQ format"
  echo "    --sfa"
  echo "        Scaffold reads are in FASTA format"
  echo "  Consensus Collapsing:"
  echo "    -d, --isoform-delta (35)"
  echo "        Maximum indel size to be 'corrected', beyond this size a new isoform is declared"
  echo "    -e, --ends-delta (35)"
  echo "        Maximum size at the ends of isoforms to 'correct' before splitting"
  echo "  Ouput:"
  echo "    -o, --output-dir <path> (conduit/)"
  echo "        <path> where corrected clusters will be written"
  echo "        NOTE: THIS WILL OVERWRITE EXISTING FILES!"
  echo "  Miscellaneous:"
  echo "    -h, --help"
  echo "        Display this help message and exit"
  echo "    -v, --version"
  echo "        Display the installed version number of CONDUIT and exit"
  echo "    --tmp-dir <path> (conduit-tmp/)"
  echo "        <path> where temporary files will be created"
  echo "    -t, --threads (4)"
  echo "        Number of threads to run in parallel (used for both Bowtie2 and Partial Order Graph correction)"
  #echo "    --cdna" #Not supported yet #TODO

proc writeHybridHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitVersion()
  echo "Usage:"
  echo "  ./conduit hybrid [options] <clusters_directory> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>}"
  echo "  <clusters_directory>   Directory containing the .fasta/.fa or .fastq/.fq files of reads separated by gene cluster"
  echo "                         NOTE: .gz support coming for nanopore scaffold data, but is not an option at this time" #TODO
  echo ""
  echo "  Illumina data is aligned with Bowtie2, therefore Illumina data is provided in the same format as Bowtie2, namely:"
  echo ""
  echo "    <m1>                   Files with #1 mates, paired with files in <m2>"
  echo "                           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)."
  echo "    <m2>                   Files with #2 mates, paired with files in <m1>"
  echo "                           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)."
  echo "    <r>                    Files with unpaired reads"
  echo "                           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)."
  echo "    <i>                    File with interleaved paired-end FASTQ/FASTA reads"
  echo "                           Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)."
  echo "    <bam>                  Files are unaligned BAM sorted by read name."
  echo ""
  echo "  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be specified many times."
  echo "  E.g. '-U file1.fq,file2.fq -U file3.fq'."
  echo ""
  echo "Options (defaults in parentheses):"
  echo "  Scaffold Type:"
  echo "    --drna (default)"
  echo "        Scaffold reads are stranded forward relative to coding strand, and may contain U characters instead of Ts"
  echo "    --cdna-rev-stranded"
  echo "        Scaffold reads are stranded reverse complemented relative to coding strand"
  echo "    --cdna"
  echo "        Scaffold reads are NOT stranded"
  echo "    --sfq (default)"
  echo "        Scaffold reads are in FASTQ format"
  echo "    --sfa"
  echo "        Scaffold reads are in FASTA format"
  echo "  Illumina Type:"
  echo "    -u, --unstranded"
  echo "        Illumina reads are unstranded"
  echo "    -f, --fwd-stranded"
  echo "        Illumina reads are stranded s.t. the first mate originates from the RNA strand"
  echo "        Ignored if scaffold reads are not stranded"
  echo "    -r, --rev-stranded (default)"
  echo "        Illumina reads are stranded s.t. the first mate is the reverse complement of the RNA strand"
  echo "        Ignored if scaffold reads are not stranded"
  echo "    --ifq (default)"
  echo "        Illumina reads are in FASTQ format; Mutually exclusive with --ifa"
  echo "    --ifa"
  echo "        Illumina reads are in FASTA format; Mutually exclusive with --ifq"
  echo "  Consensus Collapsing:"
  echo "    -m, --score-matrix <path>"
  echo "        Provide an alternative scoring matrix to use in partial order alignment"
  echo "    -d, --isoform-delta (35)"
  echo "        Maximum indel size to be 'corrected', beyond this size a new isoform is declared. Must be between 0 and 255"
  echo "    -e, --ends-delta (35)"
  echo "        Maximum size at the ends of isoforms to 'correct' before splitting. Must be between 0 and 255"
  echo "    -i, --max-iterations (5)"
  echo "        Maximum number of iterations to align to and correct scaffolds. Does not include optional final polshing step"
  echo "        Note: Providing a value of 0 will not perform any graph based illumina correction"
  echo "    -w, --illumina-weight (10)"
  echo "        Weight of illumina reads relative to nanopore reads when generating consensus"
  echo "    --final-polish (default)"
  echo "        Include a final correction of individual isoforms, not in a splice graph"
  echo "    --no-final-polish"
  echo "        Do not do a final correction of individual isoforms not in a splice graph" 
  echo "  Ouput:"
  echo "    -o, --output-dir <path> (conduit/)"
  echo "        <path> where corrected clusters will be written"
  echo "        NOTE: THIS WILL OVERWRITE EXISTING FILES!"
  echo "    -n, --no-intermediates (default)"
  echo "        Does not save FASTA file generated for intermediate rounds of polishing"
  echo "    -s, --save-intermediates"
  echo "        Saves the FASTA file generated for intermediate rounds of polishing"
  echo "  Bowtie2:"
  echo "    --end-to-end (default)"
  echo "        Align Illumina reads to ONT scaffolds in end-to-end alignment mode"
  echo "    --local"
  echo "        Align Illumina reads to ONT scaffolds in local alignment mode"
  echo "  Miscellaneous:"
  echo "    -h, --help"
  echo "        Display this help message and exit"
  echo "    -v, --version"
  echo "        Display the installed version number of CONDUIT and exit"
  echo "    --tmp-dir <path> (conduit-tmp/)"
  echo "        <path> where temporary files will be created"
  echo "    -t, --threads (4)"
  echo "        Number of threads to run in parallel (used for both Bowtie2 and Partial Order Graph correction)"
  # echo "        NOTE: Providing a value of 0 will attempt to autodetect the number of CPUs availible and use that." #TODO
  # echo "              If CPU number cannot be detected, the default of 4 threads will be used. #TODO

proc returnFalse() : bool {.thread.} = 
  return false

proc convertFASTQtoFASTA(infilepath,outfilepath:string) = 
  var infile,outfile : File
  discard open(infile,infilepath,fmRead)
  discard open(outfile,outfilepath,fmWrite)
  while true:
    try:
      let line1 = infile.readLine()
      try:
        assert line1[0] == '@'
      except AssertionError:
        echo "File not in FASTQ format"
        raise
      outfile.write(&">{line1[1..^1]}")
      outfile.write(infile.readLine())
      discard infile.readLine()
      discard infile.readLine()
    except EOFError:
      break
  infile.close()
  outfile.close()

proc runPOAandCollapsePOGraph(intuple : (string,string,string,string,uint16,uint16)) {.thread.} =
  let (infilepath,outdir,matrix_filepath,format,isoform_delta,ends_delta) = intuple 
  let trim = infilepath.split(os.DirSep)[^1].split(".")[0]
  var fasta_file : string
  if format == "fasta":
    fasta_file = infilepath
  elif format == "fastq":
    fasta_file = &"{outdir}{trim}.tmp.fa"
    convertFASTQtoFASTA(infilepath,fasta_file)
  var seq_file : PFile = fopen(cstring(fasta_file), "r")
  # let matrix_filepath : cstring = "../poaV2/myNUC3.4.4.mat"
  var po = getPOGraphFromFasta(seq_file,cast[cstring](matrix_filepath),cint(1),matrix_scoring_function)
  if format == "fastq":
    removeFile(fasta_file)
  
  var trim_po = poParser.convertPOGraphtoTrimmedPOGraph(po)
  var representative_paths = poParser.getRepresentativePaths3(addr trim_po, psi = isoform_delta, ends_delta = ends_delta) #TODO - modify to accept ends_delta as well
  let consensus_po = poParser.buildConsensusPO(addr po, representative_paths,trim)
  let outFASTAfilepath = &"{outdir}fasta/{trim}.consensus.fa"
  var outFASTAfile : File
  discard open(outFASTAfile,outFASTAfilepath,fmWrite)
  poParser.writeCorrectedReads(consensus_po,outFASTAfile)
  outFASTAfile.close()

proc runGraphBasedIlluminaCorrection(intuple : (string,string,string,uint64,uint16,uint16)) : bool {.thread.} = #TODO convert from passing tuple back to passing vars individually; relic of older threading approach 
  let (tmp_dir, trim, matrix_filepath, iter,isoform_delta,ends_delta) = intuple
  let last_fasta_dir = &"{tmp_dir}{iter-1}/fasta/"
  let this_fasta_dir = &"{tmp_dir}{iter}/fasta/"

  let last_fasta_filepath : cstring = &"{last_fasta_dir}{trim}.consensus.fa"
  var seq_file : PFile = fopen(last_fasta_filepath, "r")
  var po = getPOGraphFromFasta(seq_file,matrix_filepath,cint(1),matrix_scoring_function)

  let this_fasta_filepath = &"{this_fasta_dir}{trim}.consensus.fa"

  let bamfilepath = &"{tmp_dir}{iter-1}/alignments.bam"

  var bam : Bam
  var trim_po = convertPOGraphtoTrimmedPOGraph(po)
  discard open(bam,bamfilepath,index=true)
  illuminaPolishPOGraph(addr trim_po, bam)
  var representative_paths = getRepresentativePaths3(addr trim_po, psi = isoform_delta,ends_delta = ends_delta)
  var records = getFastaRecordsFromTrimmedPOGraph(addr trim_po, representative_paths, trim)
  var outfile : File
  discard open(outfile,this_fasta_filepath,fmWrite)
  writeCorrectedReads(records,outfile)
  outfile.close()
  result = sameFileContent(cast[string](last_fasta_filepath),this_fasta_filepath)
  removeFile(cast[string](last_fasta_filepath))

proc runLinearBasedIlluminaCorrection(intuple : (string,string,uint64,uint64,uint16)) {.thread.} = 
  let (tmp_dir, trim, converged_iter, iter,isoform_delta) = intuple
  let last_fasta_dir = &"{tmp_dir}{converged_iter}/fasta/"
  let this_fasta_dir = &"{tmp_dir}{iter}/fasta/"

  let last_fasta_filepath = &"{last_fasta_dir}{trim}.consensus.fa"
  let this_fasta_filepath = &"{this_fasta_dir}{trim}.consensus.fa"

  let bamfilepath = &"{tmp_dir}{iter-1}/alignments.bam"

  var infile : File
  var bam : Bam
  discard open(infile,last_fasta_filepath)
  var reads = parseFasta(infile)
  var corrected : seq[FastaRecord]
  discard open(bam,bamfilepath,index=true)
  for read in reads:
    var trim_po = getTrimmedGraphFromFastaRecord(read)
    illuminaPolishPOGraph(addr trim_po, bam,debug=true)
    discard getRepresentativePaths3(addr trim_po, psi = isoform_delta)
    corrected.add(FastaRecord(read_id : read.read_id, sequence : getSequenceFromPath(trim_po,trim_po.reads[0].corrected_path)))
  var outfile : File
  discard open(outfile,this_fasta_filepath,fmWrite)
  writeCorrectedReads(corrected,outfile)
  outfile.close()

proc removeFiles(files : openArray[string]) =
  for file in files:
    removeFile(file)

proc createDirs(dirs : openArray[string]) = 
  for dir in dirs:
    createDir(dir)

proc combineFiles(indirectory : string, intrims : openArray[string], outfilepath : string) = 
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  for i,trim in intrims:
    let filepath = &"{indirectory}{trim}.consensus.fa"
    var file : File
    discard open(file,filepath,fmRead)
    outfile.write(file.readAll)
  outfile.close()

proc combineFilesIntermediate(indirectory : string, intrims : openArray[string], outfilepath : string, last_correction : Table[int,int]) = 
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  for i,trim in intrims:
    if i in last_correction:
      continue
    let filepath = &"{indirectory}{trim}.consensus.fa"
    var file : File
    discard open(file,filepath,fmRead)
    outfile.write(file.readAll)
  outfile.close()

proc combineFilesFinal(tmp_directory : string,last_num : uint64, intrims : openArray[string], outfilepath : string, last_correction : Table[int,int]) = 
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  for i,trim in intrims:
    var last = last_num
    if i in last_correction:
      last = uint64(last_correction[i]) #TODO convert last_correction to uint64 types.
    let indirectory = &"{tmp_directory}{last}/fasta/"
    let filepath = &"{indirectory}{trim}.consensus.fa"
    var file : File
    discard open(file,filepath,fmRead)
    outfile.write(file.readAll)
  outfile.close()

proc getBowtie2options(opt : ConduitOptions, index_prefix, sam : string) : seq[string] = 
  var arguments : seq[string]
  arguments.add("--xeq")
  arguments.add("--no-unal")
  arguments.add("-p")
  arguments.add(&"{opt.thread_num}")
  arguments.add(opt.bowtie_strand_constraint)
  arguments.add(opt.bowtie_alignment_mode)
  arguments.add("--n-ceil")
  arguments.add("L,0,0")
  arguments.add("-x")
  arguments.add(index_prefix)
  for arg in opt.bowtie_read_inputs.split(" "):
    arguments.add(arg)
  arguments.add("-S")
  arguments.add(sam)
  return arguments

proc parseOptions() : ConduitOptions = 
  var clusters_directory = ""
  var clusters_directory_flag = false

  var mate1s : seq[string]
  var mate2s : seq[string]
  var unpaireds : seq[string]
  var interleaved : seq[string]
  var bams : seq[string]

  var illumina_strand="reverse"
  var illumina_strand_flag = false

  var illumina_format="fastq"
  var illumina_format_flag = false

  var nanopore_type="drna"
  var nanopore_type_flag = false

  var nanopore_strand="forward"
  # var nanopore_strand_flag = false

  var nanopore_format="fastq"
  var nanopore_format_flag = false

  var score_matrix_path = ""
  var score_matrix_path_flag = false

  var isoform_delta = 35'u64
  var isoform_delta_flag = false

  var ends_delta = 35'u64
  var ends_delta_flag = false

  var max_iterations = 5'u64
  var max_iterations_flag = false

  var illumina_weight = 10'u64
  var illumina_weight_flag = false

  var final_polish = true
  var final_polish_flag = false

  var output_dir = "conduit-out/"
  var output_dir_flag = false

  var intermediates = false
  var intermediates_flag = false

  var local = false
  var local_flag = false

  var help_flag = true
  var version_flag = false

  var tmp_dir = "conduit-tmp/"
  var tmp_dir_flag = false

  var thread_num = 4'u64
  var thread_num_flag = false

  var run_flag = true


  var i = 0
  var mode = ""
  var last = ""

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
      of "nano", "hybrid":
        case kind:
          of cmdEnd:
            break
          of cmdShortOption, cmdLongOption:
            if last != "":
              echo &"ERROR - Option {last} provided without an argument"
              run_flag = false
              break
            case key:
              of "1":
                if val != "":
                  mate1s.add(val)
                else:
                  last = "1"
              of "2":
                if val != "":
                  mate2s.add(val)
                else:
                  last = "2"
              of "U":
                if val != "":
                  unpaireds.add(val)
                else:
                  last = "U"
              of "interleaved":
                if val != "":
                  interleaved.add(val)
                else:
                  last = "interleaved"
              of "b":
                if val != "":
                  bams.add(val)
                else:
                  last = "b"
              of "drna":
                if not nanopore_type_flag:
                  nanopore_type_flag = true
                else:
                  run_flag = false
                  echo "ERROR - Multiple scaffold types input, choose one of \"--drna\", \"--cdna-rev-stranded\", or \"--cdna\""
              of "cdna-rev-stranded":
                if not nanopore_type_flag:
                  nanopore_type_flag = true
                  nanopore_type = "cdna"
                  nanopore_strand = "reverse"
                else:
                  run_flag = false
                  echo "ERROR - Multiple scaffold types input, choose one of \"--drna\", \"--cdna-rev-stranded\", or \"--cdna\""
                  help_flag = true
                  break
              of "cdna":
                if not nanopore_type_flag:
                  nanopore_type_flag = true
                  nanopore_type = "cdna"
                  nanopore_strand = "unstranded"
                else:
                  run_flag = false
                  echo "ERROR - Multiple scaffold types input, choose one of \"--drna\", \"--cdna-rev-stranded\", or \"--cdna\""
                  help_flag = true
                  break
              of "sfq":
                if not nanopore_format_flag:
                  #already fastq
                  nanopore_format_flag = true
                else:
                  run_flag = false
                  echo "ERROR - Multiple scaffold format input, choose one of FASTA (--sfa) or FASTQ (--sfq)"
                  help_flag = true
                  break
              of "sfa":
                if not nanopore_format_flag:
                  #already fastq
                  nanopore_format_flag = true
                  nanopore_format = "fasta"
                else:
                  run_flag = false
                  echo "ERROR - Multiple scaffold format input, choose one of FASTA (--sfa) or FASTQ (--sfq)"
                  help_flag = true
                  break
              of "m", "score-matrix":
                if not score_matrix_path_flag:
                  score_matrix_path_flag = true
                  if val != "":
                    score_matrix_path = val
                  else:
                    last = "score-matrix"
                else:
                  run_flag = false
                  echo "ERROR - Multiple score matrices provided"
                  break
              of "u", "unstranded":
                if not illumina_strand_flag:
                  illumina_strand_flag = true
                  illumina_strand = "unstranded"
                else:
                  run_flag = false
                  echo "ERROR - Multiple illumina strand input, choose one of (-u,--unstranded), (-f,--fwd-stranded), (-r,--rev-stranded)"
                  help_flag = true
                  break
              of "f", "fwd-stranded":
                if not illumina_strand_flag:
                  illumina_strand_flag = true
                  illumina_strand = "forward"
                else:
                  run_flag = false
                  echo "ERROR - Multiple illumina strand input, choose one of (-u,--unstranded), (-f,--fwd-stranded), (-r,--rev-stranded)"
                  help_flag = true
                  break
                continue
              of "r", "rev-stranded":
                if not illumina_strand_flag:
                  illumina_strand_flag = true
                else:
                  run_flag = false
                  echo "ERROR - Multiple illumina strand input, choose one of (-u,--unstranded), (-f,--fwd-stranded), (-r,--rev-stranded)"
                  help_flag = true
                  break
              of "ifq":
                if not illumina_format_flag:
                  illumina_format_flag = true
                else:
                  run_flag = false
                  echo "ERROR - Multiple illumina format input, choose one of FASTA (--ifa), or FASTQ (--ifq)"
                  help_flag = true
                  break
              of "ifa":
                if not illumina_format_flag:
                  illumina_format_flag = true
                  illumina_format = "fasta"
                else:
                  run_flag = false
                  echo "ERROR - Multiple illumina format input, choose one of FASTA (--ifa), or FASTQ (--ifq)"
                  help_flag = true
                  break
              of "d", "isoform-delta":
                if not isoform_delta_flag:
                  isoform_delta_flag = true
                  if val != "":
                    isoform_delta = parseUInt(val)
                  else:
                    last = "isoform-delta"
                else:
                  run_flag = false
                  echo "ERROR - Isoform delta specified multiple times"
                  break
              of "e", "ends-delta":
                if not ends_delta_flag:
                  ends_delta_flag = true
                  if val != "":
                    ends_delta = parseUInt(val)
                  else:
                    last = "ends-delta"
                else:
                  run_flag = false
                  echo "ERROR - Ends delta specified multiple times"
                  break
              of "i", "max-iterations":
                if not max_iterations_flag:
                  max_iterations_flag = true
                  if val != "":
                    max_iterations = parseUInt(val)
                  else:
                    last = "max-iterations"
                else:
                  run_flag = false
                  echo "ERROR - Max iterations specified multiple times"
                  break
              of "w", "illumina-weight":
                if not illumina_weight_flag:
                  illumina_weight_flag = true
                  if val != "":
                    illumina_weight = parseUInt(val)
                  else:
                    last = "illumina-weight"
                else:
                  run_flag = false
                  echo "ERROR - Multiple Illumina weights specified"
                  break
              of "final-polish":
                if not final_polish_flag:
                  final_polish_flag = true
                elif not final_polish:
                  run_flag = false
                  echo "ERROR - Conflicting flags: --final-polish and --no-final-polish both specified"
                  break
              of "no-final-polish":
                if not final_polish_flag:
                  final_polish_flag = true
                  final_polish = false
                elif final_polish:
                  run_flag = false
                  echo "ERROR - Conflicting flags: --final-polish and --no-final-polish both specified"
                  break
              of "n", "no-intermediates":
                if not intermediates_flag:
                  intermediates_flag = true
                elif intermediates:
                  run_flag = false
                  echo "ERROR - Conflicting flags: --no-intermediates and --save-intermediates both specified"
                  break
              of "s", "save-intermediates":
                if not intermediates_flag:
                  intermediates_flag = true
                  intermediates = true
                elif not intermediates:
                  run_flag = false
                  echo "ERROR - Conflicting flags: --no-intermediates and --save-intermediates both specified"
                  break
              of "o", "output-dir":
                if not output_dir_flag:
                  output_dir_flag = true
                  if val != "":
                    output_dir = val
                  else:
                    last = "output-dir"
                else:
                  run_flag = false
                  echo "ERROR - Output directory specified multiple times"
                  break
              of "end-to-end":
                if not local_flag:
                  local_flag = true
                  local = false
                else:
                  if local:
                    run_flag = false
                    echo "ERROR - Conflicting flags: --local and --end-to-end both specified"
                    break
              of "local":
                if not local_flag:
                  local_flag = true
                  local = true
                else:
                  if not local:
                    run_flag = false
                    echo "ERROR - Conflicting flags: --local and --end-to-end both specified"
                    break
              of "t", "threads":
                if not thread_num_flag:
                  if val == "":
                    last = "threads"
                  else:
                    thread_num = parseUInt(val)
                else:
                  run_flag = false
                  echo "ERROR - Threads specified multiple times"
                  break
              of "tmp-dir":
                if not tmp_dir_flag:
                  tmp_dir_flag = true
                  if val == "":
                    last = "tmp-dir"
                  else:
                    tmp_dir = val
                else:
                  run_flag = false
                  echo "ERROR - Temporary directory specified twice"
                  break
              of "h", "help":
                help_flag = true
                run_flag = false
                break
              of "v", "version":
                version_flag = true
                run_flag = false
                break
              else:
                var dash = ""
                if kind == cmdShortOption:
                  dash = "-"
                else: #kind == cmdLongOption
                  dash = "--"
                echo &"ERROR - unknown option {dash}{key} provided"
          of cmdArgument:
            case last:
              of "1":
                mate1s.add(key)
              of "2":
                mate2s.add(key)
              of "U":
                unpaireds.add(key)
              of "interleaved":
                interleaved.add(key)
              of "b":
                bams.add(key)
              of "score-matrix":
                score_matrix_path = key
              of "isoform-delta":
                isoform_delta = parseUInt(key)
              of "ends-delta":
                ends_delta = parseUInt(key)
              of "max-iterations":
                max_iterations = parseUInt(key)
              of "illumina-weight":
                illumina_weight = parseUInt(key)
              of "output-dir":
                output_dir = key
              of "threads":
                thread_num = parseUInt(key)
              of "tmp-dir":
                if key[^1] == os.DirSep:
                  tmp_dir = key
                else:
                  tmp_dir = key & os.DirSep
              of "":
                # echo "Argument: ", key
                if not clusters_directory_flag:
                  clusters_directory_flag = true
                  clusters_directory = key
                else:
                  run_flag = false
                  echo "ERROR - Argument provided without associated option; please provide one <clusters_directory> only"
                  break
              else:
                echo &"ERROR - unknown option {last} provided"
                run_flag = false
                break
            last = ""
      else:
        help_flag = true
        break
  if clusters_directory == "":
    echo "ERROR - <clusters directory> must be specified"
    help_flag = true
    run_flag = false

  if help_flag:
    case mode:
      of "nano":
        writeNanoHelp()
      of "hybrid":
        writeHybridHelp()
      else:
        echo "ERROR - first argument must specify correction mode, \"nano\" or \"hybrid\""
        writeDefaultHelp()
  if version_flag:
    echo conduitVersion()
  var files : seq[string]
  var bowtie_strand_constraint : string
  var bowtie_alignment_mode : string
  var bowtie_read_inputs : string
  if run_flag and mode == "hybrid":
    # Determine strand relationship between nanopore and illumina reads:
    if (nanopore_strand == "forward" and illumina_strand == "reverse") or (nanopore_strand == "reverse" and illumina_strand == "forward"):
      bowtie_strand_constraint = "--nofw"
    elif (nanopore_strand == "forward" and illumina_strand == "forward") or (nanopore_strand == "reverse" and illumina_strand == "reverse"):
      bowtie_strand_constraint = "--norc"
    
    # Format illumina inputs in a manner readable by Bowtie2
    var bowtie_mate1_inputs = ""
    var bowtie_mate2_inputs = ""
    var bowtie_unpaired_inputs = ""
    var bowtie_interleaved_inputs = ""
    var bowtie_bam_inputs = ""
    if mate1s.len != 0:
      bowtie_mate1_inputs = "-1 " & mate1s.join(",") & " "
    if mate2s.len != 0:
      bowtie_mate2_inputs = "-2 " & mate2s.join(",") & " "
    if unpaireds.len != 0:
      bowtie_unpaired_inputs = "-U " & unpaireds.join(",") & " "
    if interleaved.len != 0:
      bowtie_interleaved_inputs = "--interleaved " & interleaved.join(",") & " "
    if bams.len != 0:
      bowtie_bam_inputs = "-b " & bams.join(",") & " "
    bowtie_read_inputs = &"{bowtie_mate1_inputs}{bowtie_mate2_inputs}{bowtie_unpaired_inputs}{bowtie_interleaved_inputs}{bowtie_bam_inputs}"
    if bowtie_read_inputs == "":
      echo "ERROR - No Illumina data provided"
      run_flag = false
    
    if local:
      bowtie_alignment_mode = "--very-sensitive-local"
    else:
      bowtie_alignment_mode = "--very-sensitive"
    
    if nanopore_format == "fasta":
      for file in walkFiles(&"{clusters_directory}*.fa"):
        files.add(file)
      for file in walkFiles(&"{clusters_directory}*.fasta"):
        files.add(file)
      if files.len == 0:
        run_flag = false
        echo &"ERROR - No files of type .fa or .fasta found in <clusters directory> {clusters_directory}"
        echo "NOTE: We don't currently support .gzip'd or bzip2'd scaffold files, though support for these formats is coming"
    elif nanopore_format == "fastq":
      for file in walkFiles(&"{clusters_directory}*.fq"):
        files.add(file)
      for file in walkFiles(&"{clusters_directory}*.fastq"):
        files.add(file)
      if files.len == 0:
        run_flag = false
        echo &"ERROR - No files of type .fq or .fastq found in <clusters directory> {clusters_directory}"
        echo "NOTE: We don't currently support .gzip'd or bzip2'd scaffold files, though support for these formats is coming"
    if isoform_delta > 255'u64 or isoform_delta < 0'u64:
      echo "ERROR - Isoform delta must be between 0 and 255"
      run_flag = false
    elif isoform_delta < 15'u64:
      echo "WARNING - Low isoform delta values will increase the number of distinct isoforms and dramatically increase runtime"
      echo "          An isoform delta value of 15 or above is reccomended"
    if ends_delta > 255'u64 or ends_delta < 0'u64:
      echo "ERROR - Ends delta must be between 0 and 255"
      run_flag = false
    if thread_num < 1'u64:
      echo "ERROR - Must run at least 1 thread"
      run_flag = false
    if illumina_weight == 0'u64:
      echo "ERROR - Illumina weight must be greater than 0"
      run_flag = false
    elif illumina_weight < 5'u64:
      echo "WARNING - We reccomend weighing Illumina reads by at least 5x their long-read counterparts"
  var trims = newSeq[string](files.len)
  for i,infilepath in files:
    trims[i] = infilepath.split(os.DirSep)[^1].split(".")[0]
  return ConduitOptions(run_flag : run_flag,
                        final_polish : final_polish,
                        intermediates : intermediates,
                        mode : mode,
                        nanopore_format : nanopore_format,
                        bowtie_strand_constraint : bowtie_strand_constraint,
                        bowtie_alignment_mode : bowtie_alignment_mode,
                        bowtie_read_inputs : bowtie_read_inputs,
                        score_matrix_path : score_matrix_path,
                        output_dir : output_dir,
                        tmp_dir : tmp_dir,
                        files : files,
                        trims : trims,
                        isoform_delta : isoform_delta,
                        ends_delta : ends_delta,
                        max_iterations : max_iterations,
                        illumina_weight : illumina_weight,
                        thread_num : thread_num ) 

proc main() = #TODO - Add output indicating completion percentage for each iteration.
  let opt = parseOptions()
  if opt.mode == "hybrid" and opt.run_flag:
    if not existsDir(opt.output_dir):
      createDir(opt.output_dir)

    var tmp_already_existed = true
    if not existsDir(opt.tmp_dir):
      tmp_already_existed = false
      createDir(opt.tmp_dir)

    let directory_number = opt.max_iterations + uint64(opt.final_polish)
    for i in 0..directory_number:
      createDirs([&"{opt.tmp_dir}{i}/", &"{opt.tmp_dir}{i}/fasta/"]) #TODO get rid of fasta subdirectory, no longer needed.

    let p = tps.newThreadPool(int(opt.thread_num))
    for file in opt.files:
      p.spawn runPOAandCollapsePOGraph((file, &"{opt.tmp_dir}0/", opt.score_matrix_path, opt.nanopore_format, uint16(opt.isoform_delta), uint16(opt.ends_delta)))
    p.sync()

    var last_correction : Table[int,int]
    for iter in 1..opt.max_iterations:
      let last_dir = &"{opt.tmp_dir}{iter-1}/"
      let last_fasta_dir = &"{last_dir}fasta/"

      # let cur_dir = &"{opt.tmp_dir}{iter}/"

      var last_consensus : string
      if opt.intermediates:
        last_consensus = &"{opt.output_dir}conduit_consensuses_iter{iter-1}.fa"
      else:
        last_consensus = &"{opt.tmp_dir}conduit_consensuses_iter{iter-1}.fa"

      removeFile(last_consensus)
      combineFilesIntermediate(last_fasta_dir,opt.trims,last_consensus,last_correction)

      let index_prefix = &"{last_dir}bowtie2_index"
      discard execProcess("bowtie2-build", args =["--threads",&"{opt.thread_num}", last_consensus, index_prefix],options={poUsePath})
      
      if not opt.intermediates:
        removeFile(last_consensus)
      
      let sam = &"{last_dir}alignments.sam"
      let arguments = getBowtie2options(opt,index_prefix,sam)
      discard execProcess("bowtie2", args = arguments, options={poUsePath,poStdErrToStdOut})
      
      let bam = &"{last_dir}alignments.bam"
      discard execProcess(&"samtools sort -@ {opt.thread_num} {sam} > {bam}", options={poEvalCommand,poUsePath,poStdErrToStdOut})
      discard execProcess("samtools", args=["index", bam],options={poUsePath})
      
      removeFiles([sam, &"{index_prefix}.1.bt2", &"{index_prefix}.2.bt2", &"{index_prefix}.3.bt2", &"{index_prefix}.4.bt2", &"{index_prefix}.rev.1.bt2", &"{index_prefix}.rev.2.bt2"])

      let p = tps.newThreadPool(int(opt.thread_num))
      var converged = newSeq[tps.FlowVar[bool]](opt.trims.len)
      for i,trim in opt.trims:
        if i in last_correction:
          converged[i] = p.spawn returnFalse()
          continue
        converged[i] = p.spawn runGraphBasedIlluminaCorrection((opt.tmp_dir,trim,opt.score_matrix_path,iter,uint16(opt.isoform_delta),uint16(opt.ends_delta)))
      p.sync()

      for i,converge in converged:
        if converge.read():
          last_correction[i] = int(iter)

      removeFiles([bam, &"{bam}.bai"])
    if opt.final_polish:
      let iter = opt.max_iterations + 1
      let last_dir = &"{opt.tmp_dir}{iter-1}/"
      let last_fasta_dir = &"{last_dir}fasta/"
      
      var last_consensus = ""
      if opt.intermediates:
        last_consensus = &"{opt.output_dir}conduit_consensuses_iter{iter-1}.fa"
      else:
        last_consensus = &"{opt.tmp_dir}conduit_consensuses_iter{iter-1}.fa"
      removeFile(last_consensus)
      combineFilesFinal(opt.tmp_dir, iter-1, opt.trims, last_consensus, last_correction)

      let index_prefix = &"{last_dir}bowtie2_index"
      discard execProcess("bowtie2-build", args =["--threads",&"{opt.thread_num}", last_consensus, index_prefix],options={poUsePath})
      
      if not opt.intermediates:
        removeFile(last_consensus)

      let sam = &"{last_dir}alignments.sam"
      let arguments = getBowtie2options(opt,index_prefix,sam)
      discard execProcess("bowtie2", args = arguments, options={poUsePath,poStdErrToStdOut})
      
      let bam = &"{last_dir}alignments.bam"
      discard execProcess(&"samtools sort -@ {opt.thread_num} {sam} > {bam}", options={poEvalCommand,poUsePath})
      discard execProcess("samtools", args=["index", bam],options={poUsePath})
      
      removeFiles([sam, &"{index_prefix}.1.bt2", &"{index_prefix}.2.bt2", &"{index_prefix}.3.bt2", &"{index_prefix}.4.bt2", &"{index_prefix}.rev.1.bt2", &"{index_prefix}.rev.2.bt2"])

      let p = tps.newThreadPool(int(opt.thread_num))
      for i,trim in opt.trims:
        var converged_iter = iter - 1
        if i in last_correction:
          converged_iter = uint64(last_correction[i])
        p.spawn runLinearBasedIlluminaCorrection((opt.tmp_dir,trim,converged_iter,iter,uint16(opt.isoform_delta)))
      p.sync()

      removeFiles([bam, &"{bam}.bai"])
      removeDir(last_fasta_dir)
    
    let final_consensus_path = &"{opt.output_dir}conduit_final_consensuses.fa"
    if existsFile(final_consensus_path):
      removeFile(final_consensus_path)
    var last_fasta_dir = &"{opt.tmp_dir}{directory_number}/fasta/"

    combineFiles(last_fasta_dir, opt.trims, final_consensus_path)
    
    if not tmp_already_existed:
      removeDir(opt.tmp_dir)

main()