import os
import osproc
import parseopt
import strutils
import strformat
import poParser

proc conduitVersion() : string =
  return "CONDUIT Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu, https://github.com/NatPRoach/conduit/ )"

proc writeDefaultHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitVersion()
  echo "Usage:"
  echo "\t./conduit <nano | hybrid>"

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
  #echo "    --cdna" #Not supported yet

proc writeHybridHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitVersion()
  echo "Usage:"
  echo "  ./conduit hybrid [options] <clusters_directory> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>}"
  echo "  <clusters_directory>   Directory containing the .fasta/.fa or .fastq/.fq files of reads separated by gene cluster"
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
  echo "    -d, --isoform-delta (35)"
  echo "        Maximum indel size to be 'corrected', beyond this size a new isoform is declared. Must be between 0 and 255"
  echo "    -e, --ends-delta (35)"
  echo "        Maximum size at the ends of isoforms to 'correct' before splitting. Must be between 0 and 255"
  echo "    -i, --max-iterations (5)"
  echo "        Maximum number of iterations to align to and correct scaffolds. Does not include optional final polshing step"
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


# case paramStr(1):
#   of "nano":
#     writeNanoHelp()
#   of "hybrid":
#     writeHybridHelp()
#   else:
#     echo "ERROR - first argument must specify correction mode, \"nano\" or \"hybrid\""
#     writeDefaultHelp()

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

var isoform_delta = 35'u64
var isoform_delta_flag = false

var ends_delta = 35'u64
var ends_delta_flag = false

var max_iterations = 5'u64
var max_iterations_flag = false

var final_polish = true
var final_polish_flag = false

var output_dir = "conduit/"
var output_dir_flag = false

var intermediates = false
var intermediates_flag = false

var local = false
var local_flag = false

var help_flag = true
var version_flag = false

var tmp_dir = "conduit-tmp/"
var tmp_dir_flag = false

var threads = 4'u64
var threads_flag = false

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
              continue
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
              if not threads_flag:
                if val == "":
                  last = "threads"
                else:
                  threads = parseUInt(val)
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
            of "isoform-delta":
              isoform_delta = parseUInt(key)
            of "ends-delta":
              ends_delta = parseUInt(key)
            of "max-iterations":
              max_iterations = parseUInt(key)
            of "output-dir":
              output_dir = key
            of "threads":
              threads = parseUInt(key)
            of "tmp-dir":
              tmp_dir = key
            of "":
              # echo "Argument: ", key
              if not clusters_directory_flag:
                clusters_directory_flag = true
                clusters_directory = key
              else:
                run_flag = false
                echo "ERROR - Argument provided without associated option; please provide one <clusters_directory> only"
            else:
              echo &"ERROR - unknown option {last} provided"
          last = ""
    else:
      help_flag = true
      break
if clusters_directory == "":
  echo "ERROR - <clusters directory> must be specified"
  help_flag = true
  run_flag = false

if run_flag and mode == "hybrid":
  # Determine strand relationship between nanopore and illumina reads:
  var bowtie_strand_constraint = ""
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
  let bowtie_read_inputs = &"{bowtie_mate1_inputs}{bowtie_mate2_inputs}{bowtie_unpaired_inputs}{bowtie_interleaved_inputs}{bowtie_bam_inputs}"
  if bowtie_read_inputs == "":
    echo "ERROR - No Illumina data provided"
    run_flag = false

  var files : seq[string]
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
  if threads < 1'u64:
    echo "ERROR - Must run at least 1 thread"
    run_flag = false
  if run_flag:
    #TODO - write code to more seamlessly integrate POAv2 or SIMD POA from Eyras lab with C bindings. Until then - use explicit path to POAv2
    # echo execShellCmd("../poaV2/poa")
    var tmp_already_existed = true
    if not existsDir(tmp_dir):
      tmp_already_existed = false
      discard execProcess("mkdir", args = [tmp_dir],options={poUsePath})
    for file in files:
      var file2 : string
      if nanopore_format == "fasta":
        file2 = file
      elif nanopore_format == "fastq":
        discard execProcess(&"seqtk seq {file} > {tmp_dir}tmp.fa",options={poEvalCommand,poUsePath})
        file2 = &"{tmp_dir}tmp.fa"
        #POA requires FASTA format - convert
      let out_po = "tmp.po"
      discard execProcess("../poaV2/poa", args =["-do_global", "-read_fasta", file2, "-po", out_po, "../poaV2/myNUC3.4.4.mat"],options={})
      if nanopore_format == "fastq":
        discard execProcess("rm",args=[&"{tmp_dir}tmp.fa"],options = {})
    if not tmp_already_existed:
      discard execProcess("rmdir", args = [tmp_dir],options={poUsePath})
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
# echo execShellCmd("../poaV2/poa -do_global -read_fasta clusters/fasta/cluster_1423.fa -po tmp.po ../poaV2/myNUC3.4.4.mat")
# echo execProcess("echo", args = ["$PATH"], options={poEchoCmd,poUsePath})
echo execProcess("seqtk",args = ["seq", "-A", "clusters/fasta/cluster_1423.fa", ">", "tmp.fa"],options={poEchoCmd,poUsePath})