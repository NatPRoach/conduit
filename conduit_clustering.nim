import os
import times
import osproc
import parseopt
import strutils
import strformat
import tables
import heapqueue
import threadpool_simple as tps
import hts
import sets

type
  ClusteringOptions = object
    file : string
    mode : char
    cluster_zero_introns : bool
    output_dir : string
    prefix : string
    out_type : string


proc conduitClusterVersion() : string =
  return "CONDUIT Clustering Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu, https://github.com/NatPRoach/conduit/ )"

proc writeClusterHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitClusterVersion()
  echo "Usage:"
  echo "  ./conduit_cluster [options] input.bam"
  echo ""
  echo "Options (defaults in parentheses):"
  echo "  Clustering mode:"
#TODO - Figure out if there's any computationally reasonable way to do ss tolerance?
#  echo "    -s, --splice-site-tolerance (15)"
#  echo "        Allow tolerance"
  echo "    --ss (default)"
  echo "        Cluster reads with at least one splice site in common"
  echo "    --overlap"
  echo "        Cluster reads with at least one base of overlap in alignment"
  echo "    -z (default)"
  echo "        Cluster reads with zero introns by overlap"
  echo "    -d"
  echo "        Do not cluster reads with zero introns by overlap"
  echo "  Output:"
  echo "    -o, --output-dir <path> (clusters/)"
  echo "    -p, --prefix (cluster_)"
  echo "        The prefix for output clusters"
  echo "    --fq (default)"
  echo "        Output reads in FASTQ format"
  echo "    --fa"
  echo "        Output reads in FASTA format"

# proc outputTiming(outfilepath : string,time_seq : seq[Time],opts : ConduitOptions) =
#   var outfile : File
#   discard open(outfile,outfilepath,fmWrite)
#   outfile.close()

# proc outputSettings(outfilepath : string,opts : ConduitOptions) = 
#   var outfile : File
#   discard open(outfile,outfilepath,fmWrite)
#   outfile.close()

proc spliceSitesFromBamRecord( record : Record, stranded : bool = false) : (char, (uint64, uint64), seq[(uint64, uint64)]) = 
  if record.flag.unmapped:
      return
  let alignment_start = uint64(record.start)
  var ref_position = uint64(alignment_start)
  let cigar = record.cigar
  let flag = record.flag
  var strand : char
  var splice_sites : seq[(uint64, uint64)]
  if  ( not stranded ) or flag.has_flag(16'u16): # If we're not stranded treat everything as if it's (-) strand, else check for the actual strand
    strand = '-'
  else:
    strand = '+'
  for op in cigar:
    case int(op.op):
      of 0, 2, 7, 8:
        ref_position += uint64(op.len)
      of 3:
        let start = ref_position
        ref_position += uint64(op.len)
        splice_sites.add((start,ref_position))
      else:
        discard
  let alignment_end = ref_position
  result = (strand,(alignment_start,alignment_end),splice_sites)

proc spliceSitesTableFromBam( bam : Bam,  chromosome : string, stranded : bool = false) : HashSet[(char, seq[(uint64, uint64)])] = 
  for record in bam.query(chromosome):
    echo spliceSitesFromBamRecord(record, stranded)

proc parseOptions() : ClusteringOptions = 
  var file : string
  var file_flag = false

  var output_format = "fastq"
  var output_format_flag = false

  var cluster_mode = 's'
  var cluster_mode_flag = false

  var single_exon = true
  var single_exon_flag = false

  var output_dir = &"clusters{os.DirSep}"
  var output_dir_flag = false

  var prefix = &"cluster_"
  var prefix_flag = false

  var run_flag = true
  var help_flag = false
  var version_flag = false

  var last = ""

  for kind, key, val in getopt():
    case kind:
      of cmdEnd:
        break
      of cmdShortOption, cmdLongOption:
        if last != "":
          echo &"ERROR - Option {last} provided without an argument"
          run_flag = false
          break
        
        case key:
          of "ss":
            if cluster_mode_flag and ( cluster_mode != 's' ):
              echo "ERROR - Conflicting flags --ss and --overlap"
              help_flag = true
              break
            else:
              cluster_mode_flag = true
              cluster_mode = 's'
          of "overlap":
            if cluster_mode_flag and ( cluster_mode != 'o' ):
              echo "ERROR - Conflicting flags --ss and --overlap"
              help_flag = true
              break
            else:
              cluster_mode_flag = true
              cluster_mode = 'o'
          of "z":
            if single_exon_flag and not single_exon:
              echo "ERROR - Conflicting flags -z and -d"
              help_flag = true
              break
            else:
              single_exon_flag = true
              single_exon = true
          of "d":
            if single_exon_flag and single_exon:
              echo "ERROR - Conflicting flags -z and -d"
              help_flag = true
              break
            else:
              single_exon_flag = true
              single_exon = false
          of "fq":
            if output_format_flag and output_format != "fastq":
              echo "ERROR - Conflicting flags --fq and --fa"
              help_flag = true
              break
            else:
              output_format_flag = true
              output_format = "fastq"
          of "fa":
            if output_format_flag and output_format != "fasta":
              echo "ERROR - Conflicting flags --fq and --fa"
              help_flag = true
              break
            else:
              output_format_flag = true
              output_format = "fasta"
          of "o", "output-dir":
            if not output_dir_flag:
              output_dir_flag = true
              if val != "":
                output_dir = val
              else:
                last = "output"
            else:
              echo "ERROR - two output dirs provided"
              help_flag = true
              break
          of "p", "prefix":
            if not prefix_flag:
              prefix_flag = true
              if val != "":
                output_dir = val
              else:
                last = "output"
            else:
              echo "ERROR - two output dirs provided"
              help_flag = true
              break
          of "h", "help":
            help_flag = true
            break
          of "v", "version":
            help_flag = false
            version_flag = true
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
          of "output":
            output_dir = key
          of "prefix":
            prefix = key
          of "":
            if not file_flag:
              file = key
            else:
              file_flag = true
              echo "ERROR - multiple files provided"
              break
          else:
            echo &"ERROR - unknown option {last} provided"
            run_flag = false
            break
        last = ""

  if help_flag:
    writeClusterHelp()
  if version_flag:
    echo conduitClusterVersion()

  return ClusteringOptions(file : file,
                           mode : cluster_mode,
                           cluster_zero_introns : single_exon,
                           output_dir : output_dir,
                           prefix : prefix,
                           out_type : output_format)
proc main() = 
  let opt = parseOptions()
  if not existsDir(opt.output_dir):
    createDir(opt.output_dir)
  var bam : Bam
  discard open(bam,opt.file,index=true)
  discard spliceSitesTableFromBam(bam, "chrI", true)

main()