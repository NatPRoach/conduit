import hts

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

proc spliceSitesFromBamRecord( record : Record, stranded : bool = false) : (string, char, (uint64, uint64), seq[(uint64, uint64)]) = 
  if record.flag.unmapped:
      return
  var chrom = record.chrom
  let alignment_start = uint64(record.start)
  var ref_position = uint64(alignment_start)
  let cigar = record.cigar
  let flag = record.flag
  var strand : char
  var splice_sites : seq[(uint64, uint64)]
  if flag.has_flag(16'u16):
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

  result = (chrom,strand,(alignment_start,alignment_end),splice_sites)

proc spliceSitesTableFromBam( bam:Bam, stranded : bool = false) = 
  for record in bam.items():
    discard spliceSitesFromBamRecord(record, stranded)