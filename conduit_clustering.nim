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
#  echo "    -s, --splice-site-tolerance (15)"
#  echo "        Allow tolerance"
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