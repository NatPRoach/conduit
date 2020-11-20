import os
import bitops
import times
import osproc
import parseopt
import strutils
import strformat
import tables
import threadpool_simple as tps
import hts
import sets
import genomeKDE
import algorithm
import poGraphUtils
import fasta
import fastq


type
  ClusteringOptions = object
    file : string
    mode : char
    cluster_zero_introns : bool
    output_dir : string
    prefix : string
    out_type : string
  
  SpliceSiteGraph = object
    adjacencies : seq[seq[uint32]]
    ss_to_index : Table[uint32, uint32]
    index_to_index : OrderedTable[uint32,uint32]
    previously_visited : seq[uint32]

# TODO -   Clean up line wrapping where necessary (make it readable)
# TODO - but also conform to nim style guide length restrictions
# TODO - where possible

# TODO -   Possibly change uint32 to uint64 where necessary.
# TODO - Reevaluate before publication.
proc conduitClusterVersion() : string =
  return "CONDUIT Clustering Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu," &
    " https://github.com/NatPRoach/conduit/ )"


proc writeClusterHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitClusterVersion()
  echo "Usage:"
  echo "  ./conduit_cluster [options] input.bam"
  echo ""
  echo "Options (defaults in parentheses):"
  echo "  Correction:"
  echo "    -r, --reference path/to/reference"
  echo "        Used to providing a FASTA reference with corresponding .fai" &
    " index"
  echo "        If provided reads will be 'corrected' to sequence corresponding"
  echo "        to their genomic alignment before output into cluster files."
  echo "        Forces outputing reads in FASTA format (--fa) as the original"
  echo "        quality scores become meaningless after this correction"
  echo ""
  echo "  Clustering mode:"
# TODO -
#[ Figure out if there's any computationally reasonable way to do ss tolerance?
   Reduce to unique ss locations (weighted by depth?)
   Reduce ss to corresponding peaks of splice sites to
     group splice sites together easier.
   Calculate peaks at sections of the genome w/ # of ss > some threshold
   Find local maxima, correct splice sites to discovered local maxima
   Optionally allow for reference genome based correction
     of individual reads
   Optionally allow for reference genome annotation based splice site
     definition for collapsing ]#
#  echo "    -s, --splice-site-tolerance (15)"
#  echo "        Allow tolerance in splice site collapsing step"
  echo "    --ss (default, only option at the moment)"
  echo "        Cluster reads with at least one splice site in common"
  # echo "    --overlap"
  # echo "        Cluster reads with at least one base of overlap in alignment"
  echo "    -z (default)"
  echo "        Cluster reads with zero introns by overlap"
  echo "    -d"
  echo "        Do not cluster reads with zero introns by overlap"
  echo ""
  echo "  Output:"
  echo "    -o, --output-dir <path> (clusters/)"
  echo "    -p, --prefix (cluster_)"
  echo "        The prefix for output clusters"
  echo "    --fq (default)"
  echo "        Output reads in FASTQ format"
  echo "    --fa"
  echo "        Output reads in FASTA format"

proc summaryFromBamRecord( record : Record,
                           stranded : bool = false)
                          :
                          (char,
                           (uint64, uint64),
                           seq[(uint64, uint64)]) = 
  if record.flag.unmapped:
      return
  let alignment_start = uint64(record.start)
  var ref_position = uint64(alignment_start)
  let cigar = record.cigar
  let flag = record.flag
  var strand : char
  var splice_sites : seq[(uint64, uint64)]
  # If we're not stranded treat everything as if it's (-) strand
  # else check for the actual strand
  if  ( not stranded ) or flag.has_flag(16'u16):
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

# TODO -
#[Figure out if bitshifting flag for strand makes sense.
At u32 assumes no chr larger than 2147483648 bp.
# Probably safe, just takes time to change the logic. ]#
proc collapsedSummariesFromBam( bam : Bam,
                                chromosome : string,
                                stranded : bool = false) :
                                (Table[char,
                                       OrderedTable[seq[(uint64, uint64)],
                                                    seq[string]]],
                                 Table[char,
                                       OrderedTable[(uint64, uint64),
                                                    seq[string]]]) = 
  for record in bam.query(chromosome):
    # TODO -
    #[ Figure out if there's anything to be done w/ secondary / 
       supplmentary alignments - there probably is w/r/t clustering of
       overlapping splice sites especially but it may take some clever
       approaches.
    ]#
    if record.flag.supplementary or record.flag.secondary:
      echo &"WARNING - At the moment secondary / supplementary alignments " &
        "are not supported"
      continue
    let summary = summaryFromBamRecord(record, stranded)
    if summary[2].len == 0:
      if summary[0] in  result[1]:
        if summary[1] in result[1][summary[0]]:
          result[1][summary[0]][summary[1]].add(record.qname)
        else:
          result[1][summary[0]][summary[1]] = @[record.qname]
      else:
        result[1][summary[0]] = [(summary[1], @[record.qname])].toOrderedTable
    else:
      if summary[0] in  result[0]:
        if summary[2] in result[0][summary[0]]:
          result[0][summary[0]][summary[2]].add(record.qname)
        else:
          result[0][summary[0]][summary[2]] = @[record.qname]
      else:
        result[0][summary[0]] = [(summary[2], @[record.qname])].toOrderedTable


proc cmpSpliced(a,b : (seq[(uint64,uint64)], seq[string])) : int =
  for i in 0..<min(a[0].len,b[0].len):
    result = system.cmp(a[0][i],b[0][i])
    if result != 0:
      return
  result = system.cmp(a[0].len,b[0].len)


proc cmpNonSpliced(a,b : ((char,(uint64,uint64)), seq[string])) : int =
  result = system.cmp(a[0],b[0])


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


proc bfs(ssgraph : ptr SpliceSiteGraph,
         cluster_number : uint32,
         starting_node : uint32) : seq[uint32] = 
  var to_visit : seq[uint32]
  ssgraph[].previously_visited[starting_node] = cluster_number
  to_visit.add(starting_node)
  while to_visit.len != 0:
    let u = to_visit.pop()
    result.add(u)
    for adj in ssgraph[].adjacencies[u]:
      if ssgraph[].previously_visited[adj] == 0:
        ssgraph[].previously_visited[adj] = cluster_number
        to_visit.add(adj)


proc findIsoformClusters(ssgraph : ptr SpliceSiteGraph) : seq[seq[uint32]] =
  var cluster_counter = 1'u32
  for isoform_adjacency_index in ssgraph[].index_to_index.keys:
    if ssgraph[].previously_visited[isoform_adjacency_index] == 0:
      var clustered_isoforms : seq[uint32]
      let clustered_nodes = bfs(ssgraph,
                                cluster_counter,
                                isoform_adjacency_index)
      for node_id in clustered_nodes:
        if node_id in ssgraph[].index_to_index:
          clustered_isoforms.add(ssgraph[].index_to_index[node_id])
      cluster_counter += 1
      result.add(clustered_isoforms.sorted)


proc getWeightedSpliceJunctionLocations(
        sstable : ptr OrderedTable[seq[(uint64, uint64)],seq[string]]) :
        (seq[(uint32,uint32)],seq[(uint32,uint32)]) = 
  var donor_table : OrderedTable[uint32,uint32]
  var acceptor_table : OrderedTable[uint32,uint32]
  for sss, read_ids in sstable[].pairs:
    let weight = uint32(read_ids.len)
    for (donor,acceptor) in sss:
      if uint32(donor) in donor_table:
        donor_table[uint32(donor)] += weight
      else:
        donor_table[uint32(donor)] = weight
      if uint32(acceptor) in acceptor_table:
        acceptor_table[uint32(acceptor)] += weight
      else:
        acceptor_table[uint32(acceptor)] = weight
  for ss, weight in donor_table.pairs:
    result[0].add((ss,weight))
  for ss, weight in acceptor_table.pairs:
    result[1].add((ss,weight))
  echo result

proc writeFASTAsFromBAM(bam : Bam,
                        read_id_to_cluster : ptr Table[string,int],
                        cluster_sizes : ptr CountTable[int],
                        query : string,
                        cluster_prefix : string,
                        starting_count : int) : int =
  var open_files : Table[int, File]
  var written_reads : CountTable[int]
  for record in bam.query(query):
    if record.flag.supplementary or record.flag.secondary:
      echo &"WARNING - At the moment secondary / supplementary alignments " &
        "are not supported"
      continue
    if record.qname notin read_id_to_cluster[]:
      continue
    let cluster_id = read_id_to_cluster[][record.qname]
    if cluster_id notin open_files:
      var file : File
      discard open(file,
                   &"{cluster_prefix}{cluster_id + starting_count}.fa",
                   fmWrite)
      result += 1
      # echo &"Opening cluster file {cluster_id + starting_count}"
      open_files[cluster_id] = file
      file.writeFASTArecordToFile(record)
    else:
      # echo &"Appending to cluster file {cluster_id + starting_count}"
      open_files[cluster_id].writeFASTArecordToFile(record)
    written_reads.inc(cluster_id)
    echo written_reads[cluster_id], "\t", cluster_sizes[][cluster_id]
    if written_reads[cluster_id] == cluster_sizes[][cluster_id]:
      # echo &"Closing cluster file {cluster_id + starting_count}"
      open_files[cluster_id].close
      open_files.del(cluster_id)
  # Just in case of secondary / supplementary alignments
  for cluster_id, open_file in open_files.pairs:
    open_file.close


proc correctBamRecordWithGenome(record : Record, fai : Fai) : FastaRecord = 
  let summary = summaryFromBamRecord(record, true)
  var nt_sequence : string
  var start_base = summary[1][0]
  var end_base = -1
  for donor,acceptor in summary[2]:
    end_base = donor
    nt_sequence.add(fai.get(record.chrom,start_base,end_base))
    start_base = acceptor
  end_base = summary[1][1]
  nt_sequence.add(fai.get(record.chrom,start_base,end_base))
  if summary[0] == '-':
    nt_sequence = nt_sequence.revComp
  result = FastaRecord( read_id : record.qname, sequence : nt_sequence)

proc writeFASTAsFromBAM(bam : Bam,
                        read_id_to_cluster : ptr Table[string,int],
                        cluster_sizes : ptr CountTable[int],
                        query : string,
                        cluster_prefix : string,
                        starting_count : int,
                        fai : Fai) : int =
  var open_files : Table[int, File]
  var written_reads : CountTable[int]
  for record in bam.query(query):
    if record.flag.supplementary or record.flag.secondary:
      echo "WARNING - At the moment secondary / supplementary alignments " &
        "are not supported"
      continue
    if record.flag.unmapped:
      echo "WARNING - Unmapped read detected, skipping"
      continue
    let fasta_record = correctBamRecordWithGenome(record, fai)
    if fasta_record.read_id notin read_id_to_cluster[]:
      continue
    let cluster_id = read_id_to_cluster[][fasta_record.read_id]
    if cluster_id notin open_files:
      var file : File
      discard open(file,&"{cluster_prefix}{cluster_id + starting_count}.fa",
        fmWrite)
      result += 1
      # echo &"Opening cluster file {cluster_id + starting_count}"
      open_files[cluster_id] = file
      file.writeFASTArecordToFile(fasta_record)
    else:
      # echo &"Appending to cluster file {cluster_id + starting_count}"
      open_files[cluster_id].writeFASTArecordToFile(fasta_record)
    written_reads.inc(cluster_id)
    echo written_reads[cluster_id], "\t", cluster_sizes[][cluster_id]
    if written_reads[cluster_id] == cluster_sizes[][cluster_id]:
      # echo &"Closing cluster file {cluster_id + starting_count}"
      open_files[cluster_id].close
      open_files.del(cluster_id)
  # Just in case of secondary / supplementary alignments
  for cluster_id, open_file in open_files.pairs:
    open_file.close


proc main() = 
  let opt = parseOptions()
  if not dirExists(opt.output_dir):
    createDir(opt.output_dir)
  var bam : Bam
  discard open(bam,opt.file,index=true)
  var (strand_spliced_table, strand_nonspliced_table) =
    collapsedSummariesFromBam(bam, "chrI", true)
  bam.close
  var output_file_count = 0
  for strand, spliced_table in strand_spliced_table.mpairs:
    spliced_table.sort(cmpSpliced)
    ### Populate SpliceSiteGraph
    var ssgraph : SpliceSiteGraph
    var isoform_node_index = 0'u32
    # isoform_node_index.setBit(31'u8)
    for sss in spliced_table.keys:
      var isoform_adjacency_index = uint32(ssgraph.adjacencies.len)
      ssgraph.adjacencies.add(@[])
      ssgraph.previously_visited.add(0)
      ssgraph.index_to_index[isoform_adjacency_index] = isoform_node_index
      isoform_node_index += 1
      for (donor,acceptor) in sss:
        # echo donor, "\t", acceptor
        if uint32(donor) notin ssgraph.ss_to_index:
          ssgraph.ss_to_index[uint32(donor)] = uint32(ssgraph.adjacencies.len)
          ssgraph.adjacencies.add(@[isoform_adjacency_index])
          ssgraph.previously_visited.add(0)
        else:
          ssgraph.adjacencies[ssgraph.ss_to_index[uint32(donor)]].add(
            isoform_adjacency_index)
        if uint32(acceptor) notin ssgraph.ss_to_index:
          ssgraph.ss_to_index[uint32(acceptor)] = uint32(
            ssgraph.adjacencies.len)
          ssgraph.adjacencies.add(@[isoform_adjacency_index])
          ssgraph.previously_visited.add(0)
        else:
          ssgraph.adjacencies[ssgraph.ss_to_index[uint32(acceptor)]].add(
            isoform_adjacency_index)
        ssgraph.adjacencies[isoform_adjacency_index].add(
          ssgraph.ss_to_index[uint32(donor)])
        ssgraph.adjacencies[isoform_adjacency_index].add(
          ssgraph.ss_to_index[uint32(acceptor)])
    ### Fetch clusters
    let isoform_clusters = findIsoformClusters(addr ssgraph)
    # echo isoform_clusters
    var cluster_counter = 0
    var subcluster_counter = 0
    var splice_table_counter = 0
    var cluster_sizes : CountTable[int]
    var read_id_to_cluster : Table[string, int]
    for read_ids in spliced_table.values:
      for read_id in read_ids:
        read_id_to_cluster[read_id] = cluster_counter
      cluster_sizes.inc(cluster_counter,read_ids.len)
      # echo int(isoform_clusters[cluster_counter][subcluster_counter]),
      #   "\t", splice_table_counter
      assert int(isoform_clusters[cluster_counter][subcluster_counter]) ==
        splice_table_counter
      if subcluster_counter == isoform_clusters[cluster_counter].len - 1:
        subcluster_counter = 0
        cluster_counter += 1
      else:
        subcluster_counter += 1
      splice_table_counter += 1
    # echo cluster_to_read_id

    #Write cluster FASTA/Q files:
    echo cluster_sizes
    var bam : Bam
    discard open(bam,opt.file,index=true)
    if opt.reference != "":
      var fai : Fai
      discard open(fai,opt.reference)
      output_file_count += writeFASTAsFromBAM(bam,
                                              addr read_id_to_cluster,
                                              addr cluster_sizes,
                                              "chrI",
                                              "clusters_out/cluster_",
                                              output_file_count,
                                              fai)
    else:
      output_file_count += writeFASTAsFromBAM(bam,
                                              addr read_id_to_cluster,
                                              addr cluster_sizes,
                                              "chrI",
                                              "clusters_out/cluster_",
                                              output_file_count)
    #TODO - add logic to write FASTQs
    #TODO - add option to 'correct' reads based on reference genome
    bam.close

main()