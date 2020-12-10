import os
import parseopt
import strformat
import tables
import hts
import algorithm
# import times

import na
import fasta
import fastq
import version
# import genomeKDE
# TODO - Multithreading
# import threadpools/threadpool_simple as tps


type
  ClusteringOptions = object
    file : string
    mode : char
    clusterZeroIntrons : bool
    outputDir : string
    prefix : string
    outType : string
    reference : string
    run : bool
    exitCode : int
  
  SpliceSiteGraph = object
    adjacencies : seq[seq[uint32]]
    ssToIndex : Table[uint32, uint32]
    indexToIndex : OrderedTable[uint32,uint32]
    previouslyVisited : seq[uint32]


# TODO -   Possibly change uint32 to uint64 where necessary.
# TODO - Reevaluate before publication.
proc conduitClusterVersion() : string =
  return &"CONDUIT Clustering Version {version.ConduitVersion} by Nathan Roac" &
      "h\n( nroach2@jhu.edu, https://github.com/NatPRoach/conduit/ )"


proc writeClusterHelp() = 
  echo "WARNING - CONDUIT clustering is actively being developed and may not"
  echo "          be accurately reflected by this help text."
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
  # echo "    -z (default)"
  # echo "        Cluster reads with zero introns by overlap"
  # echo "    -d"
  # echo "        Do not cluster reads with zero introns by overlap"
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
                           stranded : bool = false) :
                            ( char,
                              (uint64, uint64),
                              seq[(uint64, uint64)] ) = 
  if record.flag.unmapped:
      return
  let alignmentStart = uint64(record.start)
  var refPosition = uint64(alignmentStart)
  let cigar = record.cigar
  let flag = record.flag
  var strand : char
  var spliceSites : seq[(uint64, uint64)]
  # If we're not stranded treat everything as if it's (-) strand
  # else check for the actual strand
  if  ( not stranded ) or flag.has_flag(16'u16):
    strand = '-'
  else:
    strand = '+'
  for op in cigar:
    case int(op.op):
      of 0, 2, 7, 8:
        refPosition += uint64(op.len)
      of 3:
        let start = refPosition
        refPosition += uint64(op.len)
        spliceSites.add((start,refPosition))
      else:
        discard
  let alignmentEnd = refPosition
  result = (strand,(alignmentStart,alignmentEnd),spliceSites)


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
      # echo &"WARNING - At the moment secondary / supplementary alignments " &
      #   "are not supported"
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


proc parseClusteringOptions() : ClusteringOptions = 
  var file : string
  var fileFlag = false

  var outputFormat = "fastq"
  var outputFormatFlag = false

  var clusterMode = 's'
  var clusterModeFlag = false

  var singleExon = true
  var singleExonFlag = false

  var outputDir = &"clusters{os.DirSep}"
  var outputDirFlag = false

  var prefix = &"cluster_"
  var prefixFlag = false

  var referenceFilepath = ""
  var referenceFilepathFlag = false

  var runFlag = true
  var helpFlag = false
  var versionFlag = false
  var exitCode = QuitFailure

  var last = ""

  for kind, key, val in getopt():
    case kind:
      of cmdEnd:
        break
      of cmdShortOption, cmdLongOption:
        if last != "":
          echo &"ERROR - Option {last} provided without an argument"
          runFlag = false
          break
        
        case key:
          of "ss":
            if clusterModeFlag and ( clusterMode != 's' ):
              echo "ERROR - Conflicting flags --ss and --overlap"
              helpFlag = true
              break
            else:
              clusterModeFlag = true
              clusterMode = 's'
          of "overlap":
            if clusterModeFlag and ( clusterMode != 'o' ):
              echo "ERROR - Conflicting flags --ss and --overlap"
              helpFlag = true
              break
            else:
              clusterModeFlag = true
              clusterMode = 'o'
          of "z":
            if singleExonFlag and not singleExon:
              echo "ERROR - Conflicting flags -z and -d"
              helpFlag = true
              break
            else:
              singleExonFlag = true
              singleExon = true
          of "d":
            if singleExonFlag and singleExon:
              echo "ERROR - Conflicting flags -z and -d"
              helpFlag = true
              break
            else:
              singleExonFlag = true
              singleExon = false
          of "fq":
            if outputFormatFlag and outputFormat != "fastq":
              echo "ERROR - Conflicting flags --fq and --fa"
              helpFlag = true
              break
            else:
              outputFormatFlag = true
              outputFormat = "fastq"
          of "fa":
            if outputFormatFlag and outputFormat != "fasta":
              echo "ERROR - Conflicting flags --fq and --fa"
              helpFlag = true
              break
            else:
              outputFormatFlag = true
              outputFormat = "fasta"
          of "r", "reference":
            if not referenceFilepathFlag:
              referenceFilepathFlag = true
              if val != "":
                referenceFilepath = val
              else:
                last = "reference"
            else:
              echo "ERROR - two references provided"
              helpFlag = true
              break
          of "o", "output-dir":
            if not outputDirFlag:
              outputDirFlag = true
              if val != "":
                outputDir = val
              else:
                last = "output"
            else:
              echo "ERROR - two output dirs provided"
              helpFlag = true
              break
          of "p", "prefix":
            if not prefixFlag:
              prefixFlag = true
              if val != "":
                outputDir = val
              else:
                last = "output"
            else:
              echo "ERROR - two output dirs provided"
              helpFlag = true
              break
          of "h", "help":
            helpFlag = true
            runFlag = false
            exitCode = QuitSuccess
            break
          of "v", "version":
            helpFlag = false
            versionFlag = true
            runFlag = false
            exitCode = QuitSuccess
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
            outputDir = key
          of "prefix":
            prefix = key
          of "reference":
            referenceFilepath = key
          of "":
            if not fileFlag:
              file = key
            else:
              fileFlag = true
              echo "ERROR - multiple files provided"
              break
          else:
            echo &"ERROR - unknown option {last} provided"
            runFlag = false
            break
        last = ""

  if file == "" and not helpFlag:
    echo "ERROR - clustering requires an input BAM file"
    helpFlag = true
    runFlag = false
    exitCode = QuitFailure

  if helpFlag:
    writeClusterHelp()
  if versionFlag:
    echo conduitClusterVersion()

  return ClusteringOptions(file : file,
                           mode : clusterMode,
                           clusterZeroIntrons : singleExon,
                           outputDir : outputDir,
                           prefix : prefix,
                           outType : outputFormat,
                           reference : referenceFilepath,
                           run : runFlag,
                           exitCode : exitCode)


proc bfs(ssgraph : ptr SpliceSiteGraph,
         cluster_number : uint32,
         starting_node : uint32) : seq[uint32] = 
  var toVisit : seq[uint32]
  ssgraph[].previouslyVisited[starting_node] = cluster_number
  toVisit.add(starting_node)
  while toVisit.len != 0:
    let u = toVisit.pop()
    result.add(u)
    for adj in ssgraph[].adjacencies[u]:
      if ssgraph[].previouslyVisited[adj] == 0:
        ssgraph[].previouslyVisited[adj] = cluster_number
        toVisit.add(adj)


proc findIsoformClusters(ssgraph : ptr SpliceSiteGraph) : seq[seq[uint32]] =
  var clusterCounter = 1'u32
  for isoformAdjacencyIndex in ssgraph[].indexToIndex.keys:
    if ssgraph[].previouslyVisited[isoformAdjacencyIndex] == 0:
      var clusteredIsoforms : seq[uint32]
      let clusteredNodes = bfs(ssgraph,
                                clusterCounter,
                                isoformAdjacencyIndex)
      for node_id in clusteredNodes:
        if node_id in ssgraph[].indexToIndex:
          clusteredIsoforms.add(ssgraph[].indexToIndex[node_id])
      clusterCounter += 1
      result.add(clusteredIsoforms.sorted)


proc getWeightedSpliceJunctionLocations(
        sstable : ptr OrderedTable[seq[(uint64, uint64)],seq[string]]) :
        (seq[(uint32,uint32)],seq[(uint32,uint32)]) = 
  var donorTable : OrderedTable[uint32,uint32]
  var acceptorTable : OrderedTable[uint32,uint32]
  for sss, read_ids in sstable[].pairs:
    let weight = uint32(read_ids.len)
    for (donor,acceptor) in sss:
      if uint32(donor) in donorTable:
        donorTable[uint32(donor)] += weight
      else:
        donorTable[uint32(donor)] = weight
      if uint32(acceptor) in acceptorTable:
        acceptorTable[uint32(acceptor)] += weight
      else:
        acceptorTable[uint32(acceptor)] = weight
  for ss, weight in donorTable.pairs:
    result[0].add((ss,weight))
  for ss, weight in acceptorTable.pairs:
    result[1].add((ss,weight))


proc writeFASTXsFromBAM(bam : Bam,
                        readIdToCluster : ptr Table[string,int],
                        clusterSizes : ptr CountTable[int],
                        query : string,
                        cluster_prefix : string,
                        starting_count : int,
                        output_type : string) : int =
  var openFiles : Table[int, File]
  var writtenReads : CountTable[int]
  for record in bam.query(query):
    if record.flag.supplementary or record.flag.secondary:
      # echo &"WARNING - At the moment secondary / supplementary alignments " &
      #   "are not supported"
      continue
    if record.qname notin readIdToCluster[]:
      continue
    let clusterId = readIdToCluster[][record.qname]
    if clusterId notin openFiles:
      var file : File
      if output_type == "fasta":
        discard open(file,
                     &"{cluster_prefix}{clusterId + starting_count}.fa",
                     fmWrite)
      elif output_type == "fastq":
        discard open(file,
                     &"{cluster_prefix}{clusterId + starting_count}.fq",
                     fmWrite)
      result += 1
      openFiles[clusterId] = file
      # file.writeFASTArecordToFile(record)
      if output_type == "fasta":
        file.writeBamRecordToFASTAfile(record)
      elif output_type == "fastq":
        file.writeBamRecordToFASTQfile(record)
    else:
      # echo &"Appending to cluster file {clusterId + starting_count}"
      if output_type == "fasta":
        openFiles[clusterId].writeBamRecordToFASTAfile(record)
      elif output_type == "fastq":
        openFiles[clusterId].writeBamRecordToFASTQfile(record)
    writtenReads.inc(clusterId)
    # echo writtenReads[clusterId], "\t", clusterSizes[][clusterId]
    if writtenReads[clusterId] == clusterSizes[][clusterId]:
      # echo &"Closing cluster file {clusterId + starting_count}"
      openFiles[clusterId].close
      openFiles.del(clusterId)
  # Just in case of secondary / supplementary alignments
  for clusterId, open_file in openFiles.pairs:
    open_file.close


proc correctBamRecordWithGenome(record : Record, fai : Fai) : FastaRecord = 
  let summary = summaryFromBamRecord(record, true)
  var ntSequence : string
  var startBase = int(summary[1][0])
  var endBase = -1
  for (donor,acceptor) in summary[2]:
    endBase = int(donor) - 1
    let nts = fai.get(record.chrom,
                      startBase,
                      endBase)
    # if record.qname == "m54284U_191110_105540/34605360/ccs":
    #   echo record.qname, "\t", record.chrom, "\t", startBase, "\t", endBase
    #   echo nts
    ntSequence.add(nts)
    startBase = int(acceptor)
  endBase = int(summary[1][1]) - 1
  let nts = fai.get(record.chrom,
                    startBase,
                    endBase)
  # if record.qname == "m54284U_191110_105540/34605360/ccs":
  #   echo record.qname, "\t", record.chrom, "\t", startBase, "\t", endBase
  #   echo nts
  ntSequence.add(nts)
  if summary[0] == '-':
    ntSequence = ntSequence.revComp
  result = FastaRecord( readId : record.qname, sequence : ntSequence)


proc writeFASTAsFromBAM(bam : Bam,
                        readIdToCluster : ptr Table[string,int],
                        clusterSizes : ptr CountTable[int],
                        query : string,
                        cluster_prefix : string,
                        starting_count : int,
                        fai : Fai) : int =
  var openFiles : Table[int, File]
  var writtenReads : CountTable[int]
  for record in bam.query(query):
    if record.flag.supplementary or record.flag.secondary:
      # echo "WARNING - At the moment secondary / supplementary alignments " &
      #   "are not supported"
      continue
    if record.flag.unmapped:
      echo "WARNING - Unmapped read detected, skipping"
      continue
    let fastaRecord = correctBamRecordWithGenome(record, fai)
    if fastaRecord.readId notin readIdToCluster[]:
      continue
    let clusterId = readIdToCluster[][fastaRecord.readId]
    if clusterId notin openFiles:
      var file : File
      discard open(file,&"{cluster_prefix}{clusterId + starting_count}.fa",
        fmWrite)
      result += 1
      # echo &"Opening cluster file {clusterId + starting_count}"
      openFiles[clusterId] = file
      file.writeFASTArecordToFile(fastaRecord)
    else:
      # echo &"Appending to cluster file {clusterId + starting_count}"
      openFiles[clusterId].writeFASTArecordToFile(fastaRecord)
    writtenReads.inc(clusterId)
    # echo writtenReads[clusterId], "\t", clusterSizes[][clusterId]
    if writtenReads[clusterId] == clusterSizes[][clusterId]:
      # echo &"Closing cluster file {clusterId + starting_count}"
      openFiles[clusterId].close
      openFiles.del(clusterId)
  # Just in case of secondary / supplementary alignments
  for clusterId, open_file in openFiles.pairs:
    open_file.close


proc main() = 
  let opt = parseClusteringOptions()
  if not opt.run:
    quit(opt.exitCode)
  if not dirExists(opt.outputDir):
    createDir(opt.outputDir)
  var bam : Bam
  discard open(bam,opt.file,index=true)
  var fai : Fai
  discard open(fai,opt.reference)
  var outputFileCount = 0
  for i,target in bam.hdr.targets:
    let chr = target.name
    # echo "Starting ", chr
    #TODO - iterate over chromosomes, chunks of chromosomes?
    #TODO - multithread? - see how the memory and time looks like then decide?
    var (strandSplicedTable, strandNonsplicedTable) =
      collapsedSummariesFromBam(bam, chr, true)
    for strand, spliced_table in strandSplicedTable.mpairs:
      spliced_table.sort(cmpSpliced)
      ### Populate SpliceSiteGraph
      var ssgraph : SpliceSiteGraph
      var isoformNodeIndex = 0'u32
      # isoformNodeIndex.setBit(31'u8)
      for sss in spliced_table.keys:
        var isoformAdjacencyIndex = uint32(ssgraph.adjacencies.len)
        ssgraph.adjacencies.add(@[])
        ssgraph.previouslyVisited.add(0)
        ssgraph.indexToIndex[isoformAdjacencyIndex] = isoformNodeIndex
        isoformNodeIndex += 1
        for (donor,acceptor) in sss:
          # echo donor, "\t", acceptor
          if uint32(donor) notin ssgraph.ssToIndex:
            ssgraph.ssToIndex[uint32(donor)] = uint32(ssgraph.adjacencies.len)
            ssgraph.adjacencies.add(@[isoformAdjacencyIndex])
            ssgraph.previouslyVisited.add(0)
          else:
            ssgraph.adjacencies[ssgraph.ssToIndex[uint32(donor)]].add(
              isoformAdjacencyIndex)
          if uint32(acceptor) notin ssgraph.ssToIndex:
            ssgraph.ssToIndex[uint32(acceptor)] = uint32(
              ssgraph.adjacencies.len)
            ssgraph.adjacencies.add(@[isoformAdjacencyIndex])
            ssgraph.previouslyVisited.add(0)
          else:
            ssgraph.adjacencies[ssgraph.ssToIndex[uint32(acceptor)]].add(
              isoformAdjacencyIndex)
          ssgraph.adjacencies[isoformAdjacencyIndex].add(
            ssgraph.ssToIndex[uint32(donor)])
          ssgraph.adjacencies[isoformAdjacencyIndex].add(
            ssgraph.ssToIndex[uint32(acceptor)])
      ### Fetch clusters
      let isoformClusters = findIsoformClusters(addr ssgraph)
      var clusterCounter = 0
      var subclusterCounter = 0
      var spliceTableCounter = 0
      var clusterSizes : CountTable[int]
      var readIdToCluster : Table[string, int]
      for read_ids in spliced_table.values:
        for readId in read_ids:
          readIdToCluster[readId] = clusterCounter
        clusterSizes.inc(clusterCounter,read_ids.len)
        # echo int(isoformClusters[clusterCounter][subclusterCounter]),
        #   "\t", spliceTableCounter
        if subclusterCounter == isoformClusters[clusterCounter].len - 1:
          subclusterCounter = 0
          clusterCounter += 1
        else:
          subclusterCounter += 1
        spliceTableCounter += 1
      #Write cluster FASTA/Q files:
      if opt.reference != "":
        outputFileCount += writeFASTAsFromBAM(bam,
                                                addr readIdToCluster,
                                                addr clusterSizes,
                                                chr,
                                                &"{opt.outputDir}" &
                                                  &"{os.DirSep}" &
                                                  &"{opt.prefix}",
                                                outputFileCount,
                                                fai)
      else:
        outputFileCount += writeFASTXsFromBAM(bam,
                                                addr readIdToCluster,
                                                addr clusterSizes,
                                                chr,
                                                &"{opt.outputDir}" &
                                                  &"{os.DirSep}" &
                                                  &"{opt.prefix}",
                                                outputFileCount,
                                                opt.outType)
  bam.close
  # fai.close
main()