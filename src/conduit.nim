import os
import times
import osproc
import parseopt
import strutils
import strformat
import poGraphUtils
import tables
import heapqueue
import threadpools/threadpool_simple as tps
import hts
import sets
import poaV2/header
import poaV2/poa
import fasta
import fastq

{.experimental.}

type
  ConduitOptions = object
    runFlag : bool
    finalPolish : bool
    intermediates : bool
    mode : string
    clustersDirectory : string
    nanoporeFormat : string
    illuminaFormat : string
    u2t : bool
    bowtieStrandConstraint : string
    bowtieAlignmentMode : string
    bowtieReadInputs : string
    scoreMatrixPath : string
    outputDir : string
    tmpDir : string
    files : seq[string]
    trims : seq[string]
    isoformDelta : uint64
    endsDelta : uint64
    illuminaWeight : uint64
    threadNum : uint64
    maxIterations : uint64
    stringent : bool
    stringentTolerance : int
    samtoolsMemory : string
    maxAlignments : uint64

#Minor TODOs:
#TODO - convert from passing tuple back to passing vars individually; relic of older threading approach

#Major TODOs (Future release versions?):
#TODO - Add clustering tool that runs with the benefit of a reference genome ( Cluster based on splice junctions / or overlap of alignments )
#TODO - Move filtering step to a graph based polishing step, the linear step has a problem of too low coverage or too many isoforms leading to isoforms not being completely covered by illumina reads when they should be.
#TODO - Add less stringent filtering step that extracts out the longest contiguous region covered by Illumina reads (with ends tolerance?) should be easy to do and reduce # of false negatives.
#TODO - Add support for duplicate read ID's that doesn't break everything.
#TODO - Add option to output both sensitive and stringent isoform sets in the same run
#TODO - .gz support for nanopore scaffolds. Can probably do what Trinity does and just add a decompress step for generating temp files.
#TODO - Add output indicating completion percentage for each iteration. Use https://github.com/euantorano/progress.nim ?
#TODO - Add mode that continues polishing where a previous run left off
#TODO - Rewrite poa in nim(?) - Probably faster as the C code, only advantage to the rewrite is it makes doing clever things with the poa easier down the line. (Unless we did SIMD poa, which would require learning nim SIMD, and probably step on Eyras Lab's toes)
#TODO - Clustering mode for when you DO have a reference genome? Or is that too similar to Stringtie2 to be worth doing? -- Probably too similar?
#TODO - Break poParser into smaller .nim files with more accurate and descriptive names
#TODO - Figure out if there's anything to be done about the left-aligned problem inherent to the partial order graph based correction? (thereby allowing us to get rid of linear polishing step)
#TODO - Add option to output .po files (or a new format, a multi-po file) to explain the relationships between isoforms.
#TODO - Add quantification output for each isoform / gene
#TODO - Add options that allow you to specify the path to the bowtie2 and samtools binaries


proc conduitVersion() : string =
  return "CONDUIT Version 0.1.0 by Nathan Roach ( nroach2@jhu.edu, https://github.com/NatPRoach/conduit/ )"

proc writeDefaultHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitVersion()
  echo "Usage:"
  echo "  ./conduit <nano | hybrid>"
  echo "NOTE: nano mode not yet implemented... coming soon"
  echo "      to run the equivalent of nano mode, run hybrid mode with -i:0 and --no-final-polish"
  echo "      this will still require 'illumina' files to be passed, but will not check that they actually exist"
  echo "      so including -U this_file_does_not_exist.fq should run."

proc writeNanoHelp() =
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitVersion()
  echo "Usage:"
  echo "  ./conduit nano [options] <clustersDirectory>"
  echo "  <clustersDirectory>   Directory containing the .fasta/.fa or .fastq/.fq files of reads separated by gene cluster"
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

proc writeHybridHelp() = 
  echo "CONDUIT - CONsensus Decomposition Utility In Transcriptome-assembly:"
  echo conduitVersion()
  echo "Usage:"
  echo "  ./conduit hybrid [options] <clustersDirectory> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>}"
  echo "  <clustersDirectory>   Directory containing the .fasta/.fa or .fastq/.fq files of reads separated by gene cluster"
  echo "                         NOTE: .gz support coming for nanopore scaffold data, but is not an option at this time"
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
  echo "        Scaffold reads are stranded forward relative to coding strand, enforces --UtoT"
  echo "    --cdna-rev-stranded"
  echo "        Scaffold reads are stranded reverse complemented relative to coding strand"
  echo "    --cdna"
  echo "        Scaffold reads are NOT stranded"
  echo "    --sfq (default)"
  echo "        Scaffold reads are in FASTQ format, enforces --UtoT"
  echo "    --sfa"
  echo "        Scaffold reads are in FASTA format"
  echo "    --UtoT (default)"
  echo "        Scaffold reads contain Us instead of Ts. Converts U nucleotides to Ts in the sequences"
  echo "        NOTE: This adds a bit of I/O overhead but doesn't affect things if your sequences are already U free"
  echo "    --noUtoT"
  echo "        Scaffold reads do not contain Us and do not need to be converted."
  # echo "    --duplicate-filter"
  # echo "        Scaffold reads have duplicate read IDs, filter out the duplicate reads"
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
  echo "        Example formatting for the score matrix can be found at poaV2/myNUC3.4.4.mat"
  echo "    -d, --isoform-delta (35)"
  echo "        Maximum indel size to be 'corrected', beyond this size a new isoform is declared. Must be between 2 and 255"
  echo "    -e, --ends-delta (35)"
  echo "        Maximum size at the ends of isoforms to 'correct' before splitting. Must be between 2 and 255"
  echo "    -i, --max-iterations (5)"
  echo "        Maximum number of iterations to align to and correct scaffolds. Does not include optional final polshing step"
  echo "        Note: Providing a value of 0 will not perform any graph based illumina correction"
  echo "    -w, --illumina-weight (10)"
  echo "        Weight of illumina reads relative to nanopore reads when generating consensus"
  echo "    --final-polish (default)"
  echo "        Include a final correction of individual isoforms, not in a splice graph"
  echo "    --no-final-polish"
  echo "        Do not do a final correction of individual isoforms, not in a splice graph"
  echo "    --stringent (default)"
  echo "        Enforce that every base / edge in each final reported isoform is supported by an Illumina read, excluding --stringent-tolerance bp on each end of each isoform"
  echo "    --no-stringent"
  echo "        Do not enforce that every base / edge in each final reported isoform is supported by an Illumina read"
  echo "    --stringent-tolerance (100)"
  echo "        Number of bases at the end of each isoform that do not have to have Illumina reads supporting them when run in --stringent mode; ignored when run with --no-stringents"
  # echo "    --scaffold-minimum (1)" # TODO
  # echo "        Minimum number of scaffolding reads supporting an isoform necessary to report the isoform in the final output"
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
  echo "    -k,--bowtie2-max-alignments (50)"
  echo "        Maximum number of alignments per Illumina read to be used in final polishing step"
  # echo "    --bowtie2-path (bowtie2)" # TODO
  echo "  SAMtools:"
  echo "    --samtools-thread-memory (768 MiB)"
  echo "        Memory amount to use per SAMtools thread"
  echo "        Specified either in bytes or with a K, M, or G suffix"
  # echo "    --samtools-path (samtools)" # TODO
  echo "  Miscellaneous:"
  echo "    -h, --help"
  echo "        Display this help message and exit"
  echo "    -v, --version"
  echo "        Display the installed version number of CONDUIT and exit"
  echo "    --tmp-dir <path> (conduit-tmp/)"
  echo "        <path> where temporary files will be created"
  echo "    -t, --threads (4)"
  echo "        Number of threads to run in parallel (used for both Bowtie2 and Partial Order Graph correction)"
  # echo "        NOTE: Providing a value of 0 will attempt to autodetect the number of CPUs availible and use that." # TODO
  # echo "              If CPU number cannot be detected, the default of 4 threads will be used. # TODO

proc removeFiles(files : openArray[string]) =
  for file in files:
    removeFile(file)

proc createDirs(dirs : openArray[string]) = 
  for dir in dirs:
    createDir(dir)

proc returnFalse() : bool {.thread.} = 
  return false

proc outputTiming(outfilepath : string,time_seq : seq[Time],opts : ConduitOptions) =
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  outfile.write("CONDUIT Timing Log:\n")
  for i in 0..<time_seq.len - 2:
    outfile.write(&"Iter {i}: {(time_seq[i + 1] - time_seq[i]).inSeconds} s\n")
  if opts.finalPolish:
    outfile.write(&"Final Polish: {(time_seq[^1] - time_seq[^2]).inSeconds} s\n")
  outfile.write(&"Total: {(time_seq[^1] - time_seq[0]).inSeconds} s\n")
  outfile.close()

proc outputSettings(outfilepath : string,opts : ConduitOptions) = 
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  outfile.write("CONDUIT Command Log:\n")
  outfile.write("Working directory:\n")
  outfile.write(&"    {getCurrentDir()}\n")
  outfile.write("Command:\n")
  let commandParams = commandLineParams().join(sep = " ")
  outfile.write(&"    {getAppFilename()} {commandParams}\n")
  outfile.write("CONDUIT Settings:\n")
  outfile.write(&"    mode : {opts.mode}\n")
  outfile.write(&"    rounds of graph based polishing : {opts.maxIterations}\n")
  outfile.write(&"    final round of linear polishing : {opts.finalPolish}\n")
  outfile.write(&"    stringent mode : {opts.stringent}\n")
  if opts.stringent:
    outfile.write(&"    stringent tolerance : {opts.stringentTolerance}\n")
  outfile.write(&"    save intermediates : {opts.intermediates}\n")
  outfile.write(&"    scaffold directory : {opts.clustersDirectory}\n")
  outfile.write(&"    scaffold file type : {opts.nanoporeFormat}\n")
  outfile.write(&"    Convert U to T : {opts.u2t}\n")
  outfile.write(&"    Illumina files : {opts.bowtieReadInputs}\n")
  outfile.write(&"    Illumina file type : {opts.illuminaFormat}\n")
  outfile.write(&"    isoform delta : {opts.isoformDelta}\n")
  outfile.write(&"    ends delta : {opts.endsDelta}\n")
  outfile.write(&"    Illumina weight : {opts.illuminaWeight}\n")
  outfile.write(&"    threads : {opts.threadNum}\n")
  if opts.scoreMatrixPath == "":
    outfile.write( "    score matrix : default\n")
  else:
    outfile.write(&"    score matrix : {opts.scoreMatrixPath}\n")
  outfile.write(&"    output directory : {opts.outputDir}\n")
  outfile.write(&"    temporary directory : {opts.tmpDir}\n")
  outfile.write("SAMtools Settings\n")
  outfile.write(&"    SAMtools thread memory : {opts.samtoolsMemory}\n")
  outfile.write("Bowtie2 Settings:\n")
  outfile.write(&"    bowtie2 strand constraint : {opts.bowtieStrandConstraint}\n")
  outfile.write(&"    bowtie2 alignment mode : {opts.bowtieAlignmentMode}\n")
  outfile.write(&"    bowtie2 max alignments : {opts.maxAlignments}\n")
  outfile.close()


proc mergeFiles(infilepaths : openArray[string], outfilepath : string, delete_old_files : bool = false) = 
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  for infilepath in infilepaths:
    var infile : File
    discard open(infile,infilepath,fmRead)
    outfile.write(infile.readAll)
    infile.close()
  if delete_old_files:
    removeFiles(infilepaths)
  outfile.close()


proc runPOAandCollapsePOGraph(intuple : (string,string,string,string,uint16,uint16,bool)) {.thread.} =
  let (infilepath,outdir,matrixFilepath,format,isoformDelta,endsDelta,u2t) = intuple
  let trim = infilepath.split(os.DirSep)[^1].split(".")[0]
  var fastaFile : string
  if format == "fasta":
    if u2t:
      fastaFile = &"{outdir}{trim}.tmp.fa"
      convertUtoTinFASTA(infilepath,fastaFile)
    else:
      fastaFile = infilepath
  elif format == "fastq":
    fastaFile = &"{outdir}{trim}.tmp.fa"
    convertFASTQfileToFASTAfile(infilepath,fastaFile)
  var splitNum = 200
  var (numFastas,_) = splitFASTA2(fastaFile,&"{outdir}{trim}.tmp",splitNum = splitNum)
  var totalFastas = numFastas
  var deleteFastaFlag = false
  if format == "fastq" or u2t or numFastas > 1:
    deleteFastaFlag = true
  if numFastas > 1:
    removeFile(fastaFile)
    # Cacluate representative reads for each subfasta, store each in separate consensus fasta file
    var mergeQueue : HeapQueue[(int, int)]
    for i in 0..<numFastas:
      let outFASTAfilepath = &"{outdir}{trim}.tmp_consensus{i}.fa"
      var outFASTAfile : File
      discard open(outFASTAfile,outFASTAfilepath,fmWrite)
      let tmpFasta = &"{outdir}{trim}.tmp_subfasta{i}.fa"
      var seqFile : PFile = fopen(cstring(tmpFasta))
      var po = getPOGraphFromFasta(seqFile,cstring(matrixFilepath),cint(1),matrix_scoring_function)
      removeFile(tmpFasta)
      var trimPo = poGraphUtils.convertPOGraphtoTrimmedPOGraph(po)
      var (representativePaths,readSupports) = poGraphUtils.getRepresentativePaths3(addr trimPo, psi = isoformDelta, endsDelta = endsDelta)
      let consensusPo = poGraphUtils.buildConsensusPO(addr po, representativePaths, readSupports, &"{trim}.tmp_subfasta{i}")
      poGraphUtils.writeCorrectedReads(consensusPo,outFASTAfile)
      outFASTAfile.close()
      mergeQueue.push((consensusPo.reads.len,i))
    var (_,fileIdx1) = mergeQueue.pop()
    var (_,fileIdx2) = mergeQueue.pop()
    
    var newFilepath = ""
    while true:
      # Figure out which 2 files we're merging
      let filepath1 = &"{outdir}{trim}.tmp_consensus{fileIdx1}.fa"
      let filepath2 = &"{outdir}{trim}.tmp_consensus{fileIdx2}.fa"
      
      # Merge the files
      let newTmpFilepath = &"{outdir}{trim}.tmp_subfasta{totalFastas}.fa"
      mergeFiles([filepath1,filepath2], newTmpFilepath, delete_old_files=true)

      # Decompose the new temp file
      var seqFile : PFile = fopen(cstring(newTmpFilepath))
      var po = getPOGraphFromFasta(seqFile,cstring(matrixFilepath),cint(1),matrix_scoring_function,weight_support = true)
      removeFile(newTmpFilepath)
      var trimPo = poGraphUtils.convertPOGraphtoTrimmedPOGraph(po)
      var (representativePaths,readSupports) = poGraphUtils.getRepresentativePaths3(addr trimPo, psi = isoformDelta, endsDelta = endsDelta)
      let consensusPo = poGraphUtils.buildConsensusPO(addr po, representativePaths, readSupports, &"{trim}.tmp_subfasta{totalFastas}")
      
      # Write the decomposition to a new consensus file
      newFilepath = &"{outdir}{trim}.tmp_consensus{totalFastas}.fa"
      var newFile : File
      discard open(newFile,newFilepath,fmWrite)
      poGraphUtils.writeCorrectedReads(consensusPo,newFile)
      newFile.close()
      
      if mergeQueue.len > 0:
        var tmp : int
        mergeQueue.push((consensusPo.reads.len,totalFastas))
        (tmp,fileIdx1) = mergeQueue.pop()
        (tmp,fileIdx2) = mergeQueue.pop()
        totalFastas += 1
      else:
        break
    moveFile(newFilepath, &"{outdir}fasta{os.DirSep}{trim}.consensus.fa")
  else:
    var seqFile : PFile = fopen(cstring(fastaFile), "r")
    var po2 = getPOGraphFromFasta(seqFile,cstring(matrixFilepath),cint(1),matrix_scoring_function)
    if deleteFastaFlag:
      removeFile(fastaFile)
    var trimPo2 = poGraphUtils.convertPOGraphtoTrimmedPOGraph(po2)
    var (representativePaths,readSupports) = poGraphUtils.getRepresentativePaths3(addr trimPo2, psi = isoformDelta, endsDelta = endsDelta)
    let consensusPo = poGraphUtils.buildConsensusPO(addr po2, representativePaths, readSupports, trim)
    let outFASTAfilepath = &"{outdir}fasta{os.DirSep}{trim}.consensus.fa"
    var outFASTAfile : File
    discard open(outFASTAfile,outFASTAfilepath,fmWrite)
    poGraphUtils.writeCorrectedReads(consensusPo,outFASTAfile)
    outFASTAfile.close()


proc runGraphBasedIlluminaCorrection(intuple : (string,string,string,uint64,uint16,uint16)) : bool {.thread.} =
  let (tmpDir, trim, matrixFilepath, iter,isoformDelta,endsDelta) = intuple
  let lastFastaDir = &"{tmpDir}{iter-1}{os.DirSep}fasta{os.DirSep}"
  let thisFastaDir = &"{tmpDir}{iter}{os.DirSep}fasta{os.DirSep}"

  let lastFastaFilepath = &"{lastFastaDir}{trim}.consensus.fa"
  var seqFile : PFile = fopen(cstring(lastFastaFilepath), "r")
  var po = getPOGraphFromFasta(seqFile,matrixFilepath,cint(1),matrix_scoring_function,weight_support = true)

  let thisFastaFilepath = &"{thisFastaDir}{trim}.consensus.fa"

  let bamfilepath = &"{tmpDir}{iter-1}{os.DirSep}alignments.bam"

  var bam : Bam
  var trimPo = convertPOGraphtoTrimmedPOGraph(po)
  discard open(bam,bamfilepath,index=true)
  illuminaPolishPOGraph(addr trimPo, bam)
  var (representativePaths,readSupports) = getRepresentativePaths3(addr trimPo, psi = isoformDelta,endsDelta = endsDelta)
  var records = getFastaRecordsFromTrimmedPOGraph(addr trimPo, representativePaths, readSupports, trim)
  var outfile : File
  discard open(outfile,thisFastaFilepath,fmWrite)
  writeFASTArecordsToFile(outfile,records)
  # writeCorrectedReads(records,outfile)
  outfile.close()
  result = sameFileContent(lastFastaFilepath,thisFastaFilepath)
  removeFile(lastFastaFilepath)

proc runLinearBasedIlluminaCorrection(intuple : (string,string,uint64,uint64,uint16,bool,int)) {.thread.} = 
  let (tmpDir, trim, convergedIter, iter,isoformDelta,stringent,stringentTolerance) = intuple
  let lastFastaDir = &"{tmpDir}{convergedIter}{os.DirSep}fasta{os.DirSep}"
  let thisFastaDir = &"{tmpDir}{iter}{os.DirSep}fasta{os.DirSep}"

  let lastFastaFilepath = &"{lastFastaDir}{trim}.consensus.fa"
  let thisFastaFilepath = &"{thisFastaDir}{trim}.consensus.fa"

  let bamfilepath = &"{tmpDir}{iter-1}{os.DirSep}alignments.bam"

  var infile : File
  var bam : Bam
  discard open(infile,lastFastaFilepath)
  var reads = parseFasta(infile)
  infile.close()
  var corrected : seq[FastaRecord]
  discard open(bam,bamfilepath,index=true)
  for read in reads:
    var trimPo = getTrimmedGraphFromFastaRecord(read)
    illuminaPolishPOGraph(addr trimPo, bam,debug=true)
    discard getRepresentativePaths3(addr trimPo, psi = isoformDelta)
    if (not stringent) or stringencyCheck(addr trimPo,trimPo.reads[0].correctedPath,stringentTolerance = stringentTolerance):
      corrected.add(FastaRecord(readId : read.readId, sequence : getSequenceFromPath(trimPo,trimPo.reads[0].correctedPath)))
  var outfile : File
  discard open(outfile,thisFastaFilepath,fmWrite)
  writeFASTArecordsToFile(outfile,corrected)
  # writeCorrectedReads(corrected,outfile)
  outfile.close()


proc combineFiles(indirectory : string, intrims : openArray[string], outfilepath : string) = 
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  for i,trim in intrims:
    let filepath = &"{indirectory}{trim}.consensus.fa"
    if fileExists(filepath):
      var file : File
      echo &"appending {filepath}"
      discard open(file,filepath,fmRead)
      outfile.write(file.readAll)
      file.close()
    else:
      echo &"{filepath} doesn't exist, poaV2 went wrong with that cluster"
  outfile.close()

proc combineFilesIntermediate(indirectory : string, intrims : openArray[string], outfilepath : string, lastCorrection : Table[int,int]) = 
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  for i,trim in intrims:
    if i in lastCorrection:
      continue
    let filepath = &"{indirectory}{trim}.consensus.fa"
    if fileExists(filepath):
      var file : File
      echo &"appending {filepath}"
      discard open(file,filepath,fmRead)
      outfile.write(file.readAll)
      file.close()
    else:
      echo &"{filepath} doesn't exist, poaV2 went wrong with that cluster..."
  outfile.close()

proc combineFilesFinal(tmp_directory : string,last_num : uint64, intrims : openArray[string], outfilepath : string, lastCorrection : Table[int,int]) = 
  var outfile : File
  discard open(outfile,outfilepath,fmWrite)
  for i,trim in intrims:
    var last = last_num
    if i in lastCorrection:
      last = uint64(lastCorrection[i]) #TODO convert lastCorrection to uint64 types.
    let indirectory = &"{tmp_directory}{last}{os.DirSep}fasta{os.DirSep}"
    let filepath = &"{indirectory}{trim}.consensus.fa"
    if fileExists(filepath):
      var file : File
      echo &"appending {filepath}"
      discard open(file,filepath,fmRead)
      outfile.write(file.readAll)
      file.close()
    else:
      echo &"{filepath} doesn't exist, poaV2 went wrong with that cluster..."
  outfile.close()

proc getBowtie2options(opt : ConduitOptions, indexPrefix, sam : string, finalPolish : bool = false) : seq[string] = 
  var arguments : seq[string]
  arguments.add("--xeq")
  arguments.add("--no-unal")
  arguments.add("-p")
  arguments.add(&"{opt.threadNum}")
  if finalPolish:
    arguments.add("-k")
    arguments.add(&"{opt.maxAlignments}")
  arguments.add(opt.bowtieStrandConstraint)
  arguments.add(opt.bowtieAlignmentMode)
  arguments.add("--n-ceil")
  arguments.add("L,0,0")
  arguments.add("-x")
  arguments.add(indexPrefix)
  for arg in opt.bowtieReadInputs.split(" "):
    arguments.add(arg)
  arguments.add("-S")
  arguments.add(sam)
  return arguments

proc parseOptions() : ConduitOptions = 
  var clustersDirectory = ""
  var clustersDirectoryFlag = false

  var mate1s : seq[string]
  var mate2s : seq[string]
  var unpaireds : seq[string]
  var interleaved : seq[string]
  var bams : seq[string]

  var illuminaStrand="reverse"
  var illuminaStrandFlag = false

  var illuminaFormat="fastq"
  var illuminaFormatFlag = false

  var nanoporeType="drna"
  var nanoporeTypeFlag = false

  var nanoporeStrand="forward"
  # var nanopore_strand_flag = false

  var nanoporeFormat="fastq"
  var nanoporeFormatFlag = false

  var u2t = true
  var u2tFlag = false

  var scoreMatrixPath = ""
  var scoreMatrixPathFlag = false

  var isoformDelta = 35'u64
  var isoformDeltaFlag = false

  var endsDelta = 35'u64
  var endsDeltaFlag = false

  var maxIterations = 5'u64
  var maxIterationsFlag = false

  var illuminaWeight = 10'u64
  var illuminaWeightFlag = false

  var finalPolish = true
  var finalPolishFlag = false

  var stringent = true
  var stringentFlag = false

  var stringentTolerance = 100
  var stringentToleranceFlag = false

  var outputDir = &"conduit-out{os.DirSep}"
  var outputDirFlag = false

  var intermediates = false
  var intermediatesFlag = false

  var local = false
  var localFlag = false

  var maxAlignments = 50'u64
  var maxAlignmentsFlag = false

  var samtoolsMemory = "768M"
  var samtoolsMemoryFlag = false

  var helpFlag = true
  var versionFlag = false

  var tmpDir = &"conduit-tmp{os.DirSep}"
  var tmpDirFlag = false

  var threadNum = 4'u64
  var threadNumFlag = false

  var runFlag = true


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
      helpFlag = false
    i += 1
    case mode:
      of "nano", "hybrid":
        case kind:
          of cmdEnd:
            break
          of cmdShortOption, cmdLongOption:
            if last != "":
              echo &"ERROR - Option {last} provided without an argument"
              runFlag = false
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
                if not nanoporeTypeFlag:
                  nanoporeTypeFlag = true
                else:
                  runFlag = false
                  echo "ERROR - Multiple scaffold types input, choose one of \"--drna\", \"--cdna-rev-stranded\", or \"--cdna\""
              of "cdna-rev-stranded":
                if not nanoporeTypeFlag:
                  nanoporeTypeFlag = true
                  nanoporeType = "cdna"
                  nanoporeStrand = "reverse"
                else:
                  runFlag = false
                  echo "ERROR - Multiple scaffold types input, choose one of \"--drna\", \"--cdna-rev-stranded\", or \"--cdna\""
                  helpFlag = true
                  break
              of "cdna":
                if not nanoporeTypeFlag:
                  nanoporeTypeFlag = true
                  nanoporeType = "cdna"
                  nanoporeStrand = "unstranded"
                else:
                  runFlag = false
                  echo "ERROR - Multiple scaffold types input, choose one of \"--drna\", \"--cdna-rev-stranded\", or \"--cdna\""
                  helpFlag = true
                  break
              of "sfq":
                if not nanoporeFormatFlag:
                  #already fastq
                  nanoporeFormatFlag = true
                else:
                  runFlag = false
                  echo "ERROR - Multiple scaffold format input, choose one of FASTA (--sfa) or FASTQ (--sfq)"
                  helpFlag = true
                  break
              of "sfa":
                if not nanoporeFormatFlag:
                  #already fastq
                  nanoporeFormatFlag = true
                  nanoporeFormat = "fasta"
                else:
                  runFlag = false
                  echo "ERROR - Multiple scaffold format input, choose one of FASTA (--sfa) or FASTQ (--sfq)"
                  helpFlag = true
                  break
              of "UtoT":
                if not u2tFlag:
                  u2tFlag = true
                elif not u2t:
                  runFlag = false
                  echo "ERROR - Conflicting flags: --UtoT and --noUtoT both specified"
                  break
              of "noUtoT":
                if not u2tFlag:
                  u2tFlag = true
                  u2t = false
                elif u2t:
                  runFlag = false
                  echo "ERROR - Conflicting flags: --UtoT and --noUtoT both specified"
                  break
              of "m", "score-matrix":
                if not scoreMatrixPathFlag:
                  scoreMatrixPathFlag = true
                  if val != "":
                    scoreMatrixPath = val
                  else:
                    last = "score-matrix"
                else:
                  runFlag = false
                  echo "ERROR - Multiple score matrices provided"
                  break
              of "u", "unstranded":
                if not illuminaStrandFlag:
                  illuminaStrandFlag = true
                  illuminaStrand = "unstranded"
                else:
                  runFlag = false
                  echo "ERROR - Multiple illumina strand input, choose one of (-u,--unstranded), (-f,--fwd-stranded), (-r,--rev-stranded)"
                  helpFlag = true
                  break
              of "f", "fwd-stranded":
                if not illuminaStrandFlag:
                  illuminaStrandFlag = true
                  illuminaStrand = "forward"
                else:
                  runFlag = false
                  echo "ERROR - Multiple illumina strand input, choose one of (-u,--unstranded), (-f,--fwd-stranded), (-r,--rev-stranded)"
                  helpFlag = true
                  break
                continue
              of "r", "rev-stranded":
                if not illuminaStrandFlag:
                  illuminaStrandFlag = true
                else:
                  runFlag = false
                  echo "ERROR - Multiple illumina strand input, choose one of (-u,--unstranded), (-f,--fwd-stranded), (-r,--rev-stranded)"
                  helpFlag = true
                  break
              of "ifq":
                if not illuminaFormatFlag:
                  illuminaFormatFlag = true
                else:
                  runFlag = false
                  echo "ERROR - Multiple illumina format input, choose one of FASTA (--ifa), or FASTQ (--ifq)"
                  helpFlag = true
                  break
              of "ifa":
                if not illuminaFormatFlag:
                  illuminaFormatFlag = true
                  illuminaFormat = "fasta"
                else:
                  runFlag = false
                  echo "ERROR - Multiple illumina format input, choose one of FASTA (--ifa), or FASTQ (--ifq)"
                  helpFlag = true
                  break
              of "d", "isoform-delta":
                if not isoformDeltaFlag:
                  isoformDeltaFlag = true
                  if val != "":
                    isoformDelta = parseUInt(val)
                  else:
                    last = "isoform-delta"
                else:
                  runFlag = false
                  echo "ERROR - Isoform delta specified multiple times"
                  break
              of "e", "ends-delta":
                if not endsDeltaFlag:
                  endsDeltaFlag = true
                  if val != "":
                    endsDelta = parseUInt(val)
                  else:
                    last = "ends-delta"
                else:
                  runFlag = false
                  echo "ERROR - Ends delta specified multiple times"
                  break
              of "i", "max-iterations":
                if not maxIterationsFlag:
                  maxIterationsFlag = true
                  if val != "":
                    maxIterations = parseUInt(val)
                  else:
                    last = "max-iterations"
                else:
                  runFlag = false
                  echo "ERROR - Max iterations specified multiple times"
                  break
              of "w", "illumina-weight":
                if not illuminaWeightFlag:
                  illuminaWeightFlag = true
                  if val != "":
                    illuminaWeight = parseUInt(val)
                  else:
                    last = "illumina-weight"
                else:
                  runFlag = false
                  echo "ERROR - Multiple Illumina weights specified"
                  break
              of "final-polish":
                if not finalPolishFlag:
                  finalPolishFlag = true
                elif not finalPolish:
                  runFlag = false
                  echo "ERROR - Conflicting flags: --final-polish and --no-final-polish both specified"
                  break
              of "no-final-polish":
                if not finalPolishFlag:
                  finalPolishFlag = true
                  finalPolish = false
                elif finalPolish:
                  runFlag = false
                  echo "ERROR - Conflicting flags: --final-polish and --no-final-polish both specified"
                  break
              of "stringent":
                if not stringentFlag:
                  stringentFlag = true
                elif not stringent:
                  runFlag = false
                  echo "ERROR - Conflicting flags: --stringent and --no-stringent both specified"
              of "no-stringent":
                if not stringentFlag:
                  stringentFlag = true
                  stringent = false
                elif stringent:
                  runFlag = false
                  echo "ERROR - Conflicting flags: --stringent and --no-stringent both specified"
              of "stringent-tolerance":
                if not stringentToleranceFlag:
                  stringentToleranceFlag = true
                  if val != "":
                    stringentTolerance = parseInt(val)
                  else:
                    last = "stringent-tolerance"
              of "n", "no-intermediates":
                if not intermediatesFlag:
                  intermediatesFlag = true
                elif intermediates:
                  runFlag = false
                  echo "ERROR - Conflicting flags: --no-intermediates and --save-intermediates both specified"
                  break
              of "s", "save-intermediates":
                if not intermediatesFlag:
                  intermediatesFlag = true
                  intermediates = true
                elif not intermediates:
                  runFlag = false
                  echo "ERROR - Conflicting flags: --no-intermediates and --save-intermediates both specified"
                  break
              of "o", "output-dir":
                if not outputDirFlag:
                  outputDirFlag = true
                  if val != "":
                    outputDir = val
                  else:
                    last = "output-dir"
                else:
                  runFlag = false
                  echo "ERROR - Output directory specified multiple times"
                  break
              of "end-to-end":
                if not localFlag:
                  localFlag = true
                  local = false
                else:
                  if local:
                    runFlag = false
                    echo "ERROR - Conflicting flags: --local and --end-to-end both specified"
                    break
              of "local":
                if not localFlag:
                  localFlag = true
                  local = true
                else:
                  if not local:
                    runFlag = false
                    echo "ERROR - Conflicting flags: --local and --end-to-end both specified"
                    break
              of "k","bowtie2-max-alignments":
                if not maxAlignmentsFlag:
                  maxAlignmentsFlag = true
                  if val == "":
                    last = "bowtie2-max-alignments"
                  else:
                    maxAlignments = parseUInt(val)
                else:
                  runFlag = false
                  echo "ERROR - Max bowtie2 alignments specified multiple times"
                  break
              of "samtools-thread-memory":
                if not samtoolsMemoryFlag:
                  samtoolsMemoryFlag = true
                  if val == "":
                    last = "samtools-thread-memory"
                  else:
                    samtoolsMemory = val
                else:
                  runFlag = false
                  echo "ERROR SAMtools thread memory specified multiple times"
                  break

              of "t", "threads":
                if not threadNumFlag:
                  if val == "":
                    last = "threads"
                  else:
                    threadNum = parseUInt(val)
                else:
                  runFlag = false
                  echo "ERROR - Threads specified multiple times"
                  break
              of "tmp-dir":
                if not tmpDirFlag:
                  tmpDirFlag = true
                  if val == "":
                    last = "tmp-dir"
                  else:
                    tmpDir = val
                else:
                  runFlag = false
                  echo "ERROR - Temporary directory specified twice"
                  break
              of "h", "help":
                helpFlag = true
                runFlag = false
                break
              of "v", "version":
                versionFlag = true
                runFlag = false
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
                scoreMatrixPath = key
              of "isoform-delta":
                isoformDelta = parseUInt(key)
              of "ends-delta":
                endsDelta = parseUInt(key)
              of "max-iterations":
                maxIterations = parseUInt(key)
              of "illumina-weight":
                illuminaWeight = parseUInt(key)
              of "stringent-tolerance":
                stringentTolerance = parseInt(key)
              of "output-dir":
                outputDir = key
              of "bowtie2-max-alignments":
                maxAlignments = parseUInt(key)
              of "samtools-thread-memory":
                samtoolsMemory = key
              of "threads":
                threadNum = parseUInt(key)
              of "tmp-dir":
                if key[^1] == os.DirSep:
                  tmpDir = key
                else:
                  tmpDir = key & os.DirSep
              of "":
                # echo "Argument: ", key
                if not clustersDirectoryFlag:
                  clustersDirectoryFlag = true
                  clustersDirectory = key
                else:
                  runFlag = false
                  echo "ERROR - Argument provided without associated option; please provide one <clustersDirectory> only"
                  break
              else:
                echo &"ERROR - unknown option {last} provided"
                runFlag = false
                break
            last = ""
      else:
        helpFlag = true
        break
  if clustersDirectory == "":
    echo "ERROR - <clusters directory> must be specified"
    helpFlag = true
    runFlag = false

  if helpFlag:
    case mode:
      of "nano":
        writeNanoHelp()
      of "hybrid":
        writeHybridHelp()
      else:
        echo "ERROR - first argument must specify correction mode, \"nano\" or \"hybrid\""
        writeDefaultHelp()
  if versionFlag:
    echo conduitVersion()
  var files : seq[string]
  var bowtieStrandConstraint : string
  var bowtieAlignmentMode : string
  var bowtieReadInputs : string
  if runFlag and mode == "hybrid":
    # Determine strand relationship between nanopore and illumina reads:
    if (nanoporeStrand == "forward" and illuminaStrand == "reverse") or (nanoporeStrand == "reverse" and illuminaStrand == "forward"):
      bowtieStrandConstraint = "--nofw"
    elif (nanoporeStrand == "forward" and illuminaStrand == "forward") or (nanoporeStrand == "reverse" and illuminaStrand == "reverse"):
      bowtieStrandConstraint = "--norc"
    
    # Format illumina inputs in a manner readable by Bowtie2
    var bowtieMate1Inputs = ""
    var bowtieMate2Inputs = ""
    var bowtieUnpairedInputs = ""
    var bowtieInterleavedInputs = ""
    var bowtieBamInputs = ""
    if mate1s.len != 0:
      bowtieMate1Inputs = "-1 " & mate1s.join(",") & " "
    if mate2s.len != 0:
      bowtieMate2Inputs = "-2 " & mate2s.join(",") & " "
    if unpaireds.len != 0:
      bowtieUnpairedInputs = "-U " & unpaireds.join(",") & " "
    if interleaved.len != 0:
      bowtieInterleavedInputs = "--interleaved " & interleaved.join(",") & " "
    if bams.len != 0:
      bowtieBamInputs = "-b " & bams.join(",") & " "
    bowtieReadInputs = &"{bowtieMate1Inputs}{bowtieMate2Inputs}{bowtieUnpairedInputs}{bowtieInterleavedInputs}{bowtieBamInputs}"
    if bowtieReadInputs == "":
      echo "ERROR - No Illumina data provided"
      runFlag = false
    
    if local:
      bowtieAlignmentMode = "--very-sensitive-local"
    else:
      bowtieAlignmentMode = "--very-sensitive"
    
    if nanoporeFormat == "fasta":
      for file in walkFiles(&"{clustersDirectory}*.fa"):
        files.add(file)
      for file in walkFiles(&"{clustersDirectory}*.fasta"):
        files.add(file)
      if files.len == 0:
        runFlag = false
        echo &"ERROR - No files of type .fa or .fasta found in <clusters directory> {clustersDirectory}"
        echo "NOTE: We don't currently support .gzip'd or bzip2'd scaffold files, though support for these formats is coming"
    elif nanoporeFormat == "fastq":
      for file in walkFiles(&"{clustersDirectory}*.fq"):
        files.add(file)
      for file in walkFiles(&"{clustersDirectory}*.fastq"):
        files.add(file)
      if files.len == 0:
        runFlag = false
        echo &"ERROR - No files of type .fq or .fastq found in <clusters directory> {clustersDirectory}"
        echo "NOTE: We don't currently support .gzip'd or bzip2'd scaffold files, though support for these formats is coming"
    if isoformDelta > 255'u64 or isoformDelta < 2'u64:
      echo "ERROR - Isoform delta must be between 2 and 255"
      runFlag = false
    elif isoformDelta < 15'u64:
      echo "WARNING - Low isoform delta values will increase the number of distinct isoforms and dramatically increase runtime"
      echo "          An isoform delta value of 15 or above is reccomended"
    if endsDelta > 255'u64 or endsDelta < 2'u64:
      echo "ERROR - Ends delta must be between 2 and 255"
      runFlag = false
    if threadNum < 1'u64:
      echo "ERROR - Must run at least 1 thread"
      runFlag = false
    if illuminaWeight == 0'u64:
      echo "ERROR - Illumina weight must be greater than 0"
      runFlag = false
    elif illuminaWeight < 5'u64:
      echo "WARNING - We reccomend weighing Illumina reads by at least 5x their long-read counterparts"
    if stringentTolerance < 3:
      echo "ERROR - Stringent tolerance cannot be less than 3"
      runFlag = false
    elif stringentTolerance < 25:
      echo "WARNING - Reads map poorly to the ends of long read scaffolds, we reccomend a stringency tolerance of at least 25"
    
  var trims = newSeq[string](files.len)
  for i,infilepath in files:
    trims[i] = infilepath.split(os.DirSep)[^1].split(".")[0]
  return ConduitOptions(runFlag : runFlag,
                        finalPolish : finalPolish,
                        intermediates : intermediates,
                        mode : mode,
                        clustersDirectory : clustersDirectory,
                        nanoporeFormat : nanoporeFormat,
                        illuminaFormat : illuminaFormat,
                        u2t : u2t,
                        bowtieStrandConstraint : bowtieStrandConstraint,
                        bowtieAlignmentMode : bowtieAlignmentMode,
                        bowtieReadInputs : bowtieReadInputs,
                        scoreMatrixPath : scoreMatrixPath,
                        outputDir : &"{outputDir}{os.DirSep}",
                        tmpDir : &"{tmpDir}{os.DirSep}",
                        files : files,
                        trims : trims,
                        isoformDelta : isoformDelta,
                        endsDelta : endsDelta,
                        maxIterations : maxIterations,
                        illuminaWeight : illuminaWeight,
                        threadNum : threadNum,
                        stringent : stringent,
                        stringentTolerance : stringentTolerance,
                        samtoolsMemory : samtoolsMemory,
                        maxAlignments : maxAlignments)

proc main() =
  let opt = parseOptions()
  if opt.mode == "hybrid" and opt.runFlag:
    
    var iterTimes : seq[Time]
    iterTimes.add(getTime())

    if not dirExists(opt.outputDir):
      createDir(opt.outputDir)
    outputSettings(&"{opt.outputDir}CONDUIT.settings", opt)

    var tmpAlreadyExisted = true
    if not dirExists(opt.tmpDir):
      tmpAlreadyExisted = false
      createDir(opt.tmpDir)

    let directoryNumber = opt.maxIterations + uint64(opt.finalPolish)
    for i in 0..directoryNumber:
      createDirs([&"{opt.tmpDir}{i}{os.DirSep}", &"{opt.tmpDir}{i}{os.DirSep}fasta{os.DirSep}"])

    let p = tps.newThreadPool(int(opt.threadNum))
    for file in opt.files:
      p.spawn runPOAandCollapsePOGraph((file, &"{opt.tmpDir}0/", opt.scoreMatrixPath, opt.nanoporeFormat, uint16(opt.isoformDelta), uint16(opt.endsDelta),opt.u2t))
    p.sync()
    iterTimes.add(getTime())
    var lastCorrection : Table[int,int]
    for iter in 1..opt.maxIterations:
      let lastDir = &"{opt.tmpDir}{iter-1}{os.DirSep}"
      let lastFastaDir = &"{lastDir}fasta{os.DirSep}"

      # let cur_dir = &"{opt.tmpDir}{iter}/"

      var lastConsensus : string
      if opt.intermediates:
        lastConsensus = &"{opt.outputDir}conduit_consensuses_iter{iter-1}.fa"
      else:
        lastConsensus = &"{opt.tmpDir}conduit_consensuses_iter{iter-1}.fa"

      removeFile(lastConsensus)
      combineFilesIntermediate(lastFastaDir,opt.trims,lastConsensus,lastCorrection)

      let indexPrefix = &"{lastDir}bowtie2_index"
      echo execProcess("bowtie2-build", args =["--threads",&"{opt.threadNum}", lastConsensus, indexPrefix],options={poUsePath})
      
      removeFile(lastConsensus)
      if opt.intermediates:
        combineFilesFinal(opt.tmpDir,iter-1,opt.trims,lastConsensus,lastCorrection)
      
      let sam = &"{lastDir}alignments.sam"
      let arguments = getBowtie2options(opt,indexPrefix,sam)
      echo execProcess("bowtie2", args = arguments, options={poUsePath,poStdErrToStdOut})
      
      let bam = &"{lastDir}alignments.bam"
      # echo execProcess(&"samtools sort -@ {opt.threadNum} {sam} > {bam}", options={poEvalCommand,poUsePath})
      echo execProcess("samtools", args=["sort", "-@", &"{opt.threadNum}", "-o", bam, "-m", opt.samtoolsMemory, sam], options={poUsePath,poStdErrToStdOut})
      echo execProcess("samtools", args=["index", bam],options={poUsePath,poStdErrToStdOut})
      
      removeFiles([sam, &"{indexPrefix}.1.bt2", &"{indexPrefix}.2.bt2", &"{indexPrefix}.3.bt2", &"{indexPrefix}.4.bt2", &"{indexPrefix}.rev.1.bt2", &"{indexPrefix}.rev.2.bt2"])

      let p = tps.newThreadPool(int(opt.threadNum))
      var converged = newSeq[tps.FlowVar[bool]](opt.trims.len)
      for i,trim in opt.trims:
        if i in lastCorrection:
          converged[i] = p.spawn returnFalse()
          continue
        converged[i] = p.spawn runGraphBasedIlluminaCorrection((opt.tmpDir,trim,opt.scoreMatrixPath,iter,uint16(opt.isoformDelta),uint16(opt.endsDelta)))
      p.sync()

      for i,converge in converged:
        if converge.read():
          lastCorrection[i] = int(iter)

      removeFiles([bam, &"{bam}.bai"])
      iterTimes.add(getTime())
    if opt.finalPolish:
      let iter = opt.maxIterations + 1
      let lastDir = &"{opt.tmpDir}{iter-1}{os.DirSep}"
      # let lastFastaDir = &"{lastDir}fasta{os.DirSep}"
      
      var lastConsensus = ""
      if opt.intermediates:
        lastConsensus = &"{opt.outputDir}conduit_consensuses_iter{iter-1}.fa"
      else:
        lastConsensus = &"{opt.tmpDir}conduit_consensuses_iter{iter-1}.fa"
      removeFile(lastConsensus)
      combineFilesFinal(opt.tmpDir, iter-1, opt.trims, lastConsensus, lastCorrection)

      let indexPrefix = &"{lastDir}bowtie2_index"
      echo execProcess("bowtie2-build", args =["--threads",&"{opt.threadNum}", lastConsensus, indexPrefix],options={poUsePath})
      
      if not opt.intermediates:
        removeFile(lastConsensus)

      let sam = &"{lastDir}alignments.sam"
      let arguments = getBowtie2options(opt,indexPrefix,sam,finalPolish = true)
      echo execProcess("bowtie2", args = arguments, options={poUsePath,poStdErrToStdOut})
      
      let bam = &"{lastDir}alignments.bam"
      # echo execProcess(&"samtools sort -@ {opt.threadNum} {sam} > {bam}", options={poEvalCommand,poUsePath})
      echo execProcess("samtools", args=["sort", "-@", &"{opt.threadNum}", "-o", bam, "-m", opt.samtoolsMemory, sam], options={poUsePath,poStdErrToStdOut})
      echo execProcess("samtools", args=["index", bam],options={poUsePath,poStdErrToStdOut})
      
      removeFiles([sam, &"{indexPrefix}.1.bt2", &"{indexPrefix}.2.bt2", &"{indexPrefix}.3.bt2", &"{indexPrefix}.4.bt2", &"{indexPrefix}.rev.1.bt2", &"{indexPrefix}.rev.2.bt2"])

      let p = tps.newThreadPool(int(opt.threadNum))
      for i,trim in opt.trims:
        var convergedIter = iter - 1
        if i in lastCorrection:
          convergedIter = uint64(lastCorrection[i])
        p.spawn runLinearBasedIlluminaCorrection((opt.tmpDir,trim,convergedIter,iter,uint16(opt.isoformDelta),opt.stringent,opt.stringentTolerance))
      p.sync()

      removeFiles([bam, &"{bam}.bai"])
      # removeDir(lastFastaDir)
    
      let finalConsensusPath = &"{opt.outputDir}conduit_final_consensuses.fa"
      if fileExists(finalConsensusPath):
        removeFile(finalConsensusPath)
      var finalFastaDir = &"{opt.tmpDir}{directoryNumber}{os.DirSep}fasta{os.DirSep}"

      combineFiles(finalFastaDir, opt.trims, finalConsensusPath)
    else:
      let finalConsensusPath = &"{opt.outputDir}conduit_final_consensuses.fa"
      if fileExists(finalConsensusPath):
        removeFile(finalConsensusPath)

      combineFilesFinal(opt.tmpDir, directoryNumber, opt.trims, finalConsensusPath, lastCorrection)

    
    if not tmpAlreadyExisted:
      removeDir(opt.tmpDir)
    iterTimes.add(getTime())
    outputTiming(&"{opt.outputDir}CONDUIT.timing", iterTimes, opt)
  
main()