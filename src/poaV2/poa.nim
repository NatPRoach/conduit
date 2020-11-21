###########################################################
###--------  Define Relevant Types from poaV2   --------###
###########################################################
import header
import ../poGraphUtils
import tables
import sequtils
import strutils
# import strformat

###########################################################
###------  Import Relevant Functions from poaV2   ------###
###########################################################
proc toString(str: seq[char]): string =
  result = newStringOfCap(len(str))
  for ch in str:
    add(result, ch)

proc fopen*(filename : cstring, mode : cstring = "r") : PFile {.importc: "fopen", header: "<stdio.h>".}

proc matrix_scoring_function*(i: cint;
                              j: cint;
                              seq_x: ptr LPOLetterT;
                              seq_y: ptr LPOLetterT;
                              m: ptr ResidueScoreMatrixT): LPOScoreT {.importc: "matrix_scoring_function".}

proc free_return_result(lpo_result : LPOReturnResultT ) {.importc: "free_return_result".}

proc build_lpo_from_fasta(seq_file : PFile;
               score_matrix_filepath : cstring;
               use_global_alignment : cint;
               scoring_function : proc (a1 : cint; a2 : cint; a3 : ptr LPOLetterT; a4 : ptr LPOLetterT; a5 : ptr ResidueScoreMatrixT) : LPOScoreT) : LPOReturnResultT {.importc : "build_lpo_from_fasta".}

###########################################################
###---------  Utility functions to interface   ---------###
###########################################################

proc convertPOFormats(lpoReturnResult : LPOReturnResultT, weight_support : bool = false) : POGraph =
  let lpo = lpoReturnResult.lpoSeqs
  let matrix = lpoReturnResult.matrix
  let numReads = lpo.nsourceSeq
  var seqs = cast[ptr UncheckedArray[LPOSequenceT]](lpo)
  var sourceSeqs = cast[ptr UncheckedArray[LPOSourceInfoT]](lpo.sourceSeq)
  var letters = cast[ptr UncheckedArray[LPOLetterT]](lpo.letter)
  # echo seqs[0].length
  let numNodes = seqs[0].length
  var reads : seq[Read]
  var nodes : seq[Node]
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  var paths : seq[seq[uint32]]
  for i in 0..<numReads:
    let name = toString(toSeq(sourceSeqs[i].name)).strip(chars = {'\0'})
    let length = sourceSeqs[i].length
    var support = 0'u32
    if weight_support:
      support = uint32(parseUInt(name.split('_')[^1]))
    else:
      support = 1'u32
    reads.add(Read(name : name,length : uint32(length),support : support))
    paths.add(@[])
  for i in 0..<numNodes:
    # echo matrix.symbol
    edges[uint32(i)] = @[]
    let nt = matrix.symbol[ord(letters[i].letter)]
    var link = addr (letters[i].left)
    while not isNil(link) and link.ipos >= 0:
      edges[uint32(link.ipos)].add(uint32(i))
      link = link.more
    var supports : seq[uint32]
    var source = addr(letters[i].source)
    while not isNil(source):
      paths[source.iseq].add(uint32(i))
      supports.add(uint32(source.iseq))
      source = source.more
    var alignRing = letters[i].alignRing
    var aPair = -1'i32
    if alignRing != cint(i):
      aPair = int32(alignRing)
    nodes.add(Node(nts : $nt,
                   align_ring_partner : aPair,
                   visited : false,
                   end_node_flag : false,
                   start_node_flag : false,
                   nanopore_support : uint32(supports.len)))
  for j,path in paths:
    for i in 1..<path.len:
      let u = path[i-1]
      let v = path[i]
      if (u,v) in weights:
        weights[(u,v)] += reads[j].support
      else:
        weights[(u,v)] = reads[j].support
  
  free_return_result(lpoReturnResult)

  for i in 0..<numReads:
    reads[i].path = paths[i]
    reads[i].correctedPath = paths[i]
  return POGraph(nodes : nodes,
                 reads : reads,
                 edges : edges,
                 weights : weights,
                 og_nodes : uint32(numNodes))
proc getPOGraphFromFasta*(seq_file : PFile;
                         score_matrix_filepath : cstring;
                         use_global_alignment : cint;
                         scoring_function : proc (a1 : cint; a2 : cint; a3 : ptr LPOLetterT; a4 : ptr LPOLetterT; a5 : ptr ResidueScoreMatrixT) : LPOScoreT,
                         weight_support : bool = false) : POGraph =
  let lpoReturnResult = build_lpo_from_fasta(seq_file,
                                               score_matrix_filepath,
                                               use_global_alignment,
                                               scoring_function)
  return convertPOFormats(lpoReturnResult,weight_support = weight_support)