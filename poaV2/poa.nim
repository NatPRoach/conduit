###########################################################
###--------  Define Relevant Types from poaV2   --------###
###########################################################
import header
import ../poGraphUtils
import tables
import sequtils
import strutils
import strformat

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
                              seq_x: ptr LPOLetter_T;
                              seq_y: ptr LPOLetter_T;
                              m: ptr ResidueScoreMatrix_T): LPOScore_T {.importc: "matrix_scoring_function".}

proc free_return_result(lpo_result : LPOReturnResult_T ) {.importc: "free_return_result".}

proc build_lpo_from_fasta(seq_file : Pfile;
               score_matrix_filepath : cstring;
               use_global_alignment : cint;
               scoring_function : proc (a1 : cint; a2 : cint; a3 : ptr LPOLetter_T; a4 : ptr LPOLetter_T; a5 : ptr ResidueScoreMatrix_T) : LPOScore_T) : LPOReturnResult_T {.importc : "build_lpo_from_fasta".}

###########################################################
###---------  Utility functions to interface   ---------###
###########################################################

proc convertPOFormats(lpo_return_result : LPOReturnResult_T) : POGraph =
  let lpo = lpo_return_result.lpo_seqs
  let matrix = lpo_return_result.matrix
  let num_reads = lpo.nsource_seq
  var seqs = cast[ptr UncheckedArray[LPOSequence_T]](lpo)
  var source_seqs = cast[ptr UncheckedArray[LPOSourceInfo_T]](lpo.source_seq)
  var letters = cast[ptr UncheckedArray[LPOLetter_T]](lpo.letter)
  # echo seqs[0].length
  let num_nodes = seqs[0].length
  var reads : seq[Read]
  var nodes : seq[Node]
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  var paths : seq[seq[uint32]]
  var read_id_duplicate_check : Table[string,uint32]
  for i in 0..<num_reads:
    let name = source_seqs[i].name
    let length = source_seqs[i].length
    let string_name = toString(toSeq(name)).strip(chars = {'\0'})
    var name_pad = 0'u32
    if string_name in read_id_duplicate_check:
      name_pad = read_id_duplicate_check[string_name]
      read_id_duplicate_check[string_name] += 1'u32
    else:
      read_id_duplicate_check[string_name] = 1'u32
    reads.add(Read(name : &"{string_name}_{name_pad}",length : uint32(length)))
    paths.add(@[])
  for i in 0..<num_nodes:
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
    var align_ring = letters[i].align_ring
    var a_pair = -1'i32
    if align_ring != cint(i):
      a_pair = int32(align_ring)
    nodes.add(Node(nts : $nt,
                   align_ring_partner : a_pair,
                   visited : false,
                   end_node_flag : false,
                   start_node_flag : false,
                   nanopore_support : uint32(supports.len)))
  for path in paths:
    for i in 1..<path.len:
      let u = path[i-1]
      let v = path[i]
      if (u,v) in weights:
        weights[(u,v)] += 1'u32
      else:
        weights[(u,v)] = 1'u32
  
  free_return_result(lpo_return_result)

  for i in 0..<num_reads:
    reads[i].path = paths[i]
    reads[i].corrected_path = paths[i]
  return POGraph(nodes : nodes,
                 reads : reads,
                 edges : edges,
                 weights : weights,
                 og_nodes : uint32(num_nodes))
proc getPOGraphFromFasta*(seq_file : Pfile;
                         score_matrix_filepath : cstring;
                         use_global_alignment : cint;
                         scoring_function : proc (a1 : cint; a2 : cint; a3 : ptr LPOLetter_T; a4 : ptr LPOLetter_T; a5 : ptr ResidueScoreMatrix_T) : LPOScore_T) : POGraph =
  let lpo_return_result = build_lpo_from_fasta(seq_file,
                                               score_matrix_filepath,
                                               use_global_alignment,
                                               scoring_function)
  return convertPOFormats(lpo_return_result)