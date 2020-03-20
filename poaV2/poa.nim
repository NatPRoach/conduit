###########################################################
###--------  Define Relevant Types from poaV2   --------###
###########################################################
import header
import ../poParser
import tables
import sequtils
# import strformat

###########################################################
###------  Import Relevant Functions from poaV2   ------###
###########################################################
# proc read_score_matrix*(filename: cstring; m: ptr ResidueScoreMatrix_T): cint {.importc: "read_score_matrix".}
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

proc free_lpo_sequence(sequence : ptr LPOSequence_T, please_free_holder : cint) {.importc: "free_lpo_sequence".}
proc free_return_result(lpo_result : LPOReturnResult_T ) {.importc: "free_return_result".}
# proc free_residue_matrix(lpo_matrix : ptr ResidueScoreMatrix_T) {.importc: "free_residue_matrix".}
# proc read_fasta*(seq_file: PFile;
#                  sequences: ptr ptr Sequence_T;
#                  do_switch_case: cint;
#                  comment: cstringArray): cint {.importc.}

# proc initialize_seqs_as_lpo*(nseq: cint;
#                              sequences: ptr Sequence_T;
#                              m: ptr ResidueScoreMatrix_T)

# proc buildup_progressive_lpo*(nseq : cint;
#                               # seqs : ptr ptr LPOSequence_T;
#                               seqs : ptr UncheckedArray[LPOSequence_T];
#                               score_matrix : ptr ResidueScoreMatrix_T;
#                               use_aggressive_fusion : cint;
#                               do_progressive: cint;
#                               score_file : cstring;
#                               scoring_function : proc (a1 : cint; a2 : cint; a3 : ptr LPOLetter_T; a4 : ptr LPOLetter_T; a5 : ptr ResidueScoreMatrix_T) : LPOScore_T;
#                               use_global_alignment : cint;
#                               preserve_sequence_order : cint) : ptr LPOSequence_T {.importc: "buildup_progressive_lpo".}

proc build_lpo_from_fasta*(seq_file : Pfile;
               score_matrix_filepath : cstring;
               use_global_alignment : cint;
               scoring_function : proc (a1 : cint; a2 : cint; a3 : ptr LPOLetter_T; a4 : ptr LPOLetter_T; a5 : ptr ResidueScoreMatrix_T) : LPOScore_T) : LPOReturnResult_T {.importc : "build_lpo_from_fasta".}

###########################################################
###---------  Utility functions to interface   ---------###
###########################################################

proc convertPOFormats*(lpo_return_result : LPOReturnResult_T) : POGraph = #TODO
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
  for i in 0..<num_reads:
    let name = source_seqs[i].name
    let length = source_seqs[i].length
    reads.add(Read(name:toString(toSeq(name)),length : uint32(length)))
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
  # for u in edges.keys:
  #   for v in edges[u]:
  #     assert (u,v) in weights
  # for (u,v) in  weights.keys:
  #   assert u in edges
  #   assert v in edges[u]
  free_return_result(lpo_return_result)
  # for i in 1..<num_reads:
  #   # echo &"freeing1 {(cast[int](addr seqs[num_reads - i])):#x}"
  #   free_lpo_sequence(addr seqs[num_reads - i],cint(0))
  # # echo &"freeing2 {(cast[int](addr seqs[0])):#x}"
  # free_lpo_sequence(addr seqs[0],cint(1))
  # free_residue_matrix(lpo_return_result.matrix)
  for i in 0..<num_reads:
    # echo paths[i]
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
  
  # for u in edges.keys:
  #   for v in edges[u]:
  #     assert (u,v) in weights
  # for (u,v) in weights.keys:
  #   assert u in edges
  #   assert v in edges[u]
    # echo c
# proc generateScoreMatrix() : ResidueScoreMatrix_T = #TODO
#   var x,y: array[32,ResidueScore_T]
#   x = [ResidueScore_T(16), 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 0]
#   y = [ResidueScore_T(16), 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 0]
#   return ResidueScoreMatrix_T(nsymbol : cint(16),
#                               symbol : ['A', 'T', 'G', 'C', 'U', 'S', 'W', 'R', 'Y', 'K', 'M', 'B', 'V', 'H', 'D', 'N'],
#                               score : [[ResidueScore_T( 0), ResidueScore_T(-8), ResidueScore_T(-8), ResidueScore_T(-8), ResidueScore_T(-8), ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T(-4),  ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-1),  ResidueScore_T(-2)],
#                                        [ResidueScore_T(-8), ResidueScore_T( 3), ResidueScore_T(-8), ResidueScore_T(-8), ResidueScore_T( 3), ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T( 1),  ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-1),  ResidueScore_T(-2)],
#                                        [ResidueScore_T(-8), ResidueScore_T(-8), ResidueScore_T( 3), ResidueScore_T(-8), ResidueScore_T(-8), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T( 1),  ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-4), ResidueScore_T(-1),  ResidueScore_T(-2)],
#                                        [ResidueScore_T(-8), ResidueScore_T(-8), ResidueScore_T(-8), ResidueScore_T( 3), ResidueScore_T(-8), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T(-4),  ResidueScore_T( 1), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-4),  ResidueScore_T(-2)],
#                                        [ResidueScore_T(-8), ResidueScore_T( 3), ResidueScore_T(-8), ResidueScore_T(-8), ResidueScore_T( 3), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-4),  ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-4),  ResidueScore_T(-4)],
#                                        [ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-4), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-2),  ResidueScore_T(-2), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-3), ResidueScore_T(-3),  ResidueScore_T(-1)],
#                                        [ResidueScore_T( 1), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-2),  ResidueScore_T(-2), ResidueScore_T(-3), ResidueScore_T(-3), ResidueScore_T(-1), ResidueScore_T(-1),  ResidueScore_T(-1)],
#                                        [ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-1), ResidueScore_T(-4), ResidueScore_T(-2),  ResidueScore_T(-2), ResidueScore_T(-3), ResidueScore_T(-1), ResidueScore_T(-3), ResidueScore_T(-1),  ResidueScore_T(-1)],
#                                        [ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-2),  ResidueScore_T(-2), ResidueScore_T(-1), ResidueScore_T(-3), ResidueScore_T(-1), ResidueScore_T(-3),  ResidueScore_T(-1)],
#                                        [ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-1),  ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-3), ResidueScore_T(-3), ResidueScore_T(-1),  ResidueScore_T(-1)],
#                                        [ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T( 1), ResidueScore_T(-4), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-4),  ResidueScore_T(-1), ResidueScore_T(-3), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-3),  ResidueScore_T(-1)],
#                                        [ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-3), ResidueScore_T(-3), ResidueScore_T(-1), ResidueScore_T(-1),  ResidueScore_T(-3), ResidueScore_T(-1), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-2),  ResidueScore_T(-1)],
#                                        [ResidueScore_T(-1), ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-3), ResidueScore_T(-1), ResidueScore_T(-3), ResidueScore_T(-3),  ResidueScore_T(-1), ResidueScore_T(-2), ResidueScore_T(-1), ResidueScore_T(-2), ResidueScore_T(-2),  ResidueScore_T(-1)],
#                                        [ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-4), ResidueScore_T(-3), ResidueScore_T(-1), ResidueScore_T(-3), ResidueScore_T(-1), ResidueScore_T(-3),  ResidueScore_T(-1), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-1), ResidueScore_T(-2),  ResidueScore_T(-1)],
#                                        [ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-4), ResidueScore_T(-4), ResidueScore_T(-3), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-3), ResidueScore_T(-1),  ResidueScore_T(-3), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-1),  ResidueScore_T(-1)],
#                                        [ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-2), ResidueScore_T(-4), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-1),  ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-1), ResidueScore_T(-1),  ResidueScore_T(-1)]],
#                               best_match :  [[cint( 0), cint(10), cint( 6), cint( 7), cint(13), cint(14), cint(12), cint(15), cint( 9), cint( 5), cint(11), cint( 8), cint( 1), cint(4), cint(3), cint(2)],
#                                              [cint( 4), cint( 1), cint( 9), cint( 6), cint( 8), cint(11), cint(13), cint(14), cint(15), cint(12), cint(10), cint( 5), cint( 7), cint(3), cint(2), cint(0)],
#                                              [cint( 2), cint( 5), cint( 7), cint( 9), cint(12), cint(14), cint(11), cint(15), cint( 8), cint( 6), cint(10), cint(13), cint( 4), cint(1), cint(3), cint(0)],
#                                              [cint( 3), cint( 8), cint( 5), cint(10), cint(13), cint(12), cint(11), cint(15), cint( 6), cint( 9), cint(14), cint( 7), cint( 0), cint(4), cint(2), cint(1)],
#                                              [cint( 4), cint( 1), cint(15), cint(12), cint(13), cint(14), cint( 5), cint( 6), cint( 7), cint( 8), cint( 9), cint(10), cint(11), cint(3), cint(2), cint(0)],
#                                              [cint( 2), cint( 3), cint(11), cint(15), cint(12), cint( 5), cint( 8), cint( 7), cint( 9), cint(10), cint(13), cint(14), cint( 1), cint(6), cint(4), cint(0)],
#                                              [cint( 0), cint( 1), cint(15), cint( 6), cint(13), cint(14), cint( 7), cint( 8), cint( 9), cint(10), cint(12), cint(11), cint( 5), cint(3), cint(4), cint(2)],
#                                              [cint( 2), cint( 0), cint(15), cint( 7), cint(12), cint(14), cint( 6), cint( 9), cint(10), cint( 5), cint(11), cint(13), cint( 4), cint(1), cint(8), cint(3)],
#                                              [cint( 3), cint( 1), cint(15), cint( 8), cint(11), cint(13), cint( 6), cint( 9), cint(10), cint( 5), cint(14), cint(12), cint( 4), cint(7), cint(0), cint(2)],
#                                              [cint( 2), cint( 1), cint(15), cint(11), cint( 9), cint(14), cint( 8), cint( 5), cint( 6), cint( 7), cint(12), cint(13), cint(10), cint(0), cint(4), cint(3)],
#                                              [cint( 3), cint( 0), cint(15), cint(10), cint(12), cint(13), cint( 6), cint( 7), cint( 8), cint( 5), cint(11), cint(14), cint( 1), cint(2), cint(4), cint(9)],
#                                              [cint(15), cint( 1), cint( 2), cint( 3), cint( 5), cint( 8), cint( 9), cint(11), cint(13), cint(12), cint(14), cint( 7), cint(10), cint(6), cint(0), cint(4)],
#                                              [cint( 0), cint( 2), cint( 3), cint( 5), cint( 7), cint(10), cint(12), cint(15), cint(13), cint(11), cint(14), cint( 8), cint( 9), cint(6), cint(1), cint(4)],
#                                              [cint(15), cint( 1), cint( 3), cint( 6), cint( 8), cint(10), cint(13), cint( 0), cint(14), cint(11), cint(12), cint( 7), cint( 9), cint(5), cint(4), cint(2)],
#                                              [cint( 0), cint( 1), cint( 2), cint( 6), cint( 7), cint( 9), cint(14), cint(15), cint(12), cint(11), cint(13), cint( 8), cint(10), cint(5), cint(4), cint(3)],
#                                              [cint(15), cint(10), cint(11), cint(12), cint(13), cint(14), cint( 5), cint( 6), cint( 7), cint( 8), cint( 9), cint( 1), cint( 2), cint(3), cint(0), cint(4)]],
#                               gap_penalty_set : [[ResidueScore_T(16), ResidueScore_T(6), ResidueScore_T(1)],
#                                                  [ResidueScore_T(16), ResidueScore_T(6), ResidueScore_T(1)]],
#                               trunc_gap_length : cint(20),
#                               decay_gap_length : cint(10),
#                               gap_penalty_x : addr x,
#                               gap_penalty_y : addr y,
#                               max_gap_length : cint(30))


# ###########################################################
# ###--------------  Call those functions   --------------###
# ###########################################################
# proc callPOAV2(fasta : seq[FastaRecord],matrix : = provideScoreMatrix()) : POGraph = #TODO
#   for record in fasta:
#     convertFastaRecords(record)
#   return convertPOFormats()



# var matrix : ResidueScoreMatrix_T
# var s : cstring = "../../poaV2/myNUC3.4.4.mat"
# echo read_score_matrix(s,addr matrix)
# var matrix2 = generateScoreMatrix()


# let seqs : UncheckedArray[LPOSequence_T] =
# var pointer_seqs = addr seqs


# var outfile : cstring = "../clusters/fasta/cluster_2217.fa"
# var temp = fopen(outfile,"r")
# let s : cstring = "../../poaV2/myNUC3.4.4.mat"
# var lpo_return_result = build_lpo_from_fasta(temp,
#                                              s,
#                                              use_global_alignment = cint(1),
#                                              matrix_scoring_function)
# var po = convertPOFormats(lpo_return_result)
