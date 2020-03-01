import strutils
import strformat
import sequtils
import times
import algorithm
import tables
import sets
import math
import heapqueue
import neo
import os
import hts

type
  FastaRecord = object
    read_id : string
    sequence : string

  Read = object
    name : string
    length : uint32
    path : seq[uint32]
    corrected_path : seq[uint32]
    sequence : string
  
  Node = object
    nts : string
    supporting_reads : seq[uint32]
    align_ring_partner : int32 #No partner if -1
    visited : bool
    indegree : uint16
    log_likelihood : float
    path : seq[uint32]
    corrected_path : seq[uint32]
    start_node_flag : bool
    end_node_flag : bool
  
  POGraph = object of RootObj
    nodes : seq[Node]
    reads : seq[Read]
    # edges : Table[uint32,HashSet[uint32]]
    edges : Table[uint32,seq[uint32]]
    weights : Table[(uint32,uint32),uint32]
    og_nodes : uint32
  
  TrimmedPOGraph = object of POGraph
    log_likelihoods : Table[(uint32,uint32),float64]
    node_indexes : Table[(char,uint32), uint32]
    source_nodes : seq[uint32]
    end_nodes : seq[uint32]
    deleted_nodes : HashSet[uint32]
    illumina_branches : HashSet[seq[uint32]]
  
  StringtieLikeDataStructure = object
    nodes : seq[Node]
    node_indexes : Table[(char,uint32), uint32]
    fwd_edges : Table[uint32,seq[uint32]]
    rev_edges : Table[uint32,seq[uint32]]
    weights : Table[(uint32,uint32),seq[string]]
    node_support_list : seq[uint32]

  SplitProfileHMM = object
    nodes : seq[Node]
    num_source_nodes : uint32
    source_node_order : Table[uint32,uint32]
    num_match_states : uint32
    num_delete_states : uint32
    num_insert_states : uint32
    num_phmm_nodes : uint32
    sink_nodes : seq[uint32]
    start_node : uint32
    end_node : uint32
    end_nodes : HashSet[uint32]
    in_nodes : Table[uint32,seq[uint32]]
    # in_node_array : seq[uint32]
    # in_node_array_bounds : seq[uint32]
    in_node_array : Vector[uint32]
    in_node_array_bounds : Vector[uint32]
    transitions : Matrix[float64]
    emissions : Matrix[float64]
    windows : Table[uint32,HashSet[uint32]]

const NT_TO_IDX = [
0'u8, 1'u8, 2'u8, 3'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 0'u8, 4'u8, 1'u8,  4'u8, 4'u8, 4'u8, 2'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  3'u8, 3'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 0'u8, 4'u8, 1'u8,  4'u8, 4'u8, 4'u8, 2'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  3'u8, 3'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,
4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8,  4'u8, 4'u8, 4'u8, 4'u8
]
const IDX_TO_NT = ['A','G','C','T']

proc collapseLinearStretches( po2 : TrimmedPOGraph) : TrimmedPOGraph =
  var po = po2
  var to_delete : HashSet[uint32]
  var new_end_nodes : seq[uint32]
  # for idx in sorted(toSeq(po.node_indexes.keys)):
  #   echo "idx, ", idx, " ", po.node_indexes[idx]
  for i in 0..<po.nodes.len:
    po.nodes[i].visited = false
    # if po.nodes[i].end_node_flag:
    #   echo "flag ", i
  for i in 0..<po.nodes.len:
    if po.nodes[i].visited:
      continue
    let node_idx = po.node_indexes[('f',uint32(i))]
    if po.edges[node_idx].len == 1 and not po.nodes[i].end_node_flag:
      var next_node_idx = po.edges[node_idx][0]
      var back_next_node_idx = po.node_indexes[('b',next_node_idx)]
      while po.edges[node_idx].len == 1 and po.nodes[back_next_node_idx].indegree == 1 and 
            not po.nodes[back_next_node_idx].start_node_flag and
            not po.nodes[back_next_node_idx].end_node_flag:
        po.nodes[i].nts = po.nodes[i].nts & po.nodes[back_next_node_idx].nts
        po.nodes[back_next_node_idx].visited = true
        to_delete.incl(back_next_node_idx)
        po.edges[node_idx] = po.edges[next_node_idx]
        po.edges.del(next_node_idx)
        po.log_likelihoods.del((node_idx,next_node_idx))
        po.weights.del((node_idx,next_node_idx))
        if po.edges[node_idx].len == 1:
          let k = po.edges[node_idx][0]
          po.log_likelihoods[(node_idx,k)] = po.log_likelihoods[(next_node_idx,k)]
          po.log_likelihoods.del((next_node_idx,k))
          po.weights[(node_idx,k)] = po.weights[(next_node_idx,k)]
          po.weights.del((next_node_idx,k))
          next_node_idx = po.edges[node_idx][0]
          back_next_node_idx = po.node_indexes[('b',next_node_idx)]
        elif po.edges[node_idx].len > 1:
          for k in po.edges[node_idx]:
            po.log_likelihoods[(node_idx,k)] = po.log_likelihoods[(next_node_idx,k)]
            po.log_likelihoods.del((next_node_idx,k))
            po.weights[(node_idx,k)] = po.weights[(next_node_idx,k)]
            po.weights.del((node_idx,k))
          break
      if po.nodes[back_next_node_idx].end_node_flag:
        if po.edges[node_idx].len == 1 and po.nodes[back_next_node_idx].indegree == 1 and not po.nodes[back_next_node_idx].start_node_flag:
          po.nodes[i].nts = po.nodes[i].nts & po.nodes[back_next_node_idx].nts
          po.nodes[back_next_node_idx].visited = true
          to_delete.incl(back_next_node_idx)
          po.edges[node_idx] = po.edges[next_node_idx]
          po.edges.del(next_node_idx)
          po.log_likelihoods.del((node_idx,next_node_idx))
          po.weights.del((node_idx,next_node_idx))
          if po.edges[node_idx].len > 1:
            for k in po.edges[node_idx]:
              po.log_likelihoods[(node_idx,k)] = po.log_likelihoods[(next_node_idx,k)]
              po.log_likelihoods.del((next_node_idx,k))
              po.weights[(node_idx,k)] = po.weights[(next_node_idx,k)]
              po.weights.del((next_node_idx,k))
          new_end_nodes.add(node_idx)
          # echo "here ", node_idx
        elif po.edges[node_idx].len != 0:
          new_end_nodes.add(node_idx)
          # echo "here3 ", node_idx
        else:
          new_end_nodes.add(next_node_idx)
          # echo "here2 ", next_node_idx
    else:
      if po.nodes[i].end_node_flag:
        new_end_nodes.add(node_idx)
  var new_nodes : seq[Node]
  var new_node_indexes : Table[(char,uint32),uint32]
  var counter = 0'u32
  for i in 0'u32..<uint32(po.nodes.len):
    if i in to_delete:
      continue
    new_nodes.add(po.nodes[i])
    new_node_indexes[('f',counter)] = po.node_indexes[('f',i)]
    new_node_indexes[('b',po.node_indexes[('f',i)])] = counter
    counter += 1'u32
  return TrimmedPOGraph(nodes : new_nodes,
                        node_indexes : new_node_indexes,
                        edges: po.edges,
                        weights : po.weights,
                        reads : po.reads,
                        log_likelihoods : po.log_likelihoods,
                        source_nodes : po.source_nodes,
                        end_nodes : new_end_nodes)

proc jsOutput( po:TrimmedPOGraph, highlight_path1,highlight_path2 : seq[uint32] = @[],collapse : bool = false) : seq[string] = 
  var path : seq[uint32]
  if highlight_path1.len != 0:
    path = highlight_path1
  else:
    path = po.reads[0].corrected_path
  
  var highlight_nodes1,highlight_nodes2 : HashSet[uint32]
  var highlight_edges1,highlight_edges2 : HashSet[(uint32,uint32)]
  for i,node in highlight_path1:
    highlight_nodes1.incl(node)
    if i != 0:
      highlight_edges1.incl((highlight_path1[i-1],node))
  for i,node in highlight_path2:
    highlight_nodes2.incl(node)
    if i != 0:
      highlight_edges2.incl((highlight_path2[i-1],node))

  var pathTable : Table[uint32,int]
  for i, nodeID in path:
    pathTable[nodeID] = i * 150
  var lines = @["var nodes = ["]
  var count = 0
  var last_j_in_path_table = 999999
  let f = '{'
  let e = '}'
  for i in pathTable.keys():
    if int(i) < last_j_in_path_table:
      last_j_in_path_table = int(i)
  for i in 0..<po.nodes.len:
    let j = uint32(i)
    let node = po.nodes[po.node_indexes[('b',uint32(i))]]
    var line : seq[string]
    if collapse:
      line = @[&"    {f}id: {j}, label: \"{j}\""]
    else:
      line = @[&"    {f}id: {j}, label: \"{j}, {node.nts}\""]
    var highlight : string
    if j in highlight_nodes1:
      if j in highlight_nodes2:
        highlight = "color: \'purple\', "
      else:
        highlight = "color: \'red\', "
    elif j in highlight_nodes2:
      highlight = "color: \'blue\', "
    else:
      highlight = "color: \'#D3D3D3\', "
    if (j in pathTable):
      if((count %% 5) == 0):
        line.add(&", allowedToMoveX: true, x: {pathTable[j]} , y: 0, {highlight} allowedToMoveY: true{e},")
        last_j_in_path_table = int(j)
      else:
        line.add(&", allowedToMoveX: true, x: {pathTable[uint32(last_j_in_path_table)]} , y: 0, {highlight} allowedToMoveY: true{e},")
      count += 1
    else:
      line.add(&", allowedToMoveX: true, x: {pathTable[uint32(last_j_in_path_table)]} , y: 0, {highlight} allowedToMoveY: true{e},")
    lines.add(line.join(""))
  lines[^1] = lines[^1][0..^1]
  lines.add("];")
  lines.add(" ")
  lines.add("var edges = [")

  var written_set : HashSet[(uint32,uint32)]
  for i,node in po.nodes:
    let k = po.node_indexes[('f',uint32(i))]
    for edge in po.edges[k]:
      var highlight : string
      if not((k,edge) in written_set):
        if (k,edge) in highlight_edges1:
          if (k,edge) in highlight_edges2:
            highlight = ", color: \'purple\'"
          else:
            highlight = ", color: \'red\'"
        elif (k,edge) in  highlight_edges2:
          highlight = ", color: \'blue\'"
        else:
          highlight = ", color: \'black\'"
        lines.add(&"    {f}from: {k}, to: {edge}{highlight}, label: '{po.weights[(k,edge)]}',     font: {f}align: 'middle'{e}, value: {po.weights[(k,edge)]}, arrows:\'to\'{e},")
        written_set.incl((k,edge))
  lines[^1] = lines[^1][0..^1]
  lines.add("];\n")
  return lines

proc htmlOutput( po:TrimmedPOGraph, outfile : File, highlight_path1,highlight_path2 : seq[uint32] = @[],collapse = false) = 
  let header = "<!doctype html>\n<html>\n<head>\n<title>POA Graph Alignment</title>\n<script type=\"text/javascript\" src=\"https://visjs.github.io/vis-network/standalone/umd/vis-network.min.js\"></script>\n</head>\n<body>\n<div id=\"mynetwork\"></div>\n<script type=\"text/javascript\">\n// create a network\n"
  outfile.write(header)
  var lines :seq[string]
  if collapse:
    lines = jsOutput(collapseLinearStretches(po),highlight_path1,highlight_path2,collapse)
  else:
    lines = jsOutput(po,highlight_path1,highlight_path2)
  outfile.write(lines.join("\n"))
  let footer = "var container = document.getElementById('mynetwork');\nvar data= {\nnodes: nodes,\nedges: edges,\n};\nvar options = {\nwidth: '100%',\nheight: '800px'\n};\nvar network = new vis.Network(container, data, options);\n</script>\n</body>\n</html>\n"
  outfile.write(footer)
  outfile.close()

proc writeCorrectedReads( po : TrimmedPOGraph,outfile : File) = 
  for read in po.reads:
    var sequence : seq[string]
    for idx in read.corrected_path:
      sequence.add(po.nodes[po.node_indexes[('b',idx)]].nts)
    outfile.write(">",read.name,"\n")
    outfile.writeLine(sequence.join())

proc writeCorrectedReads( records : seq[FastaRecord],outfile : File) = 
  for record in records:
    outfile.write(">",record.read_id,"\n")
    outfile.writeLine(record.sequence)

proc initPOGraph( file : File) : POGraph = 
  for _ in 0..2:
    discard file.readLine()
  let num_nodes = uint32(parseUInt(file.readLine().split('=')[1]))
  let num_reads = uint32(parseUInt(file.readLine().split('=')[1]))
  var reads : seq[Read]
  var nodes : seq[Node]
  var source_node_indices : seq[uint32]
  var end_node_indices : seq[uint32]
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  for _ in 0..<num_reads:
    let name = file.readLine().strip().split('=')[1]
    let length = uint32(parseUInt(file.readLine().strip().split('=')[1].split(' ')[0]))
    reads.add(Read(name:name,length:length))
  for i in 0..<num_nodes:
    let tmp = file.readLine().strip().split(':')
    let nt = tmp[0]
    var info = tmp[1]
    let a_pos = info.find("A")
    var a_pair : int32
    if a_pos != -1:
      a_pair = int32(parseUInt(info[a_pos+1..^1]))
      info = info[0..<a_pos]
    else:
      a_pair = -1
    nodes.add(Node(nts : nt,
                   align_ring_partner : a_pair,
                   visited : false,
                   end_node_flag : false,
                   start_node_flag : false))
    let s_split = info.split('S')
    info = s_split[0]
    var seq_idxs : seq[uint32]
    for s in s_split[1..^1]:
      seq_idxs.add(uint32(parseUInt(s)))
    for seq_idx in seq_idxs:
      if reads[seq_idx].path.len == 0:
        nodes[i].start_node_flag = true
        source_node_indices.add(uint32(i))
      reads[seq_idx].path.add(uint32(i))
      reads[seq_idx].sequence.add(nt)
      nodes[i].supporting_reads.add(seq_idx)
    if info != "":
      let l_split = info.split('L')
      for l in l_split[1..^1]:
        let in_node = uint32(parseUInt(l))
        if in_node in edges:
          edges[in_node].add(uint32(i))
        else:
          edges[in_node] = @[uint32(i)]
        weights[(uint32(in_node),uint32(i))] = 0'u32
  for read in reads:
    # nodes[read.path[^1]].end_node_flag = true
    end_node_indices.add(read.path[^1])
  for i in 0'u32..<uint32(nodes.len):
    if not (i in edges):
      edges[i] = @[]
    else:
      for j in edges[i]:
        nodes[j].indegree += 1
  for read in reads:
    for i in 0..<read.path.len - 1:
      weights[((read.path[i],read.path[i+1]))] += 1'u32
  # echo weights
  return POGraph(nodes:nodes,reads:reads,edges:edges,og_nodes:uint32(nodes.len),weights:weights)

proc parseFasta(file : File) : seq[FastaRecord] = 
  var records : seq[FastaRecord]
  var read_id : string
  var sequence : string
  var count = 0
  while true:
    try:
      let line = file.readLine()
      if line[0] == '>':
        if count != 0:
          records.add(FastaRecord( read_id : read_id, sequence : sequence))
          sequence = ""
        count += 1
        read_id = line.strip(leading=true,trailing=false,chars = {'>'}).strip(leading=false,trailing=true,chars = {'\n'})
      else:
        sequence = sequence & line.strip(leading=false,trailing=true,chars = {'\n'})
    except EOFError:
      break
  return records

proc writePOGraph( po : TrimmedPOGraph, outfile : File,graphname : string ="default") = 
  outfile.write("VERSION=UNTITLEDCORRECTIONALGORITHM.0.1\n")
  outfile.write(&"NAME={graphname}\n")
  outfile.write("TITLE=untitled\n")
  outfile.write(&"LENGTH={po.nodes.len}\n")
  outfile.write(&"SOURCECOUNT={po.reads.len}\n")
  for read in po.reads:
    outfile.write(&"SOURCENAME={read.name}\n")
    outfile.write(&"SOURCEINFO={read.corrected_path.len} 0 1 -1 untitled\n")
  var reverse_edges : Table[uint32,seq[uint32]]
  for u in po.edges.keys():
    for v in po.edges[u]:
      if v in reverse_edges:
        reverse_edges[v].add(u)
      else:
        reverse_edges[v] = @[u]
  var supporting_reads_table : Table[uint32,seq[uint32]]
  for i,read in po.reads:
    # echo read.corrected_path
    for node in read.corrected_path:
      if not(node in supporting_reads_table):
        supporting_reads_table[node] = @[uint32(i)]
      else:
        supporting_reads_table[node].add(uint32(i))
  for i,node in po.nodes:
    var edges : seq[string]
    if po.node_indexes[('f',uint32(i))] in reverse_edges:
      for u in reverse_edges[po.node_indexes[('f',uint32(i))]]:
        edges.add(&"L{po.node_indexes[('b',u)]}")
    var supporting_reads : seq[string]
    for read_num in supporting_reads_table[po.node_indexes[('f',uint32(i))]]:
      supporting_reads.add(&"S{read_num}")
    outfile.write(&"{node.nts}:{edges.join()}{supporting_reads.join()}\n")
  outfile.close()

proc constructNewGraph( po : ptr POGraph,reads:seq[Read],previously_deleted_nodes:HashSet[uint32] = initHashSet[uint32]() ) : TrimmedPOGraph = 
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  var u_set,v_set : HashSet[uint32]
  for read in reads:
    for i in 0..<read.corrected_path.len - 1:
      if (read.corrected_path[i],read.corrected_path[i+1]) in weights:
        weights[(read.corrected_path[i],read.corrected_path[i+1])] += 1'u32
      else:
        weights[(read.corrected_path[i],read.corrected_path[i+1])] = 1'u32
        if read.corrected_path[i] in edges:
          edges[read.corrected_path[i]].add(read.corrected_path[i+1])
        else:
          edges[read.corrected_path[i]] = @[read.corrected_path[i+1]]
        u_set.incl(read.corrected_path[i])
        v_set.incl(read.corrected_path[i+1])
    # for i in 0..<read.corrected_path.len:
    #   echo read.corrected_path[i]
  let source_nodes = u_set - v_set
  let sink_nodes = v_set - u_set
  for sink_node in sink_nodes:
    edges[sink_node] = @[]
  var nodes : seq[Node]
  var node_indexes = sorted(toSeq(u_set + v_set))
  var node_indexes2 : Table[(char,uint32),uint32]
  var log_likelihoods : Table[(uint32,uint32),float64]
  for i,j in node_indexes:
    nodes.add(po[].nodes[j])
    node_indexes2[('f',uint32(i))] = j
    node_indexes2[('b',j)] = uint32(i)
    if j in edges:
      var total_read_count = 0'u32
      for k in edges[j]:
        total_read_count += weights[(j,k)]
      let tmp = float(total_read_count)
      for k in edges[j]:
        log_likelihoods[(j,k)] = - ln(float(weights[(j,k)]) / tmp)
  var deleted_nodes = previously_deleted_nodes
  for i in 0..<po[].og_nodes:
    if uint32(i) notin u_set and uint32(i) notin v_set:
      deleted_nodes.incl(uint32(i))
  for i in 0..<nodes.len:
    nodes[i].indegree = 0'u16
  for j in edges.keys:
    for k in edges[j]:
      nodes[node_indexes2[('b',k)]].indegree += 1'u16
  # for j in sorted(toSeq(edges.keys)):
  #   for k in sorted(toSeq(edges[j])):
  #     echo j, " ", k, " ", weights[(j,k)]
  # for key in sorted(toSeq(node_indexes2.keys)):
  #   echo "idx: ", key, " ", node_indexes2[key]
  return TrimmedPOGraph(nodes : nodes,
                        reads : reads,
                        edges : edges,
                        og_nodes : po[].og_nodes,
                        weights : weights,
                        log_likelihoods : log_likelihoods,
                        source_nodes : toSeq(source_nodes),
                        end_nodes : toSeq(sink_nodes),
                        node_indexes : node_indexes2,
                        deleted_nodes : deleted_nodes)

proc constructNewGraph2( po : ptr TrimmedPOGraph ) : TrimmedPOGraph = 
  let reads = po[].reads
  let previously_deleted_nodes = po.deleted_nodes
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  var u_set,v_set : HashSet[uint32]
  for read in reads:
    for i in 0..<read.corrected_path.len - 1:
      if (read.corrected_path[i],read.corrected_path[i+1]) in weights:
        weights[(read.corrected_path[i],read.corrected_path[i+1])] += 1'u32
      else:
        weights[(read.corrected_path[i],read.corrected_path[i+1])] = 1'u32
        if read.corrected_path[i] in edges:
          edges[read.corrected_path[i]].add(read.corrected_path[i+1])
        else:
          edges[read.corrected_path[i]] = @[read.corrected_path[i+1]]
        u_set.incl(read.corrected_path[i])
        v_set.incl(read.corrected_path[i+1])
    # for i in 0..<read.corrected_path.len:
    #   echo read.corrected_path[i]
  let source_nodes = u_set - v_set
  let sink_nodes = v_set - u_set
  for sink_node in sink_nodes:
    edges[sink_node] = @[]
  var nodes : seq[Node]
  var node_indexes = sorted(toSeq(u_set + v_set))
  var node_indexes2 : Table[(char,uint32),uint32]
  var log_likelihoods : Table[(uint32,uint32),float64]
  for i,j in node_indexes:
    nodes.add(po[].nodes[po.node_indexes[('b',j)]])
    node_indexes2[('f',uint32(i))] = j
    node_indexes2[('b',j)] = uint32(i)
    if j in edges:
      var total_read_count = 0'u32
      for k in edges[j]:
        total_read_count += weights[(j,k)]
      let tmp = float(total_read_count)
      for k in edges[j]:
        log_likelihoods[(j,k)] = - ln(float(weights[(j,k)]) / tmp)
  var deleted_nodes = previously_deleted_nodes
  for i in 0..<po[].og_nodes:
    if uint32(i) notin u_set and uint32(i) notin v_set:
      deleted_nodes.incl(uint32(i))
  for i in 0..<nodes.len:
    nodes[i].indegree = 0'u16
  for j in edges.keys:
    for k in edges[j]:
      nodes[node_indexes2[('b',k)]].indegree += 1'u16
  # for j in sorted(toSeq(edges.keys)):
  #   for k in sorted(toSeq(edges[j])):
  #     echo j, " ", k, " ", weights[(j,k)]
  # for key in sorted(toSeq(node_indexes2.keys)):
  #   echo "idx: ", key, " ", node_indexes2[key]
  return TrimmedPOGraph(nodes : nodes,
                        reads : reads,
                        edges : edges,
                        og_nodes : po[].og_nodes,
                        weights : weights,
                        log_likelihoods : log_likelihoods,
                        source_nodes : toSeq(source_nodes),
                        end_nodes : toSeq(sink_nodes),
                        node_indexes : node_indexes2,
                        deleted_nodes : deleted_nodes)

proc getStringtieLikeDataStructures( po : ptr POGraph,reads:seq[Read]) : StringtieLikeDataStructure =
  var fwd_edges,rev_edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),seq[string]]
  var u_set,v_set : HashSet[uint32]
  var node_support_counts : CountTable[uint32]
  for read in reads:
    for i in 0..<read.corrected_path.len-1:
      if (read.corrected_path[i],read.corrected_path[i+1]) in weights:
        weights[(read.corrected_path[i],read.corrected_path[i+1])].add(read.name)
      else:
        weights[(read.corrected_path[i],read.corrected_path[i+1])] = @[read.name]
        if read.corrected_path[i] in fwd_edges:
          fwd_edges[read.corrected_path[i]].add(read.corrected_path[i+1])
        else:
          fwd_edges[read.corrected_path[i]] = @[read.corrected_path[i+1]]
        if read.corrected_path[i+1] in rev_edges:
          rev_edges[read.corrected_path[i+1]].add(read.corrected_path[i])
        else:
          rev_edges[read.corrected_path[i+1]] = @[read.corrected_path[i]]
        u_set.incl(read.corrected_path[i])
        v_set.incl(read.corrected_path[i+1])
      if read.corrected_path[i] in node_support_counts:
        node_support_counts.inc(read.corrected_path[i])
      else:
        node_support_counts[read.corrected_path[i]] = 1
  let source_nodes = u_set - v_set
  let sink_nodes = v_set - u_set
  for sink_node in sink_nodes:
    fwd_edges[sink_node] = @[]
  for source_node in source_nodes:
    rev_edges[source_node] = @[]
  var nodes : seq[Node]
  let node_indexes = sorted(toSeq(u_set + v_set))
  var node_indexes2 : Table[(char,uint32),uint32]
  var node_support_list : seq[uint32]
  for i,j in node_indexes:
    nodes.add(po[].nodes[j])
    node_indexes2[('f',uint32(i))] = j
    node_indexes2[('b',j)] = uint32(i)
    node_support_list.add(uint32(node_support_counts[j]))
  return StringtieLikeDataStructure(nodes : nodes,
                                node_indexes : node_indexes2,
                                fwd_edges : fwd_edges,
                                rev_edges : rev_edges,
                                weights : weights,
                                node_support_list: node_support_list)

proc maxIdx[T]( s : seq[T]) : int = 
  var current_max = T.low
  var current_max_idx = 0
  for i,j in  s:
    if j > current_max:
      current_max = j
      current_max_idx = i
  return current_max_idx

proc orderedIdx[T]( s : seq[T], i : int) : int = 
  var list : seq[(T,int)]
  for k,j in s:
    list.add((j,k))
  # var current_max = T.low
  # var current_max_idx = 0
  let tmp = sorted(list,Descending)
  return tmp[i][1]

proc alignPaths( x : seq[uint32], y : seq[uint32]) : seq[uint8] = 
  ### Takes in two paths and aligns them, leveraging the fact that the nodes are labeled in topological sort order to speed alignment
  ### x - first path, query for the purposes of defining insertion / deletion
  ### y - second path, reference for the purposes of defining insertion / deletion
  var alignment : seq[uint8]
  var i = 0
  var j = 0
  # echo x
  # echo y
  while true:
    if x[i] == y[j]:
      alignment.add(0'u8)
      i += 1
      j += 1
    elif x[i] < y[j]:
      alignment.add(1'u8)
      i += 1
    else: # x[i] > y[j]
      alignment.add(2'u8)
      j += 1
    if i == x.len or j == y.len:
      break
  for _ in i..<x.len:
    alignment.add(1'u8)
  for _ in j..<y.len:
    alignment.add(2'u8)
  return alignment

proc compareQopTuples( x,y:(char,int,int)) : int = 
  if x[1] < y[1]:
    return -1
  elif x[1] > y[1]:
    return 1
  if x[0] < y[0]:
    return -1
  elif x[0] > y[0]:
    return 1
  if x[2] < y[2]:
    return -1
  elif x[2] > y[2]:
    return 1
  return 0

proc correctPaths( q_path : var seq[uint32],
                  r_path : seq[uint32],
                  nodes: ptr seq[Node],
                  node_indexes : Table[(char,uint32),uint32],
                  psi:uint16,
                  correct_ends : bool = false,
                  debug : bool = false) : (bool,seq[uint32]) =
  ### Takes in two paths, a query path and a reference path, and a number psi that indicates the min indel size needed to retain that indel, indels in the query relative to the reference smaller than that number are 'corrected', to reflect the reference path. Also takes in a bool indicating whether or not to correct the ends of alignments, a reference to nodes, which get modified to reflect new start and end nodes based on correction.
  # echo "here"
  for i in 0..<q_path.len - 1:
    assert q_path[i] < q_path[i+1]
  for i in 0..<r_path.len - 1:
    assert r_path[i] < r_path[i+1]
  let alignment = alignPaths(q_path,r_path)
  var q_idx = 0
  var r_idx = 0
  var continuous_ins = 0
  var continuous_del = 0
  var diff_flag = false
  var q_ops : seq[(char,int,int)]
  var start_flag = true
  var ignore_deletion_flag = false
  for op in alignment:
    if op == 0'u8:
      if debug:
        echo q_path[q_idx], "\t", r_path[r_idx]
      if continuous_ins >= int(psi):
        diff_flag = true
        ignore_deletion_flag = true
      elif continuous_ins > 0:
        # echo continuous_ins
        for i in q_idx - continuous_ins..<q_idx:
          # echo "here"
          q_ops.add(('i',i,0))
      if continuous_del >= int(psi):
        diff_flag = true
      elif continuous_del > 0 and not ignore_deletion_flag:
        for i in r_idx - continuous_del..<r_idx:
          q_ops.add(('d',q_idx-continuous_ins,i))
      ignore_deletion_flag = false
      start_flag = false
      q_idx += 1
      r_idx += 1
      continuous_ins = 0
      continuous_del = 0
    elif op == 1'u8:
      if debug:
        echo q_path[q_idx], "\t----"
      if correct_ends or not start_flag:
        continuous_ins += 1
      q_idx += 1
    else: #op == 2
      if debug:
        echo "----\t", r_path[r_idx]
      if not start_flag:
        continuous_del += 1
      r_idx += 1
  if correct_ends:
    if continuous_ins >= int(psi):
      diff_flag = true
    elif continuous_ins > 0:
      for i in q_idx - continuous_ins..<q_idx:
        q_ops.add(('i',i,0))
  assert q_idx == q_path.len
  assert r_idx == r_path.len
  let before = q_path
  for (op,q_idx,r_idx) in sorted(q_ops,compareQopTuples,Descending):
    if op == 'i':
      if q_idx != len(q_path) - 1 and q_idx != 0:
        assert q_path[q_idx - 1] < q_path[q_idx + 1]
      # let before2 = q_path
      # echo q_idx
      # echo q_path.len
      q_path.delete(q_idx)
      # echo before2 == q_path
    else: # op == 'd'
      if q_idx != 0:
        assert q_path[q_idx - 1] < r_path[r_idx] and r_path[r_idx] < q_path[q_idx]
      q_path.insert(@[r_path[r_idx]],q_idx)
  # if correct_ends:
  #   nodes[][node_indexes[('b',q_path[0])]].start_node_flag = true
  #   nodes[][node_indexes[('b',q_path[^1])]].end_node_flag = true
  if debug:
    # echo q_path == before
    echo q_path.len
  return (diff_flag,q_path)

proc walkHeaviestPaths( po : ptr POGraph,psi = 15) : seq[seq[uint32]] = 
  var representative_paths : seq[seq[uint32]]
  var remaining_reads = toSeq(0..<po[].reads.len)
  while remaining_reads.len != 0:
    # echo "here"
    var reads : seq[Read] = @[]
    for i in remaining_reads:
      reads.add(po[].reads[i])
    
    let st = getStringtieLikeDataStructures(po,reads)

    let max_node_idx = maxIdx(st.node_support_list)
    let fwd_max_node_idx = st.node_indexes[('f',uint32(max_node_idx))]
    var start_idx,end_idx : int
    if st.fwd_edges[fwd_max_node_idx].len == 0:
      end_idx = -1
    else:
      end_idx = max_node_idx
    if st.rev_edges[fwd_max_node_idx].len == 0:
      start_idx = -1
    else:
      start_idx = max_node_idx
    
    var read_idxs = newSeqWith(po[].reads.len, [0,0])
    var excluded_reads : HashSet[string]
    var all_reads : HashSet[string]
    var fwd_end_reads : HashSet[string]
    var rev_end_reads : HashSet[string]
    var read_id_to_idx : Table[string,int]
    var read_idx_to_id : Table[int,string]
    for i,read in reads:
      read_id_to_idx[read.name] = i
      read_idx_to_id[i] = read.name
      all_reads.incl(read.name)
      let j = binarySearch(read.corrected_path, fwd_max_node_idx)
      read_idxs[i] = [j,j] #If not found it'll be [-1,-1], else it'll be the index of the matched node
      if j == -1:
        excluded_reads.incl(read.name)
    var representative_path = @[fwd_max_node_idx]
    while start_idx != -1 or end_idx != -1:
      # echo "here2"
      # echo "start - ", start_idx
      # echo "end   - ", end_idx
      var next_node_id : uint32
      if start_idx == -1:
        # Reached end of movement in reverse direction, look in forward direction
        let fwd_end_idx = st.node_indexes[('f',uint32(end_idx))]
        var fwd_weights : seq[uint32]
        for fwd_node in st.fwd_edges[fwd_end_idx]:
          var weight = 0'u32
          for read_id in st.weights[(fwd_end_idx,fwd_node)]:
            if (read_id notin excluded_reads) and (read_id notin fwd_end_reads):
              weight += 1'u32
          fwd_weights.add(weight)
        # echo fwd_weights
        if max(fwd_weights) == 0 and fwd_end_reads.len != 0:
          fwd_weights = @[]
          for fwd_node in st.fwd_edges[fwd_end_idx]:
            var weight = 0'u32
            for read_id in st.weights[(fwd_end_idx,fwd_node)]:
              if read_id notin excluded_reads:
                weight += 1'u32
            fwd_weights.add(weight)
          # echo fwd_weights
          if max(fwd_weights) == 0:
            break
        # echo "current - ", fwd_end_idx
        # echo "options - ", st.fwd_edges[fwd_end_idx]
        next_node_id = st.fwd_edges[fwd_end_idx][maxIdx(fwd_weights)]
        # echo "next - ", next_node_id
        for i in 0..<reads.len:
          let read_idx = read_idxs[i][1]
          if read_idx == -1:
            continue
          if read_idx + 1 != len(reads[i].corrected_path):
            if reads[i].corrected_path[read_idx + 1] == next_node_id:
              read_idxs[i][1] += 1
              if (reads[i].corrected_path.len - read_idxs[i][1]) < psi:
                fwd_end_reads.incl(reads[i].name)
                read_idxs[i][1] = -1
            else:
              read_idxs[i][1] = -1
              if true:
              # if reads[i].name notin fwd_end_reads:
                excluded_reads.incl(reads[i].name)
          else:
            read_idxs[i][1] = -1
        end_idx = int(st.node_indexes[('b',next_node_id)])
        # echo "back node", end_idx
        # echo "end - ", end_idx
        if st.fwd_edges[next_node_id].len == 0:
          end_idx = -1
      
      elif end_idx == -1:
        # Reached end of movement in forward direction, look in reverse direction
        let fwd_start_idx = st.node_indexes[('f',uint32(start_idx))]
        var rev_weights : seq[uint32]
        for rev_node in st.rev_edges[fwd_start_idx]:
          var weight = 0'u32
          for read_id in st.weights[(rev_node,fwd_start_idx)]:
            if (read_id notin excluded_reads) and (read_id notin rev_end_reads):
              weight += 1'u32
          rev_weights.add(weight)
        if max(rev_weights) == 0 and rev_end_reads.len != 0:
          rev_weights = @[]
          for rev_node in st.rev_edges[fwd_start_idx]:
            var weight = 0'u32
            for read_id in st.weights[(rev_node,fwd_start_idx)]:
              if (read_id notin excluded_reads):
                weight += 1'u32
            rev_weights.add(weight)
          if max(rev_weights) == 0:
            break
        next_node_id = st.rev_edges[fwd_start_idx][maxIdx(rev_weights)]
        for i in 0..<reads.len:
          let read_idx = read_idxs[i][0]
          if read_idx == -1:
            continue
          if read_idx != 0:
            if reads[i].corrected_path[read_idx - 1] == next_node_id:
              read_idxs[i][0] = read_idx - 1
              if read_idxs[i][0] < psi:
                rev_end_reads.incl(reads[i].name)
                read_idxs[i][0] = -1
            else:
              read_idxs[i][0] = -1
              if true:
              # if reads[i].name notin rev_end_reads:
                excluded_reads.incl(reads[i].name)
          else:
            read_idxs[i][0] = -1
        start_idx = int(st.node_indexes[('b',next_node_id)])
        # echo "start - ", start_idx
        if st.rev_edges[next_node_id].len == 0:
          start_idx = -1
      
      else:
        # Decide whether to move forward of reverse based on max weight in next step on the path
        let fwd_start_idx = st.node_indexes[('f',uint32(start_idx))]
        let fwd_end_idx = st.node_indexes[('f',uint32(end_idx))]

        var rev_weights,fwd_weights : seq[uint32]
        for rev_node in st.rev_edges[fwd_start_idx]:
          var weight = 0'u32
          for read_id in st.weights[(rev_node,fwd_start_idx)]:
            if (read_id notin excluded_reads) and (read_id notin rev_end_reads):
              weight += 1'u32
          rev_weights.add(weight)
        
        for fwd_node in st.fwd_edges[fwd_end_idx]:
          var weight = 0'u32
          for read_id in st.weights[(fwd_end_idx,fwd_node)]:
            if (read_id notin excluded_reads) and (read_id notin fwd_end_reads):
              weight += 1'u32
          fwd_weights.add(weight)
        
        var max_rev_weights = max(rev_weights)
        if max_rev_weights == 0'u32 and rev_end_reads.len != 0:
          rev_weights = @[]
          for rev_node in st.rev_edges[fwd_start_idx]:
            var weight = 0'u32
            for read_id in st.weights[(rev_node,fwd_start_idx)]:
              if read_id notin excluded_reads:
                weight += 1'u32
            rev_weights.add(weight)
          max_rev_weights = max(rev_weights)

        var max_fwd_weights = max(fwd_weights)
        if max_fwd_weights == 0'u32 and fwd_end_reads.len != 0:
          fwd_weights = @[]
          for fwd_node in st.fwd_edges[fwd_end_idx]:
            var weight = 0'u32
            for read_id in st.weights[(fwd_end_idx,fwd_node)]:
              if read_id notin excluded_reads:
                weight += 1'u32
            fwd_weights.add(weight)
          max_fwd_weights = max(fwd_weights)
        if max_rev_weights < max_fwd_weights:
          # Move in fwd direction
          if max_fwd_weights == 0:
            break
          next_node_id = st.fwd_edges[fwd_end_idx][maxIdx(fwd_weights)]
          for i in 0..<reads.len:
            let read_idx = read_idxs[i][1]
            if read_idx == -1:
              continue
            if read_idx + 1 != reads[i].corrected_path.len:
              if reads[i].corrected_path[read_idx + 1] == next_node_id:
                read_idxs[i][1] = read_idx + 1
                if (reads[i].corrected_path.len - read_idxs[i][1]) < psi:
                  read_idxs[i][1] = -1
                  fwd_end_reads.incl(reads[i].name)
              else:
                read_idxs[i][1] = -1
                if true:
                # if reads[i].name notin fwd_end_reads:
                  excluded_reads.incl(reads[i].name)
            else:
              read_idxs[i][1] = -1
          end_idx = int(st.node_indexes[('b',next_node_id)])
          # echo "end - ", end_idx
          if st.fwd_edges[next_node_id].len == 0:
            end_idx = -1
        else:
          # Move in rev direction
          if max_rev_weights == 0:
            break
          next_node_id = st.rev_edges[fwd_start_idx][maxIdx(rev_weights)]
          for i in 0..<reads.len:
            let read_idx = read_idxs[i][0]
            if read_idx == -1:
              continue
            if read_idx != 0:
              if reads[i].corrected_path[read_idx - 1] == next_node_id:
                read_idxs[i][0] = read_idx - 1
                if read_idxs[i][0] < psi:
                  read_idxs[i][0] = -1
                  rev_end_reads.incl(reads[i].name)
              else:
                read_idxs[i][0] = -1
                if true:
                # if reads[i].name notin rev_end_reads:
                  excluded_reads.incl(reads[i].name)
            else:
              read_idxs[i][0] = -1
          start_idx = int(st.node_indexes[('b',next_node_id)])
          # echo "start - ", start_idx
          if st.rev_edges[next_node_id].len == 0:
            start_idx = -1
      representative_path.add(next_node_id)

    let accounted_for_reads = all_reads - excluded_reads
    representative_paths.add(sorted(representative_path))
    var to_delete_reads : seq[int]
    for read_id in accounted_for_reads:
      to_delete_reads.add(read_id_to_idx[read_id])
    for idx in sorted(to_delete_reads,order = SortOrder.Descending):
      remaining_reads.delete(idx)
  return representative_paths

proc updateGraph( rep_po : var TrimmedPOGraph, new_path,old_path : seq[uint32] ) : TrimmedPOGraph =
  for i in 1..<new_path.len:
    let u = new_path[i-1]
    let v = new_path[i]
    if (u,v) in rep_po.weights:
      rep_po.weights[(u,v)] += 1'u32
    else:
      rep_po.weights[(u,v)] = 1'u32
      rep_po.edges[u].add(v)
  for i in 1..<old_path.len:
    let u = old_path[i-1]
    let v = old_path[i]
    rep_po.weights[(u,v)] -= 1'u32
    if rep_po.weights[(u,v)] == 0'u32:
      rep_po.weights.del((u,v))
      for i,j in rep_po.edges[u]:
        if j == v:
          rep_po.edges[u].delete(i,i)
          break
  return rep_po

proc getRepresentativePaths3( rep_po : var TrimmedPOGraph, psi : uint16 = 10) : seq[seq[uint32]] =
  var time1 = cpuTime()
  ########################################################################################
  ###------------------- Collect greedy walks from each source node -------------------###
  ########################################################################################
  for i in 0..<rep_po.nodes.len:
    rep_po.nodes[i].visited = false
  var source_nodes : seq[uint32]
  for read in rep_po.reads:
    source_nodes.add(read.corrected_path[0])
  var greedy_walks : seq[seq[uint32]]
  for i,fwd_node_idx in sorted(source_nodes):
    var path = @[fwd_node_idx]
    var bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx)]
    if rep_po.nodes[bck_node_idx].visited:
      continue
    rep_po.nodes[bck_node_idx].visited = true
    var fwd_node_idx2 = fwd_node_idx
    while rep_po.edges[fwd_node_idx2].len != 0:
      var weight_list : seq[uint32]
      for j in rep_po.edges[fwd_node_idx2]:
        weight_list.add(uint32(rep_po.weights[(fwd_node_idx2,uint32(j))]))
      fwd_node_idx2 = rep_po.edges[fwd_node_idx2][maxIdx(weight_list)]
      path.add(fwd_node_idx2)
      let bck_node_idx2 = rep_po.node_indexes[('b',fwd_node_idx2)]
      if rep_po.nodes[bck_node_idx2].visited:
        break
      rep_po.nodes[bck_node_idx2].visited = true
    greedy_walks.add(path)
  for i in 0..<rep_po.nodes.len:
    rep_po.nodes[i].visited = false
  echo "walk greedy paths - ", cpuTime() - time1
  time1 = cpuTime()
  ########################################################################################
  ###-------------------- Correct reads based on these greedy walks -------------------###
  ########################################################################################
  var differed_reads : seq[int]
  for i in 0..<rep_po.reads.len:
    var flag = true
    for j in 0..<greedy_walks.len:
      var differed = false
      let before_path = rep_po.reads[i].corrected_path
      (differed,rep_po.reads[i].corrected_path) = correctPaths(rep_po.reads[i].corrected_path,greedy_walks[j],addr rep_po.nodes,rep_po.node_indexes,psi=psi,correct_ends = false)
      rep_po = updateGraph(rep_po, rep_po.reads[i].corrected_path, before_path)
      flag = differed and flag
    if flag:
      differed_reads.add(i)
  
  echo "correct paths - ", cpuTime() - time1
  time1 = cpuTime()

  ########################################################################################
  ###--------------- Trim mini Illumina only paths based on greedy walks --------------###
  ########################################################################################

  for i in 0..<greedy_walks.len:
    var to_remove : seq[seq[uint32]]
    for minipath in rep_po.illumina_branches:
      ###### TODO: reduce the number of copying of minipath that needs to happen, redundant and slows us down.
      var minipath2 = minipath
      var (_,after_path) = correctPaths(minipath2, greedy_walks[i],addr rep_po.nodes, rep_po.node_indexes,psi = psi,correct_ends = false)
      if after_path != minipath:
        to_remove.add(minipath)
      for j in 1..<minipath.len - 1:
        rep_po.deleted_nodes.incl(minipath[j])
        rep_po.weights.del((minipath[j-1],minipath[j]))
        for l,k in rep_po.edges[minipath[j-1]]:
          if k == minipath[j]:
            rep_po.edges[minipath[j-1]].del(l)
            break
      rep_po.weights.del((minipath[^2],minipath[^1]))
      for l,k in rep_po.edges[minipath[^2]]:
        if k == minipath[^1]:
          rep_po.edges[minipath[^2]].del(l)
          break
    for s in to_remove:
      rep_po.illumina_branches.excl(s)
  echo "trim minipaths - ", cpuTime() - time1
  time1 = cpuTime()
  ########################################################################################
  ###---             Update log likelihoods based on illumina weights               ---###
  ########################################################################################
  for u in rep_po.edges.keys():
    var total : float64
    for v in rep_po.edges[u]:
      total += float64(rep_po.weights[(u,v)])
    for v in rep_po.edges[u]:
      rep_po.log_likelihoods[(u,v)] = ln(float(rep_po.weights[(u,v)]) / total)
  echo "update ll - ", cpuTime() - time1
  ########################################################################################
  ###--- Rewalk through the graph for each source node, constructing a p. queue of  ---###
  ###--- branches not walked.                                                       ---###
  ########################################################################################
  time1 = cpuTime()
  var alt_paths : HeapQueue[(int, float64, uint32, uint32)]
  # var alt_paths : HeapQueue[(int, uint32, uint32)]
  var alt_path_set : HashSet[(uint32,uint32)]
  for greedy_walk in greedy_walks:
    for i in 1..<greedy_walk.len:
      let u = greedy_walk[i-1]
      let bck_node_idx = rep_po.node_indexes[('b',u)]
      if rep_po.nodes[bck_node_idx].visited:
        break
      var first_node = 0
      if i - 1 > int(psi):
        first_node = int(i) - 1 - int(psi)
      rep_po.nodes[bck_node_idx].path = greedy_walk[first_node..i-1]
      rep_po.nodes[bck_node_idx].visited = true
      let v1 = greedy_walk[i]
      alt_path_set.incl((u,v1))
      for v2 in rep_po.edges[u]:
        if (u,v2) notin alt_path_set:
          alt_paths.push((-int(rep_po.weights[(u,v2)]),-rep_po.log_likelihoods[(u,v2)], v2, u))
          alt_path_set.incl((u,v2))
    # assert rep_po.edges[greedy_walk[^1]].len == 0
  echo "Construct P. queue - ", cpuTime() - time1
  
  ##########################################################################################################
  ###--- While the p. queue is not empty:                                                             ---###
  ###--- If the node has been deleted (meaning it was not part of a delta > psi in len), ignore it    ---###
  ###--- Similarly, if the transition from the preceeding node to the node have been deleted, ignore  ---###
  ###--- If the neither node nor the transition have been deleted:                                    ---###
  ###---     i - build a path from -1 of the node to the next node previously seen in the greedy walk ---###
  ###---         correct reads based on this short path                                               ---###
  ###---    ii - store new alternative paths to the priority queue                                    ---###
  ###---   iii - reconstruct the graph with the corrected paths, keeping track of the deleted nodes   ---###
  ###TODO Check the logic for this was ported correctly
  ##########################################################################################################
  time1 = cpuTime()
  var iterations = 0
  while alt_paths.len != 0:
    iterations += 1 
    var (_,_,fwd_node_idx,node1) = alt_paths.pop()
    if fwd_node_idx in rep_po.deleted_nodes or (node1,fwd_node_idx) notin rep_po.weights:
      continue
    var path = rep_po.nodes[rep_po.node_indexes[('b',node1)]].path
    path = path[^min(int(psi),path.len)..^1]  & @[fwd_node_idx]
    var bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx)]
    while rep_po.edges[fwd_node_idx].len != 0 and not rep_po.nodes[bck_node_idx].visited:
      var weight_list : seq[uint32]
      for j in rep_po.edges[fwd_node_idx]:
        weight_list.add(uint32(rep_po.weights[(fwd_node_idx,j)]))
      let last_fwd_node_idx = fwd_node_idx
      fwd_node_idx = rep_po.edges[fwd_node_idx][maxIdx(weight_list)]
      for i,next_node in rep_po.edges[last_fwd_node_idx]:
        if next_node == fwd_node_idx:
          alt_path_set.incl((next_node,
                             last_fwd_node_idx))
          continue
        if (next_node, last_fwd_node_idx) notin alt_path_set:
          alt_paths.push((-int(rep_po.weights[(last_fwd_node_idx,next_node)]), - rep_po.log_likelihoods[(last_fwd_node_idx,next_node)], next_node, last_fwd_node_idx))
          alt_path_set.incl((
                             next_node,
                             last_fwd_node_idx))
      rep_po.nodes[bck_node_idx].visited = true
      # po.nodes[fwd_node_idx].visited = true
      bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx)]
      path.add(fwd_node_idx)
    if rep_po.nodes[bck_node_idx].visited:
      for _ in 0..<int(psi):
        if rep_po.edges[fwd_node_idx].len == 0:
          break
        var weight_list : seq[uint32]
        for j in rep_po.edges[fwd_node_idx]:
          weight_list.add(uint32(rep_po.weights[(fwd_node_idx,j)]))
        fwd_node_idx = rep_po.edges[fwd_node_idx][maxIdx(weight_list)]
        bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx)]
        path.add(fwd_node_idx)
    # for idx in path:
    #   echo "path ", idx
    var reads : seq[Read]
    for i in differed_reads:
      var j : bool
      let before_path = rep_po.reads[i].corrected_path
      (j,rep_po.reads[i].corrected_path) = correctPaths(rep_po.reads[i].corrected_path,path,addr rep_po.nodes, rep_po.node_indexes, psi,false)
      rep_po = updateGraph(rep_po, rep_po.reads[i].corrected_path, before_path)
      reads.add(rep_po.reads[i])
    var to_remove : seq[seq[uint32]]
    for minipath in rep_po.illumina_branches:
      ###### TODO: reduce the number of copying of minipath that needs to happen, redundant and slows us down.
      var minipath2 = minipath
      var (_,after_path) = correctPaths(minipath2, path,addr rep_po.nodes, rep_po.node_indexes,psi = psi,correct_ends = false)
      if after_path != minipath:
        to_remove.add(minipath)
      for j in 1..<minipath.len - 1:
        rep_po.deleted_nodes.incl(minipath[j])
        rep_po.weights.del((minipath[j-1],minipath[j]))
        for l,k in rep_po.edges[minipath[j-1]]:
          if k == minipath[j]:
            rep_po.edges[minipath[j-1]].del(l)
            break
      rep_po.weights.del((minipath[^2],minipath[^1]))
      for l,k in rep_po.edges[minipath[^2]]:
        if k == minipath[^1]:
          rep_po.edges[minipath[^2]].del(l)
          break
    for s in to_remove:
      rep_po.illumina_branches.excl(s)


  echo "Navigate P. queue - ", cpuTime() - time1
  ########################################################################################
  ###------ Run approach similar to Stringtie approach for resolving splice graph -----###
  ###--- i.e. start at highest coverage node, walk heaviest path based on reads     ---###
  ###--- compatible with the walk up to that point. Remove reads compatible with the---###
  ###--- full walk, repeat.                                                         ---###
  ########################################################################################
  return walkHeaviestPaths( addr rep_po )

proc getRepresentativePaths2( rep_po : var TrimmedPOGraph, psi : uint16 = 10) : seq[seq[uint32]] =
  var source_nodes : seq[uint32]
  for read in rep_po.reads:
    source_nodes.add(read.corrected_path[0])
  var greedy_walks : seq[seq[uint32]]
  for i,fwd_node_idx in sorted(source_nodes):
    var path = @[fwd_node_idx]
    var bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx)]
    var fwd_node_idx2 = fwd_node_idx
    while rep_po.edges[fwd_node_idx2].len != 0:
      var weight_list : seq[uint32]
      for j in rep_po.edges[fwd_node_idx2]:
        weight_list.add(uint32(rep_po.weights[(fwd_node_idx2,uint32(j))]))
      fwd_node_idx2 = rep_po.edges[fwd_node_idx2][maxIdx(weight_list)]
      path.add(fwd_node_idx2)
    greedy_walks.add(path)
  ########################################################################################
  ###-------------------- Correct reads based on these greedy walks -------------------###
  ########################################################################################
  var differed_reads : seq[int]
  for i in 0..<rep_po.reads.len:
    var flag = true
    for j in 0..<greedy_walks.len:
      var differed = false
      # for k in 1..<greedy_walks[j].len:
      #   if greedy_walks[j][k-1] == 15813'u32 and greedy_walks[j][k] == 15834'u32:
      #     echo "What the actual living fuck"
      # for k in 1..<rep_po.reads[i].corrected_path.len:
      #   if rep_po.reads[i].corrected_path[k-1] == 15813'u32 and rep_po.reads[i].corrected_path[k] == 15834'u32:
      #     echo "What the actual living fuck2"
      let before_path = rep_po.reads[i].corrected_path
      (differed,rep_po.reads[i].corrected_path) = correctPaths(rep_po.reads[i].corrected_path,greedy_walks[j],addr rep_po.nodes,rep_po.node_indexes,psi=psi,correct_ends = false)
      rep_po = updateGraph(rep_po, rep_po.reads[i].corrected_path, before_path)
      flag = differed and flag
      # for k in 1..<rep_po.reads[i].corrected_path.len:
      #   if rep_po.reads[i].corrected_path[k-1] == 15813'u32 and rep_po.reads[i].corrected_path[k] == 15834'u32:
      #     echo "What the actual living fuck3"
    if flag:
      differed_reads.add(i)
  ########################################################################################
  ###----- For reads that don't match a greedy walk: use as template and fix! ---------###
  ########################################################################################
  var count = 0
  while differed_reads.len != 0:
    count += 1
    # echo "here!!!"
    # echo differed_reads
    var still_differed : seq[bool]
    for _ in differed_reads:
      still_differed.add(false)
    let template_read = rep_po.reads[differed_reads[0]].corrected_path
    # echo rep_po.reads[differed_reads[0]].corrected_path
    var read_set : Table[uint32,uint32]
    for j,k in template_read:
      read_set[uint32(k)] = uint32(j)
    var node = template_read[0]
    var new_path : seq[uint32]
    var branch_count = 0
    while rep_po.edges[node].len != 0:
      # echo node
      if branch_count == 0:
        new_path.add(node)
        # echo "n - ", node, " ", rep_po.edges[node]
      if rep_po.edges[node].len == 1:
        node = rep_po.edges[node][0]
        # if node in read_set:
        #   if read_set[node] == uint32(template_read.len - 1):
        #     echo "got here"
          #   break
      else:
        # var end_flag = false
        var correct_flag = false
        var minipath : seq[uint32]
        minipath.add(node)
        if node notin read_set:
          break
        # if node in read_set:
        #   if read_set[node] == uint32(template_read.len - 1):
        #     echo "got here"
        let last_node_j = read_set[node]
        var weight_list : seq[uint32]
        for j in rep_po.edges[node]:
          weight_list.add(uint32(rep_po.weights[(node,uint32(j))]))
        # echo "edges - ", rep_po.edges[node]
        # echo "weights -", weight_list
        if branch_count < weight_list.len:
          node = rep_po.edges[node][orderedIdx(weight_list,branch_count)]
        else:
          # echo "here2"
          # break
          node = rep_po.edges[node][maxIdx(weight_list)]
          correct_flag = true
          # echo "got here"
        minipath.add(node)
        if node in read_set:
          if read_set[node] - last_node_j < uint32(psi):
            correct_flag = true
          # if read_set[node] == uint32(template_read.len - 1):
          #   echo "got here"
        if not correct_flag:
          for k in 1'u16..<psi:
            # echo "k - ", k
            var weight_list : seq[uint32]
            if rep_po.edges[node].len == 0:
              break
            for j in rep_po.edges[node]:
              weight_list.add(uint32(rep_po.weights[(node,uint32(j))]))
            # echo rep_po.edges[node]
            # echo weight_list
            node = rep_po.edges[node][maxIdx(weight_list)]
            minipath.add(node)
            if node in read_set:
              if read_set[node] - last_node_j < uint32(psi):
                correct_flag = true
              # if read_set[node] == uint32(template_read.len - 1):
              #   echo "got here"
              break
        if correct_flag:
          # echo "correcting"
          for i in 1..<minipath.len-1:
            # echo "c - ", minipath[i]
            new_path.add(minipath[i])
          branch_count = 0
        else:
          node = template_read[read_set[minipath[0]]]
          branch_count += 1
    # var outfile : File
    # discard open(outfile,&"test_po_trim_{count}.html",fmWrite)
    # htmlOutput(rep_po, outfile,rep_po.reads[differed_reads[0]].corrected_path,new_path)
    var new_differed_reads : seq[int]
    # echo new_path
    for i in differed_reads:
      var differed = false
      let before_path = rep_po.reads[i].corrected_path
      (differed,rep_po.reads[i].corrected_path) = correctPaths(rep_po.reads[i].corrected_path, new_path, addr rep_po.nodes, rep_po.node_indexes, psi=psi,correct_ends = false)
      rep_po = updateGraph(rep_po, rep_po.reads[i].corrected_path, before_path)
      if differed:
        new_differed_reads.add(i)
    assert new_differed_reads != differed_reads
    # for i,j in still_differed:
    #   if j:
    #     new_differed_reads.add(differed_reads[i])
    differed_reads = new_differed_reads
  # echo "here"
  return walkHeaviestPaths( addr rep_po )

proc getRepresentativePaths(po : var POGraph, psi : uint16 = 10) : seq[seq[uint32]] =
  #[
  From a given Partial Order graph, fetches a set of representative paths for that graph, ie paths that capture the major variations between isoforms
  ]#
  var reads : seq[Read]
  for i in 0..<po.reads.len:
    po.reads[i].corrected_path = po.reads[i].path
    reads.add(po.reads[i])
  echo reads.len
  var rep_po = constructNewGraph(addr po,reads)
  
  ########################################################################################
  ###------------------- Collect greedy walks from each source node -------------------###
  ########################################################################################
  ### TODO? Wrap this in a function? More readable but might slow down excecution since we'd need to return rep_po since nim is making deep copies when passing to functions... Or just use pointers you dummy!
  ### Come back to this depending on what speed / memory usage looks like when we're done
  
  var greedy_walks : seq[seq[uint32]]
  for i,fwd_node_idx in sorted(rep_po.source_nodes):
    var path = @[fwd_node_idx]
    var bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx)]
    rep_po.nodes[bck_node_idx].visited = true
    po.nodes[fwd_node_idx].visited = true
    var fwd_node_idx2 = fwd_node_idx
    while rep_po.edges[fwd_node_idx2].len != 0:
      var weight_list : seq[uint32]
      for j in rep_po.edges[fwd_node_idx2]:
        weight_list.add(uint32(rep_po.weights[(fwd_node_idx2,uint32(j))]))
      fwd_node_idx2 = rep_po.edges[fwd_node_idx2][maxIdx(weight_list)]
      bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx2)]
      rep_po.nodes[bck_node_idx].visited = true
      po.nodes[fwd_node_idx2].visited = true
      path.add(fwd_node_idx2)
    greedy_walks.add(path)
  ########################################################################################
  ###-------------------- Correct reads based on these greedy walks -------------------###
  ########################################################################################
  var differed_reads : seq[int]
  for i in 0..<po.reads.len:
    var flag = true
    for j in 0..<greedy_walks.len:
      var differed = false
      (differed,po.reads[i].corrected_path) = correctPaths(po.reads[i].corrected_path,greedy_walks[j],addr rep_po.nodes,rep_po.node_indexes,psi=psi,correct_ends = true)
      # echo po.reads[i].corrected_path.len
      flag = differed and flag
    # if flag and (po.reads[i].corrected_path.len != 0):
    #   differed_reads.add(i)
    if flag:
      differed_reads.add(i)
  reads = @[]
  for i in differed_reads:
    reads.add(po.reads[i])
  ########################################################################################
  ###--- Reconstruct graph based on corrected walks, keeping track of deleted nodes ---###
  ########################################################################################
  rep_po = constructNewGraph(addr po,reads,rep_po.deleted_nodes)
  ########################################################################################
  ###--- Rewalk through the graph for each source node, constructing a p. queue of  ---###
  ###--- branches not walked. (NOTE: could do this above, but it would create a bu- ---###
  ###--- nch of branches that end up being deleted after correction, so I think its ---###
  ###--- faster this way. Need to test that theory during optimization;             ---###
  ####################################### TODO ###########################################
  ########################################################################################
  var alt_paths : HeapQueue[(int, float64, uint32, uint32)]
  var alt_path_set : HashSet[(uint32,uint32)]
  for i,fwd_node_idx in rep_po.source_nodes:
    let bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx)]
    rep_po.nodes[bck_node_idx].path = @[fwd_node_idx]
    po.nodes[fwd_node_idx].path = @[fwd_node_idx]
    rep_po.nodes[bck_node_idx].visited = true
    po.nodes[fwd_node_idx].visited = true
    var fwd_node_idx2 = fwd_node_idx
    var last_fwd_node_idx = fwd_node_idx2
    while rep_po.edges[fwd_node_idx2].len != 0:
      var weight_list : seq[uint32]
      for j in rep_po.edges[fwd_node_idx2]:
        weight_list.add(uint32(rep_po.weights[(fwd_node_idx2,j)]))
      last_fwd_node_idx = fwd_node_idx2
      fwd_node_idx2 = rep_po.edges[fwd_node_idx2][maxIdx(weight_list)]
      for i,next_node in rep_po.edges[last_fwd_node_idx]:
        rep_po.nodes[rep_po.node_indexes[('b',next_node)]].path = rep_po.nodes[rep_po.node_indexes[('b',last_fwd_node_idx)]].path & @[next_node]
        po.nodes[next_node].path = rep_po.nodes[rep_po.node_indexes[('b',last_fwd_node_idx)]].path & @[next_node]
        if next_node == fwd_node_idx2:
          alt_path_set.incl((next_node,
                             last_fwd_node_idx))
          continue
        if (next_node, last_fwd_node_idx) notin alt_path_set:
          alt_paths.push((-int(rep_po.weights[(last_fwd_node_idx,next_node)]), - rep_po.log_likelihoods[(last_fwd_node_idx,next_node)], next_node, last_fwd_node_idx))
          alt_path_set.incl((next_node,
                             last_fwd_node_idx))
  
  ##########################################################################################################
  ###--- While the p. queue is not empty:                                                             ---###
  ###--- If the node has been deleted (meaning it was not part of a delta > psi in len), ignore it    ---###
  ###--- Similarly, if the transition from the preceeding node to the node have been deleted, ignore  ---###
  ###--- If the neither node nor the transition have been deleted:                                    ---###
  ###---     i - build a path from -1 of the node to the next node previously seen in the greedy walk ---###
  ###---         correct reads based on this short path                                               ---###
  ###---    ii - store new alternative paths to the priority queue                                    ---###
  ###---   iii - reconstruct the graph with the corrected paths, keeping track of the deleted nodes   ---###
  ##########################################################################################################
  var iterations = 0
  while alt_paths.len != 0:
    iterations += 1 
    var (_,_,fwd_node_idx,node1) = alt_paths.pop()
    # echo test1," ",fwd_node_idx," ",node1
    if fwd_node_idx in rep_po.deleted_nodes or (node1,fwd_node_idx) notin rep_po.weights:
      continue
    # echo "iter ", iterations
    var path = rep_po.nodes[rep_po.node_indexes[('b',node1)]].path
    path = path[^min(int(psi),path.len)..^1]  & @[fwd_node_idx]
    var bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx)]
    while rep_po.edges[fwd_node_idx].len != 0 and not rep_po.nodes[bck_node_idx].visited:
      var weight_list : seq[uint32]
      for j in rep_po.edges[fwd_node_idx]:
        weight_list.add(uint32(rep_po.weights[(fwd_node_idx,j)]))
      let last_fwd_node_idx = fwd_node_idx
      fwd_node_idx = rep_po.edges[fwd_node_idx][maxIdx(weight_list)]
      for i,next_node in rep_po.edges[last_fwd_node_idx]:
        if next_node == fwd_node_idx:
          alt_path_set.incl((next_node,
                             last_fwd_node_idx))
          continue
        if (next_node, last_fwd_node_idx) notin alt_path_set:
          alt_paths.push((-int(rep_po.weights[(last_fwd_node_idx,next_node)]), - rep_po.log_likelihoods[(last_fwd_node_idx,next_node)], next_node, last_fwd_node_idx))
          alt_path_set.incl((
                             next_node,
                             last_fwd_node_idx))
      rep_po.nodes[bck_node_idx].visited = true
      po.nodes[fwd_node_idx].visited = true
      bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx)]
      path.add(fwd_node_idx)
    if rep_po.nodes[bck_node_idx].visited:
      for _ in 0..<int(psi):
        if rep_po.edges[fwd_node_idx].len == 0:
          break
        var weight_list : seq[uint32]
        for j in rep_po.edges[fwd_node_idx]:
          weight_list.add(uint32(rep_po.weights[(fwd_node_idx,j)]))
        fwd_node_idx = rep_po.edges[fwd_node_idx][maxIdx(weight_list)]
        bck_node_idx = rep_po.node_indexes[('b',fwd_node_idx)]
        path.add(fwd_node_idx)
    # for idx in path:
    #   echo "path ", idx
    reads = @[]
    for i in differed_reads:
      var j : bool
      (j,po.reads[i].corrected_path) = correctPaths(po.reads[i].corrected_path,path,addr po.nodes, rep_po.node_indexes, psi,false)
      # echo po.reads[i].name, " - " ,po.reads[i].corrected_path.len
      # if po.reads[i].corrected_path.len > int(psi):
      #   reads.add(po.reads[i])
      reads.add(po.reads[i])
    rep_po = constructNewGraph(addr po, reads, rep_po.deleted_nodes)
  
  ########################################################################################
  ###------ Run approach similar to Stringtie approach for resolving splice graph -----###
  ###--- i.e. start at highest coverage node, walk heaviest path based on reads     ---###
  ###--- compatible with the walk up to that point. Remove reads compatible with the---###
  ###--- full walk, repeat.                                                         ---###
  ########################################################################################
  return walkHeaviestPaths(addr po)

proc scorePaths( x : seq[uint32],y : seq[uint32]) : int =
  var i,j = 0
  while true:
    if x[i] == y[j]:
      result += 1
      i += 1
      j += 1
    elif x[i] < y[j]:
      result -= 1
      i += 1
    elif x[i] > y[j]:
      result -= 1
      j += 1
    if i == len(x) or j == len(y):
      break

proc assignReadsToPaths( reads : seq[Read], paths : seq[seq[uint32]]) : seq[int] = 
  var closest_paths : seq[int]
  for read in reads:
    var path_scores : seq[int]
    for path in paths:
      path_scores.add(scorePaths(read.path,path))
    closest_paths.add(maxIdx(path_scores))
  return closest_paths

proc constructGraphFromPaths( po : ptr POGraph, paths : seq[seq[uint32]]) : TrimmedPOGraph = 
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  var u_set,v_set : HashSet[uint32]
  for path in paths:
    for i in 0..<path.len - 1:
      if (path[i],path[i+1]) in weights:
        weights[(path[i],path[i+1])] += 1'u32
      else:
        weights[(path[i],path[i+1])] = 1'u32
        if path[i] in edges:
          edges[path[i]].add(path[i+1])
        else:
          edges[path[i]] = @[path[i+1]]
        u_set.incl(path[i])
        v_set.incl(path[i+1])
    # for i in 0..<read.corrected_path.len:
    #   echo read.corrected_path[i]
  var source_nodes1,source_nodes2 : HashSet[uint32]
  for i,node in po[].nodes:
    if node.start_node_flag:
      source_nodes1.incl(uint32(i))
  let sink_nodes = v_set - u_set
  for sink_node in sink_nodes:
    edges[sink_node] = @[]
  var nodes : seq[Node]
  var node_indexes = sorted(toSeq(u_set + v_set))
  var node_indexes2 : Table[(char,uint32),uint32]
  var log_likelihoods : Table[(uint32,uint32),float64]
  for i,j in node_indexes:
    nodes.add(po[].nodes[j])
    if j in source_nodes1:
      source_nodes2.incl(j)
    node_indexes2[('f',uint32(i))] = j
    node_indexes2[('b',j)] = uint32(i)
    if j in edges:
      var total_read_count = 0'u32
      for k in edges[j]:
        total_read_count += weights[(j,k)]
      let tmp = float(total_read_count)
      for k in edges[j]:
        log_likelihoods[(j,k)] = - ln(float(weights[(j,k)]) / tmp)
  var end_nodes : seq[uint32]
  for i in 0..<nodes.len:
    nodes[i].indegree = 0'u16
    if nodes[i].end_node_flag:
      end_nodes.add(node_indexes2[('f',uint32(i))])
  for j in edges.keys:
    for k in edges[j]:
      nodes[node_indexes2[('b',k)]].indegree += 1'u16
  # for j in sorted(toSeq(edges.keys)):
  #   for k in sorted(toSeq(edges[j])):
  #     echo j, " ", k, " ", weights[(j,k)]
  return TrimmedPOGraph(nodes : nodes,
                        edges : edges,
                        weights : weights,
                        reads : po[].reads,
                        og_nodes : po[].og_nodes,
                        log_likelihoods : log_likelihoods,
                        source_nodes : toSeq(source_nodes2),
                        end_nodes : end_nodes,
                        node_indexes : node_indexes2)

proc reducePath( path : seq[uint32],nodes : seq[Node], node_indexes : Table[(char,uint32), uint32]) : (seq[uint32],string) = 
  var new_path : seq[uint32]
  var i = 0
  var seq_list : seq[string]
  while i < len(path):
    let node = nodes[node_indexes[('b',path[i])]]
    new_path.add(node_indexes[('b',path[i])])
    seq_list.add(node.nts)
    i += node.nts.len
  # echo new_path
  return (new_path, seq_list.join(""))

proc trimAndCollapsePOGraph2( po : var POGraph, psi:uint16 = 10) : seq[seq[uint32]] = 
  let paths = getRepresentativePaths(po, psi = psi)
  
  var selected_paths : seq[seq[uint32]]
  for read in po.reads:
    selected_paths.add(read.corrected_path)
  
  var rep_po = constructGraphFromPaths(addr po, selected_paths)
  return paths

proc trimAndCollapsePOGraph( po : var POGraph, psi:uint16 = 10) : (TrimmedPOGraph, seq[seq[uint32]]) = 
  let paths = getRepresentativePaths(po, psi = psi)
  # getRepresentativePaths(po,psi=psi)
  # let closest = assignReadsToPaths(po.reads,paths)
  # var selected_paths : seq[seq[uint32]]
  # for c in closest:
  #   selected_paths.add(paths[c])
  #   # for idx in paths[c]:
  #   #   echo idx
  var selected_paths : seq[seq[uint32]]
  for read in po.reads:
    selected_paths.add(read.corrected_path)
  
  var rep_po = constructGraphFromPaths(addr po, selected_paths)
  # for idx in sorted(toSeq(rep_po.node_indexes.keys)):
  #   echo "idx, ", idx, " ", rep_po.node_indexes[idx]
  # for read in rep_po.reads:
  #   rep_po.nodes[rep_po.node_indexes[('b',read.corrected_path[0])]].start_node_flag = true
  #   rep_po.nodes[rep_po.node_indexes[('b',read.corrected_path[^1])]].end_node_flag = true
    
  # rep_po = collapseLinearStretches(rep_po)
  # for idx in sorted(toSeq(rep_po.node_indexes.keys)):
  #   echo "idx, ", idx, " ", rep_po.node_indexes[idx]
  return (rep_po,selected_paths)

proc insertBasesToRead( path : var seq[uint32],pos : uint32, inserts : seq[uint32]) : (seq[uint32],uint32,bool) = 
  var next : uint32
  var changed = false
  for i,node in path:
    if node == pos:
      next = path[i+1]
      # echo next
      # for insert in sorted(inserts,Descending):
      path.insert(inserts,i)
      changed = true
      break
  return (path,next,changed)

proc mismatchBasesToGraph( po : ptr TrimmedPOGraph, pos : uint32, base : string) : int = 
  var new_node_index : uint32
  var node_index = po[].nodes.len
  var node_number = po[].og_nodes
  let new_node = Node(nts : base,
                      align_ring_partner : int32(-1),
                      visited : false,
                      indegree : 1'u16,
                      log_likelihood : 0.0,
                      start_node_flag : false,
                      end_node_flag : false)

  po[].nodes.add(new_node)
  new_node_index = node_number
  # if new_node_index == 5614:
  #   echo "here1"
  #   echo base
  po[].node_indexes[('f',uint32(node_index))] = uint32(node_number)
  po[].node_indexes[('b',uint32(node_number))] = uint32(node_index)
  po[].edges[new_node_index] = @[]
  node_index += 1
  node_number += 1

  if po[].nodes[po[].node_indexes[('b', pos)]].align_ring_partner == -1'i32:
    po[].nodes[po[].node_indexes[('b', pos)]].align_ring_partner = int32(new_node_index)
    po[].nodes[po[].node_indexes[('b',new_node_index)]].align_ring_partner = int32(pos)
  else:
    var next_node = po[].node_indexes[('b',uint32(po[].nodes[po.node_indexes[('b',pos)]].align_ring_partner))]
    while po[].nodes[next_node].align_ring_partner != int32(pos):
      next_node = po[].node_indexes[('b',uint32(po[].nodes[next_node].align_ring_partner))]
    po[].nodes[next_node].align_ring_partner = int32(new_node_index)
    po[].nodes[po[].node_indexes[('b',new_node_index)]].align_ring_partner = int32(pos)

  po[].og_nodes += 1'u32

  return int(new_node_index)

proc insertBasesToGraph( po : ptr TrimmedPOGraph, pos : uint32, insert : string) : seq[uint32] = 
  var new_node_indexes : seq[uint32]
  var node_index = po[].nodes.len
  var node_number = po[].og_nodes
  for base in insert:
    let new_node = Node(nts : $base,
                        align_ring_partner : int32(-1),
                        visited : false,
                        indegree : 1'u16,
                        log_likelihood : 0.0,
                        start_node_flag : false,
                        end_node_flag : false)

    po[].nodes.add(new_node)
    new_node_indexes.add(node_number)
    # if node_number == 5614:
    #   echo "here2"
    #   echo insert
    po[].node_indexes[('f',uint32(node_index))] = uint32(node_number)
    po[].node_indexes[('b',uint32(node_number))] = uint32(node_index)
    po[].edges[node_number] = @[]
    node_index += 1
    node_number += 1
  po[].og_nodes += uint32(insert.len)

  return new_node_indexes

proc topologicalSort( po : ptr TrimmedPOGraph) =
  var nodes_with_no_incoming_edges : seq[uint32]
  var topo_ordering : seq[uint32]
  var indegrees : seq[uint16]
  for _ in po[].nodes:
    indegrees.add(0'u16)
  for i,_ in po[].nodes:
    # echo po.node_indexes[('f',uint32(i))]
    for fwd_neighbor in po[].edges[po.node_indexes[('f',uint32(i))]]:
      let bck_neighbor = po[].node_indexes[('b',fwd_neighbor)]
      indegrees[bck_neighbor] += 1'u16
  # echo indegrees
  # echo po.edges
  for i in 0..<po[].nodes.len:
    if indegrees[i] == 0'u16:
      nodes_with_no_incoming_edges.add(uint32(i))
  while nodes_with_no_incoming_edges.len > 0:
    let node = nodes_with_no_incoming_edges.pop()
    topo_ordering.add(node)
    # echo "node\t", node
    for fwd_neighbor in po[].edges[po[].node_indexes[('f',node)]]:
      # echo "fwd\t", fwd_neighbor
      let bck_neighbor = po[].node_indexes[('b',fwd_neighbor)]
      # echo "bck\t", bck_neighbor
      indegrees[bck_neighbor] -= 1'u16
      if indegrees[bck_neighbor] == 0'u16:
        nodes_with_no_incoming_edges.add(bck_neighbor)
  # for i,indegree in indegrees:
  #   if indegree != 0'u16:
  #     echo &"Improper Node: {i}"
  #     echo indegree
  assert topo_ordering.len == po[].nodes.len

  var new_node_indexes : Table[(char,uint32),uint32]
  for i,j in topo_ordering:
    let bck_index = j
    let old_fwd_index = po[].node_indexes[('f',bck_index)]
    let new_fwd_index = uint32(i)
    # if ((3000'u32 < old_fwd_index) and (old_fwd_index < 3030'u32)) or (old_fwd_index == 6456'u32):
    #   echo old_fwd_index, "\t", new_fwd_index
    new_node_indexes[('b',new_fwd_index)] = bck_index
    new_node_indexes[('f',bck_index)] = new_fwd_index
  var new_edges : Table[uint32,seq[uint32]]
  var new_weights : Table[(uint32,uint32),uint32]
  for old_u in po[].edges.keys:
    let new_u = new_node_indexes[('f',po.node_indexes[('b',old_u)])]
    new_edges[new_u] = @[]
    for old_v in po[].edges[old_u]:
      let new_v = new_node_indexes[('f',po[].node_indexes[('b',old_v)])]
      new_edges[new_u].add(new_v)
      # echo old_u, " - ", old_v
      new_weights[(new_u,new_v)] = po[].weights[(old_u,old_v)]

  for i in 0..<po.reads.len:
    var new_path : seq[uint32]
    var print_flag = false
    for old_node in po[].reads[i].corrected_path:
      let new_node = new_node_indexes[('f',po[].node_indexes[('b',old_node)])]
      if new_node == 15813'u32:
        print_flag = true
      new_path.add(new_node)
    # if print_flag:
    #   echo "old - ", po.reads[i].corrected_path
    #   echo "new - ", new_path
    po[].reads[i].corrected_path = new_path
  var new_illumina_branches : HashSet[seq[uint32]]
  for s in po[].illumina_branches:
    var new_path : seq[uint32]
    for u in s:
      new_path.add(new_node_indexes[('f',po[].node_indexes[('b',u)])])
    new_illumina_branches.incl(new_path)
  po[].node_indexes = new_node_indexes
  po[].edges = new_edges
  po[].weights = new_weights
  po[].illumina_branches = new_illumina_branches

proc illuminaPolishPOGraph( po : ptr TrimmedPOGraph, bam:Bam, illumina_weight:uint32 = 10,debug=true)=
  var read_id_to_idx : Table[string,uint32]
  var read_paths : seq[seq[uint32]]
  for i in 0..<po.nodes.len:
    po[].nodes[i].align_ring_partner = -1'i32
  for i,read in po.reads:
    read_id_to_idx[read.name] = uint32(i)
    read_paths.add(read.corrected_path.deepCopy())
  
  var inserts : Table[(uint32,string),seq[uint32]]
  var deletes : HashSet[seq[uint32]]
  var illumina_branches : HashSet[seq[uint32]]
  # var total_first_mapping_time = 0.0
  # var total_secondary_mapping_time = 0.0
  for read in po.reads:
    # echo read.name
    for record in bam.query(read.name):
      # let time1 = cpuTime()
      if record.flag.unmapped:
        continue
      # if record.chrom == "":
      #   continue
      # echo ""
      # echo record.cigar
      # if not (record.chrom in read_id_to_idx):
      #   continue
      # if record.mapping_quality == 0:
      #   continue
      # let path = read_paths[read_id_to_idx[record.chrom]]
      let path = po[].reads[read_id_to_idx[record.chrom]].corrected_path
      let cigar = record.cigar
      var ref_index = record.start
      var que_index = 0
      var traveled_nodes : seq[uint32]
      var new_edges : seq[(uint32,uint32)]
      var last_delete = (0,0)
      var illumina_walks : seq[(uint16,uint16)]

      var new_edge_flag = false
      var report_flag = false
      var delete_flag = false
      var ambiguous_flag = false
      var new_block_flag = false
      var new_block_start = -1
      var new_block_offset = 0
      # var s : string
      # record.sequence(s)
      # echo  s
      # echo ""
      for op in cigar:
        case int(op.op):
          of 1: #Insert
            if not new_block_flag:
              new_block_start = que_index - 1 - new_block_offset
              new_block_flag = true
            var sequence : seq[string]
            for i in 0..<op.len:
              if record.base_at(que_index) == 'N':
                ambiguous_flag = true
                break
              sequence.add($record.base_at(que_index))
              que_index += 1
            if ambiguous_flag:
              break
            var insert_walk : seq[uint32]
            if not ((traveled_nodes[^1],sequence.join()) in inserts):
              insert_walk = insertBasesToGraph(po,traveled_nodes[^1],sequence.join())
              inserts[(traveled_nodes[^1],sequence.join())] = insert_walk
              new_edges.add((traveled_nodes[^1],insert_walk[0]))
              for i in 1..<insert_walk.len:
                new_edges.add((insert_walk[i-1], insert_walk[i]))
              new_edge_flag = true
            else:
              insert_walk = inserts[(traveled_nodes[^1],sequence.join())]
            if delete_flag:
              let delete_path = @[traveled_nodes[^1],insert_walk[0]]
              if delete_path notin deletes:
                deletes.incl(delete_path)
                new_edges.add((traveled_nodes[^1],insert_walk[0]))
              delete_flag = false
            for node in insert_walk:
              traveled_nodes.add(node)
              # illumina_walk.add(true)
          
          of 2: #Delete
            if not new_block_flag:
              new_block_flag = true
              new_block_start = que_index - 1 - new_block_offset
            # let start_path = traveled_nodes[^1]
            ref_index += op.len
            delete_flag = true
            # let end_path = path[uint32(ref_index + 1)]
            
            # if new_edge_flag:
            #   new_edges.add((traveled_nodes[^1],path[ref_index]))
            #   new_edge_flag = false

            # let delete_path = @[start_path,end_path]
            # if delete_path notin deletes:
            #   new_edge_flag = true
            #   deletes.incl(delete_path)
          
          # of 3: #Intron # Shouldn't be possible for us

          of 4: #Soft-clip
            # echo "Soft-clip"
            que_index += op.len
            new_block_offset += op.len
          
          # of 5: #Hard-clip # Shouldn't be possible for us

          # of 6: #Pad # Shouldn't be possible for us
          
          of 7: #Equal
            # if debug:
            #   echo "Equal"
            if new_block_flag:
              illumina_walks.add((uint16(new_block_start),uint16(que_index - new_block_offset)))
              new_block_start = -1
              new_block_flag = false
            if delete_flag:
              let delete_path = @[traveled_nodes[^1],path[ref_index]]
              if delete_path notin deletes:
                deletes.incl(delete_path)
                new_edges.add((traveled_nodes[^1],path[ref_index]))
            if new_edge_flag:
              new_edges.add((traveled_nodes[^1], path[ref_index]))
              new_edge_flag = false
            for i in 0..<op.len:
              if debug:
                # echo "Match!"
                # echo po.nodes[po[].node_indexes[('b',path[ref_index])]].nts
                # echo $record.base_at(que_index)
                # echo ref_index
                assert po[].nodes[po[].node_indexes[('b',path[ref_index])]].nts == $record.base_at(que_index)
              traveled_nodes.add(path[ref_index])
              ref_index += 1
              que_index += 1
          of 8: #Diff
            # if debug:
            #   echo "Diff"
            if not new_block_flag:
              new_block_start = que_index - 1 - new_block_offset
              new_block_flag = true
            for i in 0..<op.len:
              if record.base_at(que_index) == 'N':
                ambiguous_flag = true
                break
              # if test.nodes[path[ref_index]].nts != $record.base_at(que_index):
              if debug:
                # echo "Mismatch!"
                # echo po.nodes[po.node_indexes[('b',path[ref_index])]].nts
                # echo $record.base_at(que_index)
                # if path[ref_index] == 3017:
                #   echo "here7"
                #   echo po.nodes[po.node_indexes[('b',path[ref_index])]].nts
                #   echo $record.base_at(que_index)
                assert po[].nodes[po[].node_indexes[('b',path[ref_index])]].nts != $record.base_at(que_index)
                # echo "Diff - Mismatch! - ", po.nodes[po.node_indexes[('b',path[ref_index])]].nts, $record.base_at(que_index)
                # test.nodes[path[ref_index]].nts = $record.base_at(que_index)
                # po.nodes[po.node_indexes[('b',path[ref_index])]].nts = $record.base_at(que_index)
              var next_node = po[].node_indexes[('b',path[ref_index])]
              var correct_node = -1
              for i in 0..<3:
                if po[].nodes[next_node].align_ring_partner == -1'i32:
                  break
                if po[].nodes[next_node].align_ring_partner == int32(path[ref_index]):
                  break
                if po[].nodes[po[].node_indexes[('b',uint32(po[].nodes[next_node].align_ring_partner))]].nts == $record.base_at(que_index):
                  correct_node = po[].nodes[next_node].align_ring_partner
                  break
                else:
                  next_node = po[].node_indexes[('b',uint32(po[].nodes[next_node].align_ring_partner))]
              if correct_node == -1:
                correct_node = mismatchBasesToGraph(po,path[ref_index],$record.base_at(que_index))
                if traveled_nodes.len > 0:
                  new_edges.add((traveled_nodes[^1],uint32(correct_node)))
                new_edge_flag = true
              elif new_edge_flag:
                new_edges.add((traveled_nodes[^1],uint32(correct_node)))
              if delete_flag:
                let delete_path = @[traveled_nodes[^1],uint32(correct_node)]
                if delete_path notin deletes:
                  deletes.incl(delete_path)
                  new_edges.add((traveled_nodes[^1],uint32(correct_node)))
              traveled_nodes.add(uint32(correct_node))
              ref_index += 1
              que_index += 1
            if ambiguous_flag:
              break
          # of 9: #Back # Shouldn't be possible for us
          else:
            echo "Something went wrong"
      
      if ambiguous_flag:
        continue
      var block_count = 0
      # echo record.cigar
      # echo traveled_nodes
      for (s,e) in illumina_walks:
        illumina_branches.incl(traveled_nodes[s..e])
      for (u,v) in new_edges:
        if u in po[].edges:
          po[].edges[u].add(v)
      for i in 1..<traveled_nodes.len:
        let u = traveled_nodes[i-1]
        let v = traveled_nodes[i]
        if (u,v) in po[].weights:
          po[].weights[(u,v)] += illumina_weight
        else:
          po[].weights[(u,v)] = illumina_weight
      # po = topologicalSort(po)
      # if debug:
      #   assert que_index == s.len
      #   var count = 0
      #   let ops = [8,1,7,1]
      #   let lens = [1,3,2,2]
      #   var flag = false
      #   for op in record.cigar:
      #     if int(op.op) == ops[count] and op.len == lens[count]:
      #       if count < 3:
      #         count += 1
      #       else:
      #         flag = true
      #         break
      #   if flag:
      #     echo "here5"
      #     echo record.cigar
      #     # echo traveled_nodes
      #     for i in 1..<traveled_nodes.len:
      #       let node1 = traveled_nodes[i-1]
      #       let node2 = traveled_nodes[i]
      #       echo node1, "\t", po.nodes[po.node_indexes[('b',node1)]].nts
      #       echo po.edges[node1]
      #       for v in po.edges[node1]:
      #         echo node1, "\t", v, "\t",po.nodes[po.node_indexes[('b',v)]].nts, "\t", po.weights[(node1,v)] 
      #     echo traveled_nodes[^1], "\t", po.nodes[po.node_indexes[('b',traveled_nodes[^1])]].nts
      #   if report_flag:
      #     echo "here4"
      #     echo record.cigar
      #     # echo traveled_nodes
      #     for i in 1..<traveled_nodes.len:
      #       let node1 = traveled_nodes[i-1]
      #       let node2 = traveled_nodes[i]
      #       echo node1, "\t", po.nodes[po.node_indexes[('b',node1)]].nts
      #       echo po.edges[node1]
      #       for v in po.edges[node1]:
      #         echo node1, "\t", v, "\t", po.nodes[po.node_indexes[('b',v)]].nts, "\t", po.weights[(node1,v)] 
      #     echo traveled_nodes[^1], "\t", po.nodes[po.node_indexes[('b',traveled_nodes[^1])]].nts
      # if record.mapping_quality == 0:
      #   total_secondary_mapping_time += (cpuTime() - time1)
      # else:
      #   total_first_mapping_time += (cpuTime() - time1)

  # echo "Primary alignments - ", total_first_mapping_time
  # echo "Secondary alignments - ", total_secondary_mapping_time
  # var outfile3 : File
  # discard open(outfile3,"illumina_weighted.html",fmWrite)
  # htmlOutput(po, outfile3)
  # echo po.weights
  po.illumina_branches = illumina_branches
  topologicalSort(po)
  # for read in po.reads:
  #   if 15813 in read.corrected_path:
  #     echo read.corrected_path

  # echo "EDGES:"
  # for u in sorted(toSeq(po.edges.keys())):
  #   echo u, " - ", po.edges[u]
  # echo "WEIGHTS:"
  # for (u,v) in sorted(toSeq(po.weights.keys())):
  #   echo u, " - ", v
  # echo po.weights
  var outfile3 : File
  discard open(outfile3,"illumina_weighted.html",fmWrite)
  htmlOutput(po[], outfile3)

proc convertPOGraphtoTrimmedPOGraph( po : var POGraph) : TrimmedPOGraph = 
  var log_likelihoods : Table[(uint32,uint32),float64]
  var node_indexes : Table[(char,uint32), uint32]
  var source_nodes : seq[uint32]
  var end_nodes : seq[uint32]
  var deleted_nodes : HashSet[uint32]
  for i in 0'u32..<uint32(po.nodes.len):
    node_indexes[('f',i)] = i
    node_indexes[('b',i)] = i
  for i in 0..<po.reads.len:
    po.reads[i].corrected_path = po.reads[i].path
  return TrimmedPOGraph( nodes : po.nodes,
                         reads : po.reads, 
                         edges : po.edges,
                         weights : po.weights,
                         og_nodes : po.og_nodes,
                         log_likelihoods : log_likelihoods,
                         node_indexes : node_indexes,
                         source_nodes : source_nodes,
                         end_nodes : end_nodes,
                         deleted_nodes : deleted_nodes )

proc buildConsensusPO( po : ptr POGraph, paths : seq[seq[uint32]], read_prefix : string) : TrimmedPOGraph = 
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  var u_set,v_set : HashSet[uint32]
  var reads : seq[Read]
  for j,path in paths:
    var sequence : seq[string]
    for i in 0..<path.len - 1:
      if (path[i],path[i+1]) in weights:
        weights[(path[i],path[i+1])] += 1'u32
      else:
        weights[(path[i],path[i+1])] = 1'u32
        if path[i] in edges:
          edges[path[i]].add(path[i+1])
        else:
          edges[path[i]] = @[path[i+1]]
        u_set.incl(path[i])
        v_set.incl(path[i+1])
      sequence.add(po[].nodes[path[i]].nts)
    sequence.add(po[].nodes[path[^1]].nts)
    reads.add(Read(name : &"{read_prefix}_{j}",
                   length : uint32(path.len),
                   path : path,
                   corrected_path : path,
                   sequence : sequence.join()))
    # for i in 0..<read.corrected_path.len:
    #   echo read.corrected_path[i]
  var source_nodes1,source_nodes2 : HashSet[uint32]
  for i,node in po[].nodes:
    if node.start_node_flag:
      source_nodes1.incl(uint32(i))
  let sink_nodes = v_set - u_set
  for sink_node in sink_nodes:
    edges[sink_node] = @[]
  var nodes : seq[Node]
  var node_indexes = sorted(toSeq(u_set + v_set))
  var node_indexes2 : Table[(char,uint32),uint32]
  var log_likelihoods : Table[(uint32,uint32),float64]
  for i,j in node_indexes:
    nodes.add(po[].nodes[j])
    if j in source_nodes1:
      source_nodes2.incl(j)
    node_indexes2[('f',uint32(i))] = j
    node_indexes2[('b',j)] = uint32(i)
    if j in edges:
      var total_read_count = 0'u32
      for k in edges[j]:
        total_read_count += weights[(j,k)]
      let tmp = float(total_read_count)
      for k in edges[j]:
        log_likelihoods[(j,k)] = - ln(float(weights[(j,k)]) / tmp)
  var end_nodes : seq[uint32]
  for i in 0..<nodes.len:
    nodes[i].indegree = 0'u16
    if nodes[i].end_node_flag:
      end_nodes.add(node_indexes2[('f',uint32(i))])
  for j in edges.keys:
    for k in edges[j]:
      nodes[node_indexes2[('b',k)]].indegree += 1'u16
  var trim = TrimmedPOGraph(nodes : nodes,
                        edges : edges,
                        weights : weights,
                        reads : reads,
                        og_nodes : po[].og_nodes,
                        log_likelihoods : log_likelihoods,
                        source_nodes : toSeq(source_nodes2),
                        end_nodes : end_nodes,
                        node_indexes : node_indexes2)
  topologicalSort(addr trim)
  return trim

proc getTrimmedGraphFromFastaRecord(record : FastaRecord) : TrimmedPOGraph = 
  var nodes : seq[Node]
  var path : seq[uint32]
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  var node_indexes : Table[(char,uint32), uint32]
  for i,base in record.sequence:
    nodes.add(Node( nts : $base,
                    supporting_reads : @[0'u32],
                    align_ring_partner : -1'i32))
    path.add(uint32(i))
    let v = uint32(i)
    node_indexes[('f',v)] = v
    node_indexes[('b',v)] = v
    if i != 0:
      let u = uint32(i - 1)
      edges[u] = @[v]
      weights[(u,v)] = 1'u32
  edges[uint32(record.sequence.len - 1)] = @[]
  # assert record.sequence.len == nodes.len
  # echo edges
  let reads = @[Read(name : record.read_id, 
                     length : uint32(record.sequence.len),
                     path : path,
                     corrected_path : path,
                     sequence : record.sequence)]
  return TrimmedPOGraph(nodes : nodes,
                        reads : reads,
                        edges : edges,
                        weights : weights,
                        og_nodes : uint32(nodes.len),
                        node_indexes : node_indexes)

proc getSequenceFromPath(po : TrimmedPOGraph, path : seq[uint32]) : string = 
  var sequence1 : seq[string]
  for idx in path:
    sequence1.add(po.nodes[po.node_indexes[('b',idx)]].nts)
  return sequence1.join()

if paramStr(1) == "trim_po":
  var infile : File
  discard open(infile,paramStr(2))
  var time1 = cpuTime()
  var test =  initPOGraph(infile)
  var time2 = cpuTime()
  echo "POGraph Initial Assembly:    ", time2 - time1

  time1 = cpuTime()
  var representative_paths = trimAndCollapsePOGraph2(test,psi = 15)
  time2 = cpuTime()
  echo "POGraph Trimming & Collapsing:    ", time2 - time1
  echo representative_paths.len
  let label = paramStr(3).split(sep=os.DirSep)[^1].split(sep='.')[0]
  let consensus_po = buildConsensusPO(addr test, representative_paths, label)

  var outfile2 : File
  discard open(outfile2,paramStr(3),fmWrite)
  writeCorrectedReads(consensus_po,outfile2)

elif paramStr(1) == "illumina":
  var time1 = cpuTime()
  var infile : File
  discard open(infile,paramStr(2))
  var test =  initPOGraph(infile)
  var time2 = cpuTime()
  var trim = convertPOGraphtoTrimmedPOGraph(test)
  var outfile2 : File
  discard open(outfile2,"testing_repo.html",fmWrite)
  htmlOutput(trim,outfile2)
  # trim = topologicalSort(trim)
  echo "POGraph Initial Assembly:    ", time2 - time1

  time1 = cpuTime()
  var bam:Bam
  discard open(bam,paramStr(3),index=true)
  illuminaPolishPOGraph( addr trim, bam, debug=true )
  time2 = cpuTime()
  echo "POGraph Illumina Polishing:    ", time2 - time1
  # var outfile3 : File
  # discard open(outfile3,"testing_repo2.html",fmWrite)
  # htmlOutput(trim2,outfile3)
  time1 = cpuTime()
  var representative_paths =  getRepresentativePaths3(trim,psi = 35)
  echo representative_paths.len
  time2 = cpuTime()
  echo "Representative Paths:          ", time2 - time1
  var outfile1 : File
  discard open(outfile1,paramStr(4),fmWrite)
  writeCorrectedReads(trim,outfile1)
  let label = paramStr(2).split(sep=os.DirSep)[^1].split(sep='.')[0]
  let consensus_po = buildConsensusPO(addr trim, representative_paths,label)
  # var outfile3 : File
  # discard open(outfile3,paramStr(5),fmWrite)
  # writeCorrectedReads(trim2,outfile3)

elif paramStr(1) == "final_illumina":
  var time1 = cpuTime()
  var infile : File
  discard open(infile,paramStr(2))
  var reads =  parseFasta(infile)
  var time2 = cpuTime()
  echo "Fasta Parse:    ", time2 - time1
  var bam:Bam
  discard open(bam,paramStr(3),index=true)
  var corrected : seq[FastaRecord]
  for read in reads:
    var trim = getTrimmedGraphFromFastaRecord(read)
    topologicalSort(addr trim)
    illuminaPolishPOGraph( addr trim, bam, debug=true)
    var representative_paths = getRepresentativePaths2(trim,psi = 35)
    corrected.add(FastaRecord(read_id : read.read_id, sequence : getSequenceFromPath(trim,trim.reads[0].corrected_path)))
  var outfile1 : File
  discard open(outfile1,paramStr(4),fmWrite)
  writeCorrectedReads(corrected,outfile1)

  # var outfile3 : File
  # discard open(outfile3,paramStr(5),fmWrite)
  # writeCorrectedReads(trim2,outfile3)

  # let consensus_po = buildConsensusPO(addr trim2, representative_paths, paramStr(4))
  # var outfile1 : File
  # discard open(outfile1,paramStr(4),fmWrite)
  # writePOGraph(consensus_po, outfile1)

  # var outfile2 : File
  # discard open(outfile2,paramStr(5),fmWrite)
  # writeCorrectedReads(consensus_po,outfile2)

