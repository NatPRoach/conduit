import strutils
import strformat
import sequtils
import times
import algorithm
import tables
import sets
import math
import heapqueue
import hts
import fasta

type
  Read* = object
    name* : string
    length* : uint32
    path* : seq[uint32]
    correctedPath* : seq[uint32]
    sequence* : string
    support* : uint32
  
  Node* = object
    nts* : string
    # supportingReads : seq[uint32]
    nanoporeSupport* : uint32
    illuminaSupport* : uint32
    illuminaSupported* : bool
    alignRingPartner* : int32 #No partner if -1
    visited* : bool
    indegree* : uint16
    logLikelihood* : float
    path* : seq[uint32]
    correctedPath* : seq[uint32]
    startNodeFlag* : bool
    endNodeFlag* : bool
  
  POGraph* = object of RootObj
    nodes* : seq[Node]
    reads* : seq[Read]
    edges* : Table[uint32,seq[uint32]]
    weights* : Table[(uint32,uint32),uint32]
    ogNodes* : uint32
  
  TrimmedPOGraph* = object of POGraph
    logLikelihoods* : Table[(uint32,uint32),float64]
    nodeIndexes* : Table[(char,uint32), uint32]
    sourceNodes* : seq[uint32]
    endNodes* : seq[uint32]
    deletedNodes* : HashSet[uint32]
    illuminaBranches* : Table[seq[uint32],seq[uint32]]
    nanoporeCounts* : Table[(uint32,uint32),uint32]
    illuminaCounts* : Table[(uint32,uint32),uint32]
    illuminaSupported* : HashSet[(uint32,uint32)]

  
  IsoformExtractionDataStructure = object
    nodes : seq[Node]
    nodeIndexes : Table[(char,uint32), uint32]
    fwdEdges : Table[uint32,seq[uint32]]
    revEdges : Table[uint32,seq[uint32]]
    weights : Table[(uint32,uint32),seq[string]]
    nodeSupportList : seq[uint32]

#TODO: clean up redundancies, remove dead code blocks
#TODO: Nim style guide variable name conformity
#TODO: Nim style guide 80 char limit per line conformity
#TODO: break up this into several .nim files, each with a more descriptive /
#TODO: accurate name
#TODO: remove unused functions

proc collapseLinearStretches( po2 : TrimmedPOGraph) : TrimmedPOGraph =
  var po = po2
  var toDelete : HashSet[uint32]
  var newEndNodes : seq[uint32]
  # for idx in sorted(toSeq(po.nodeIndexes.keys)):
  #   echo "idx, ", idx, " ", po.nodeIndexes[idx]
  for i in 0..<po.nodes.len:
    po.nodes[i].visited = false
    # if po.nodes[i].endNodeFlag:
    #   echo "flag ", i
  for i in 0..<po.nodes.len:
    if po.nodes[i].visited:
      continue
    let nodeIdx = po.nodeIndexes[('f',uint32(i))]
    if po.edges[nodeIdx].len == 1 and not po.nodes[i].endNodeFlag:
      var nextNodeIdx = po.edges[nodeIdx][0]
      var backNextNodeIdx = po.nodeIndexes[('b',nextNodeIdx)]
      while po.edges[nodeIdx].len == 1 and
            po.nodes[backNextNodeIdx].indegree == 1 and 
            not po.nodes[backNextNodeIdx].startNodeFlag and
            not po.nodes[backNextNodeIdx].endNodeFlag:
        po.nodes[i].nts = po.nodes[i].nts & po.nodes[backNextNodeIdx].nts
        po.nodes[backNextNodeIdx].visited = true
        toDelete.incl(backNextNodeIdx)
        po.edges[nodeIdx] = po.edges[nextNodeIdx]
        po.edges.del(nextNodeIdx)
        po.logLikelihoods.del((nodeIdx,nextNodeIdx))
        po.weights.del((nodeIdx,nextNodeIdx))
        if po.edges[nodeIdx].len == 1:
          let k = po.edges[nodeIdx][0]
          po.logLikelihoods[(nodeIdx,k)] = po.logLikelihoods[(nextNodeIdx,k)]
          po.logLikelihoods.del((nextNodeIdx,k))
          po.weights[(nodeIdx,k)] = po.weights[(nextNodeIdx,k)]
          po.weights.del((nextNodeIdx,k))
          nextNodeIdx = po.edges[nodeIdx][0]
          backNextNodeIdx = po.nodeIndexes[('b',nextNodeIdx)]
        elif po.edges[nodeIdx].len > 1:
          for k in po.edges[nodeIdx]:
            po.logLikelihoods[(nodeIdx,k)] = po.logLikelihoods[(nextNodeIdx,k)]
            po.logLikelihoods.del((nextNodeIdx,k))
            po.weights[(nodeIdx,k)] = po.weights[(nextNodeIdx,k)]
            po.weights.del((nodeIdx,k))
          break
      if po.nodes[backNextNodeIdx].endNodeFlag:
        if po.edges[nodeIdx].len == 1 and
           po.nodes[backNextNodeIdx].indegree == 1 and
           not po.nodes[backNextNodeIdx].startNodeFlag:
          po.nodes[i].nts = po.nodes[i].nts & po.nodes[backNextNodeIdx].nts
          po.nodes[backNextNodeIdx].visited = true
          toDelete.incl(backNextNodeIdx)
          po.edges[nodeIdx] = po.edges[nextNodeIdx]
          po.edges.del(nextNodeIdx)
          po.logLikelihoods.del((nodeIdx,nextNodeIdx))
          po.weights.del((nodeIdx,nextNodeIdx))
          if po.edges[nodeIdx].len > 1:
            for k in po.edges[nodeIdx]:
              po.logLikelihoods[(nodeIdx,k)] =
                po.logLikelihoods[(nextNodeIdx,k)]
              po.logLikelihoods.del((nextNodeIdx,k))
              po.weights[(nodeIdx,k)] = po.weights[(nextNodeIdx,k)]
              po.weights.del((nextNodeIdx,k))
          newEndNodes.add(nodeIdx)
          # echo "here ", nodeIdx
        elif po.edges[nodeIdx].len != 0:
          newEndNodes.add(nodeIdx)
          # echo "here3 ", nodeIdx
        else:
          newEndNodes.add(nextNodeIdx)
          # echo "here2 ", nextNodeIdx
    else:
      if po.nodes[i].endNodeFlag:
        newEndNodes.add(nodeIdx)
  var newNodes : seq[Node]
  var newNodeIndexes : Table[(char,uint32),uint32]
  var counter = 0'u32
  for i in 0'u32..<uint32(po.nodes.len):
    if i in toDelete:
      continue
    newNodes.add(po.nodes[i])
    newNodeIndexes[('f',counter)] = po.nodeIndexes[('f',i)]
    newNodeIndexes[('b',po.nodeIndexes[('f',i)])] = counter
    counter += 1'u32
  return TrimmedPOGraph(nodes : newNodes,
                        nodeIndexes : newNodeIndexes,
                        edges: po.edges,
                        weights : po.weights,
                        reads : po.reads,
                        logLikelihoods : po.logLikelihoods,
                        sourceNodes : po.sourceNodes,
                        endNodes : newEndNodes)

proc jsOutput( po:TrimmedPOGraph,
               highlight_path1,highlight_path2 : seq[uint32] = @[],
               collapse : bool = false) : seq[string] = 
  # Logic for jsOutput and htmlOutput procs largely borrowed, modified, and
  # converted from the Simpson lab simple Partial Order Alignment python
  # implementation: https://github.com/ljdursi/poapy
  var path : seq[uint32]
  if highlight_path1.len != 0:
    path = highlight_path1
  else:
    path = po.reads[0].correctedPath
  
  var highlightNodes1,highlightNodes2 : HashSet[uint32]
  var highlightEdges1,highlightEdges2 : HashSet[(uint32,uint32)]
  for i,node in highlight_path1:
    highlightNodes1.incl(node)
    if i != 0:
      highlightEdges1.incl((highlight_path1[i-1],node))
  for i,node in highlight_path2:
    highlightNodes2.incl(node)
    if i != 0:
      highlightEdges2.incl((highlight_path2[i-1],node))

  var pathTable : Table[uint32,int]
  for i, nodeID in path:
    pathTable[nodeID] = i * 150
  var lines = @["var nodes = ["]
  var count = 0
  var lastJInPathTable = 999999
  let f = '{'
  let e = '}'
  for i in pathTable.keys():
    if int(i) < lastJInPathTable:
      lastJInPathTable = int(i)
  for i in 0..<po.nodes.len:
    let j = uint32(i)
    if j in po.deletedNodes:
      continue
    let node = po.nodes[po.nodeIndexes[('b',uint32(i))]]
    var line : seq[string]
    if collapse:
      line = @[&"    {f}id: {j}, label: \"{j}\""]
    else:
      line = @[&"    {f}id: {j}, label: \"{j}, {node.nts}\""]
    var highlight : string
    if j in highlightNodes1:
      if j in highlightNodes2:
        highlight = "color: \'purple\', "
      else:
        highlight = "color: \'red\', "
    elif j in highlightNodes2:
      highlight = "color: \'blue\', "
    # elif node.visited:
    #   highlight = "color: \'#008000\', "
    else:
      highlight = "color: \'#D3D3D3\', "
    if (j in pathTable):
      if((count %% 5) == 0):
        line.add(&", allowedToMoveX: true, x: {pathTable[j]} , y: 0, " &
          &"{highlight} allowedToMoveY: true{e},")
        lastJInPathTable = int(j)
      else:
        line.add(&", allowedToMoveX: true, x: " &
          &"{pathTable[uint32(lastJInPathTable)]} , y: 0, {highlight} " &
          &"allowedToMoveY: true{e},")
      count += 1
    else:
      line.add(&", allowedToMoveX: true, x: " &
        &"{pathTable[uint32(lastJInPathTable)]} , y: 0, {highlight} " &
        &"allowedToMoveY: true{e},")
    lines.add(line.join(""))
  lines[^1] = lines[^1][0..^1]
  lines.add("];")
  lines.add(" ")
  lines.add("var edges = [")

  var writtenSet : HashSet[(uint32,uint32)]
  for i in 0..<po.nodes.len:
    let k = po.nodeIndexes[('f',uint32(i))]
    if k in po.deletedNodes:
      continue
    for edge in po.edges[k]:
      if edge in po.deletedNodes:
        continue
      var highlight : string
      if not((k,edge) in writtenSet):
        if (k,edge) in highlightEdges1:
          if (k,edge) in highlightEdges2:
            highlight = ", color: \'purple\'"
          else:
            highlight = ", color: \'red\'"
        elif (k,edge) in  highlightEdges2:
          highlight = ", color: \'blue\'"
        else:
          highlight = ", color: \'black\'"
        lines.add(&"    {f}from: {k}, to: {edge}{highlight}, label: " &
          &"'{po.weights[(k,edge)]}',     font: {f}align: 'middle'{e}," &
          &" value: {po.weights[(k,edge)]}, arrows:\'to\'{e},")
        writtenSet.incl((k,edge))
    if po.nodes[i].alignRingPartner != -1:
      lines.add(&"    {f}from: {k}, to: {po.nodes[i].alignRingPartner}," &
        &"color: \'grey\'{e},")
  lines[^1] = lines[^1][0..^1]
  lines.add("];\n")
  return lines

proc htmlOutput*( po:TrimmedPOGraph,
                 outfile : File,
                 highlight_path1,highlight_path2 : seq[uint32] = @[],
                 collapse = false) = 
  let open = '{'
  let close = '}'
  let header = "<!doctype html>\n<html>\n<head>\n<title>POA Graph Alignment" &
    "</title>\n<script type=\"text/javascript\" src=\"https://visjs.github." &
    "io/vis-network/standalone/umd/vis-network.min.js\"></script>\n<script " &
    "type=\"text/javascript\" src=\"https://cdnjs.cloudflare.com/ajax/libs/" &
    "jspdf/1.5.3/jspdf.min.js\"></script>\n</head>\n<body>\n<button id=\"do" &
    "wnloadPDF\">Download PDF</button>\n<div id=\"mynetwork\"></div>\n<scri" &
    "pt type=\"text/javascript\">\n// create a network\n"
  outfile.write(header)
  let header2 =  &"    document.getElementById(\"downloadPDF\").addEventList" &
    &"ener(\"click\", () => {open}\n        const canvas = document.getEleme" &
    &"ntsByTagName(\"canvas\")[0];\n        const myImage = canvas.toDataURL" &
    &"(\"image/png,1.0\");\n\n//Adjust width and height\n        const imgWi" &
    &"dth = (canvas.width * 20) / 240;\n        const imgHeight = (canvas.he" &
    &"ight * 20) / 240;\n\n        const pdf = new jsPDF('p', 'mm', 'a4');\n" &
    &"        pdf.addImage(myImage, 'PNG', 15, 2, imgWidth, imgHeight); // 2" &
    &": 19\n        pdf.save('Download.pdf');\n    {close});"
  outfile.write(header2)
  var lines :seq[string]
  if collapse:
    lines = jsOutput(collapseLinearStretches(po),
                     highlight_path1,
                     highlight_path2,
                     collapse)
  else:
    lines = jsOutput(po,
                     highlight_path1,
                     highlight_path2)
  outfile.write(lines.join("\n"))
  let footer = "var container = document.getElementById('mynetwork');\nvar " &
    "data= {\nnodes: nodes,\nedges: edges,\n};\nvar options = {\nwidth: '10" &
    "0%',\nheight: '800px'\n};\nvar network = new vis.Network(container, da" &
    "ta, options);\n</script>\n</body>\n</html>\n"
  outfile.write(footer)
  outfile.close()

proc writeCorrectedReads*( po : TrimmedPOGraph,
                           outfile : File) = 
  for read in po.reads:
    var sequence : seq[string]
    for idx in read.correctedPath:
      sequence.add(po.nodes[po.nodeIndexes[('b',idx)]].nts)
    outfile.write(">",read.name,&"_{read.support}","\n")
    outfile.writeLine(sequence.join())

# proc writeCorrectedReads*( records : seq[FastaRecord],outfile : File) = 
#   for record in records:
#     outfile.write(">",record.readId,"\n")
#     outfile.writeLine(record.sequence)

proc initPOGraph*( file : File) : POGraph = 
  for _ in 0..2:
    discard file.readLine()
  let numNodes = uint32(parseUInt(file.readLine().split('=')[1]))
  let numReads = uint32(parseUInt(file.readLine().split('=')[1]))
  var reads : seq[Read]
  var nodes : seq[Node]
  var sourceNodeIndices : seq[uint32]
  var endNodeIndices : seq[uint32]
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  for _ in 0..<numReads:
    let name = file.readLine().strip().split('=')[1]
    let length = uint32(
                   parseUInt(
                     file.readLine().strip().split('=')[1].split(' ')[0]))
    reads.add(Read(name:name,length:length))
  for i in 0..<numNodes:
    let tmp = file.readLine().strip().split(':')
    let nt = tmp[0]
    var info = tmp[1]
    let aPos = info.find("A")
    var aPair : int32
    if aPos != -1:
      aPair = int32(parseUInt(info[aPos+1..^1]))
      info = info[0..<aPos]
    else:
      aPair = -1
    nodes.add(Node(nts : nt,
                   alignRingPartner : aPair,
                   visited : false,
                   endNodeFlag : false,
                   startNodeFlag : false,
                   nanoporeSupport : 0'u32))
    let sSplit = info.split('S')
    info = sSplit[0]
    var seqIdxs : seq[uint32]
    for s in sSplit[1..^1]:
      seqIdxs.add(uint32(parseUInt(s)))
    for seq_idx in seqIdxs:
      if reads[seq_idx].path.len == 0:
        nodes[i].startNodeFlag = true
        sourceNodeIndices.add(uint32(i))
      reads[seq_idx].path.add(uint32(i))
      reads[seq_idx].sequence.add(nt)
      # nodes[i].supportingReads.add(seq_idx)
      nodes[i].nanoporeSupport += 1'u32
    if info != "":
      let lSplit = info.split('L')
      for l in lSplit[1..^1]:
        let inNode = uint32(parseUInt(l))
        if inNode in edges:
          edges[inNode].add(uint32(i))
        else:
          edges[inNode] = @[uint32(i)]
        weights[(uint32(inNode),uint32(i))] = 0'u32
  for read in reads:
    # nodes[read.path[^1]].endNodeFlag = true
    endNodeIndices.add(read.path[^1])
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
  return POGraph(nodes:nodes,
                 reads:reads,
                 edges:edges,
                 ogNodes:uint32(nodes.len),
                 weights:weights)


proc writePOGraph*( po : TrimmedPOGraph,
                    outfile : File,
                    graphname : string = "default") = 
  outfile.write("VERSION=UNTITLEDCORRECTIONALGORITHM.0.1\n")
  outfile.write(&"NAME={graphname}\n")
  outfile.write("TITLE=untitled\n")
  outfile.write(&"LENGTH={po.nodes.len}\n")
  outfile.write(&"SOURCECOUNT={po.reads.len}\n")
  for read in po.reads:
    outfile.write(&"SOURCENAME={read.name}\n")
    outfile.write(&"SOURCEINFO={read.correctedPath.len} 0 1 -1 untitled\n")
  var reverseEdges : Table[uint32,seq[uint32]]
  for u in po.edges.keys():
    for v in po.edges[u]:
      if v in reverseEdges:
        reverseEdges[v].add(u)
      else:
        reverseEdges[v] = @[u]
  var supportingReadsTable : Table[uint32,seq[uint32]]
  for i,read in po.reads:
    # echo read.correctedPath
    for node in read.correctedPath:
      if not(node in supportingReadsTable):
        supportingReadsTable[node] = @[uint32(i)]
      else:
        supportingReadsTable[node].add(uint32(i))
  for i,node in po.nodes:
    var edges : seq[string]
    if po.nodeIndexes[('f',uint32(i))] in reverseEdges:
      for u in reverseEdges[po.nodeIndexes[('f',uint32(i))]]:
        edges.add(&"L{po.nodeIndexes[('b',u)]}")
    var supportingReads : seq[string]
    for read_num in supportingReadsTable[po.nodeIndexes[('f',uint32(i))]]:
      supportingReads.add(&"S{read_num}")
    outfile.write(&"{node.nts}:{edges.join()}{supportingReads.join()}\n")
  outfile.close()

proc stringencyCheck*(po : ptr TrimmedPOGraph,
                      path : seq[uint32],
                      stringent_tolerance : int = 100) : bool = 
  result = true
  for i in stringent_tolerance..<path.len-stringent_tolerance:
    let u = path[i-1]
    let v = path[i]
    let bckNodeIdx = po[].nodeIndexes[('b',v)]
    result = 
      result and
      po[].nodes[bckNodeIdx].illuminaSupported  and
      ((u,v) in po[].illuminaSupported)
  # result =
  #   result and
  #   ((path[^(stringent_tolerance + 2)],
  #     path[^(stringent_tolerance + 1)]) in po[].total_illuminaSupport)


proc getIsoformExtractionDataStructures( po : ptr POGraph) :
                                           IsoformExtractionDataStructure =
  # var time1 = cpuTime()
  var fwdEdges,revEdges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),seq[string]]
  var uSet,vSet : HashSet[uint32]
  var nodeSupportCounts : CountTable[uint32]
  for read in po[].reads:
    for i in 0..<read.correctedPath.len-1:
      if (read.correctedPath[i],read.correctedPath[i+1]) in weights:
        weights[(read.correctedPath[i],read.correctedPath[i+1])].add(read.name)
      else:
        weights[(read.correctedPath[i],read.correctedPath[i+1])] = @[read.name]
        if read.correctedPath[i] in fwdEdges:
          fwdEdges[read.correctedPath[i]].add(read.correctedPath[i+1])
        else:
          fwdEdges[read.correctedPath[i]] = @[read.correctedPath[i+1]]
        if read.correctedPath[i+1] in revEdges:
          revEdges[read.correctedPath[i+1]].add(read.correctedPath[i])
        else:
          revEdges[read.correctedPath[i+1]] = @[read.correctedPath[i]]
        uSet.incl(read.correctedPath[i])
        vSet.incl(read.correctedPath[i+1])
      if read.correctedPath[i] in nodeSupportCounts:
        nodeSupportCounts.inc(read.correctedPath[i])
      else:
        nodeSupportCounts[read.correctedPath[i]] = 1
    if read.correctedPath[^1] in nodeSupportCounts:
      nodeSupportCounts.inc(read.correctedPath[^1])
    else:
      nodeSupportCounts[read.correctedPath[^1]] = 1
  # echo "init dicts - ", (cpuTime() - time1)
  # time1 = cpuTime()
  let sourceNodes = uSet - vSet
  let sinkNodes = vSet - uSet
  for sink_node in sinkNodes:
    fwdEdges[sink_node] = @[]
  for source_node in sourceNodes:
    revEdges[source_node] = @[]
  # var nodes : seq[Node]
  let nodeIndexes = sorted(toSeq(uSet + vSet))
  # echo "node sorting - ", (cpuTime() - time1)
  # time1 = cpuTime()
  var nodeIndexes2 : Table[(char,uint32),uint32]
  var nodeSupportList : seq[uint32]
  for i,j in nodeIndexes:
    # nodes.add(po[].nodes[j])
    nodeIndexes2[('f',uint32(i))] = j
    nodeIndexes2[('b',j)] = uint32(i)
    nodeSupportList.add(uint32(nodeSupportCounts[j]))
  # echo "nodeIndexes - ", cpuTime() - time1

  return IsoformExtractionDataStructure(nodeIndexes : nodeIndexes2,
                                        fwdEdges : fwdEdges,
                                        revEdges : revEdges,
                                        weights : weights,
                                        nodeSupportList: nodeSupportList)

proc updateIsoformExtractionDataStructures(
    ie : ptr IsoformExtractionDataStructure,
    removedReads : ptr seq[Read]) = 
  var removedReadNames : HashSet[string]
  var visitedEdges : HashSet[(uint32,uint32)]
  for read in removedReads[]:
    removedReadNames.incl(read.name)
  for read in removedReads[]:
    for i in 0..<read.correctedPath.len-1:
      let u = read.correctedPath[i]
      let v = read.correctedPath[i+1]
      ie[].nodeSupportList[ie[].nodeIndexes[('b',u)]] -= 1'u32
      if (u,v) in visitedEdges:
        continue
      visitedEdges.incl((u,v))
      var newReadNames : seq[string]
      for name in ie[].weights[(u,v)]:
        if name in removedReadNames:
          continue
        newReadNames.add(name)
      if newReadNames.len != 0:
        ie[].weights[(u,v)] = newReadNames
      else:
        ie[].weights.del((u,v))
        for i,v2 in ie[].fwdEdges[u]:
          if v2 == v:
            ie[].fwdEdges[u].del(i)
            break
        for i,u2 in ie[].revEdges[v]:
          if u2 == u:
            ie[].revEdges[v].del(i)
            break
    ie[].nodeSupportList[ie[].nodeIndexes[('b',
                                           read.correctedPath[^1])]] -= 1'u32

proc maxIdx[T]( s : seq[T]) : int = 
  var currentMax = T.low
  var currentMaxIdx = 0
  for i,j in  s:
    if j > currentMax:
      currentMax = j
      currentMaxIdx = i
  return currentMaxIdx


proc alignPaths( x : seq[uint32],
                 y : seq[uint32]) : seq[uint8] = 
  ### Takes in two paths and aligns them, leveraging the fact that the nodes
  ### are labeled in topological sort order to speed alignment
  ### x - first path
  ###     query for the purposes of defining insertion / deletion
  ### y - second path,
  ###     reference for the purposes of defining insertion / deletion
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
    else:  # x[i] > y[j]
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
                    nodeIndexes : Table[(char,uint32),uint32],
                    psi:uint16,
                    correct_ends : bool = false,
                    detect_ends : bool = false,
                    debug : bool = false) : (bool,seq[uint32]) =
  ### Takes in two paths, a query path and a reference path, and a number psi
  ### that indicates the min indel size needed to retain that indel, indels in
  ### the query relative to the reference smaller than that number are 
  ### 'corrected', to reflect the reference path. Also takes in a bool 
  ### indicating whether or not to correct the ends of alignments, a reference 
  ### to nodes, which get modified to reflect new start and end nodes based on 
  ### correction.
  for i in 0..<q_path.len - 1:
    # echo q_path[i], " ", q_path[i+1]
    assert q_path[i] < q_path[i+1]
  for i in 0..<r_path.len - 1:
    assert r_path[i] < r_path[i+1]
  let alignment = alignPaths(q_path,r_path)
  var qIdx = 0
  var rIdx = 0
  var continuousIns = 0
  var continuousDel = 0
  var diffFlag = false
  var endDiffFlag = false
  var qOps : seq[(char,int,int)]
  var startFlag = true
  var ignoreDeletionFlag = false
  for op in alignment:
    if op == 0'u8:
      if debug:
        echo q_path[qIdx], "\t", r_path[rIdx]
      if continuousIns >= int(psi):
        diffFlag = true
        ignoreDeletionFlag = true
      elif continuousIns > 0 and (correct_ends or not startFlag):
        for i in qIdx - continuousIns..<qIdx:
          qOps.add(('i',i,0))
      if continuousDel >= int(psi):
        diffFlag = true
        endDiffFlag = true
      elif continuousDel > 0 and not ignoreDeletionFlag:
        for i in rIdx - continuousDel..<rIdx:
          qOps.add(('d',qIdx-continuousIns,i))
      ignoreDeletionFlag = false
      startFlag = false
      qIdx += 1
      rIdx += 1
      continuousIns = 0
      continuousDel = 0
    elif op == 1'u8: # Insertion
      if debug:
        echo q_path[qIdx], "\t----"
      if detect_ends or correct_ends or not startFlag:
        continuousIns += 1

      qIdx += 1
    else: #op == 2; Deletion
      if debug:
        echo "----\t", r_path[rIdx]
      if not startFlag:
        continuousDel += 1
      rIdx += 1
  if detect_ends:
    if continuousIns >= int(psi):
      diffFlag = true
  if correct_ends:
    if continuousIns >= int(psi):
      diffFlag = true
    elif continuousIns > 0:
      for i in qIdx - continuousIns..<qIdx:
        qOps.add(('i',i,0))
  
  assert qIdx == q_path.len
  assert rIdx == r_path.len
  # let before = q_path
  for (op,qIdx,rIdx) in sorted(qOps,compareQopTuples,Descending):
    if op == 'i':
      if qIdx != len(q_path) - 1 and qIdx != 0:
        assert q_path[qIdx - 1] < q_path[qIdx + 1]
      q_path.delete(qIdx)
    else: # op == 'd'
      if qIdx != 0:
        assert q_path[qIdx - 1] < r_path[rIdx] and r_path[rIdx] < q_path[qIdx]
      q_path.insert(@[r_path[rIdx]],qIdx)
  if debug:
    echo q_path.len
  return (diffFlag,q_path)

proc walkHeaviestPaths( po : ptr POGraph,
                        ends_delta = 15) : (seq[seq[uint32]],seq[uint32]) = 
  var representativePaths : seq[seq[uint32]]
  var readSupports : seq[uint32]
  var remainingReads = toSeq(0..<po[].reads.len)
  # var total_st_time = 0.0
  # var total_sort_time = 0.0
  # var total_walk_time = 0.0
  # var total_search_time = 0.0
  # let time = cpuTime()
  var st = getIsoformExtractionDataStructures(po)
  # total_st_time += (cpuTime() - time)
  # echo "ST init time - ", total_st_time
  # total_st_time = 0.0
  # var iterations = 0 
  var removedReads : seq[Read]
  while remainingReads.len != 0:
    # var time1 = cpuTime()
    updateIsoformExtractionDataStructures(addr st,
                                          addr removedReads)
    # total_st_time += (cpuTime() - time1)
    # echo "here"
    # time1 = cpuTime()
    var reads : seq[Read] = @[]
    for i in remainingReads:
      reads.add(po[].reads[i])
    # let st2 = getIsoformExtractionDataStructures(po,reads)
    # echo st2.nodeSupportList
    # echo st.nodeSupportList
    # assert st2.weights == st.weights
    # assert st2.nodeSupportList == st.nodeSupportList
    let maxNodeIdx = maxIdx(st.nodeSupportList)
    let fwdMaxNodeIdx = st.nodeIndexes[('f',uint32(maxNodeIdx))]
    # echo fwdMaxNodeIdx
    var startIdx,endIdx : int
    if st.fwdEdges[fwdMaxNodeIdx].len == 0:
      endIdx = -1
    else:
      endIdx = maxNodeIdx
    if st.revEdges[fwdMaxNodeIdx].len == 0:
      startIdx = -1
    else:
      startIdx = maxNodeIdx
    
    
    var readIdxs = newSeqWith(po[].reads.len, [0,0])
    var excludedReads : HashSet[string]
    var allReads : HashSet[string]
    var fwdEndReads : HashSet[string]
    var revEndReads : HashSet[string]
    var readIdToIdx : Table[string,int]
    var readIdxToId : Table[int,string]
    for i,read in reads:
      readIdToIdx[read.name] = i
      readIdxToId[i] = read.name
      allReads.incl(read.name)
      let j = binarySearch(read.correctedPath, fwdMaxNodeIdx)
      readIdxs[i] = [j,j]
      # If not found it'll be [-1,-1]
      # else it'll be the index of the matched node
      if j == -1:
        excludedReads.incl(read.name)
      if j < ends_delta:
        revEndReads.incl(read.name)
      if read.correctedPath.len - j < ends_delta:
        fwdEndReads.incl(read.name)
    var representativePath = @[fwdMaxNodeIdx]
    # total_search_time += (cpuTime() - time1)
    # time1 = cpuTime()
    while startIdx != -1 or endIdx != -1:
      # echo "here2"
      # echo "start - ", startIdx
      # echo "end   - ", endIdx
      var nextNodeId : uint32
      if startIdx == -1:
        let fwdEndIdx = st.nodeIndexes[('f',uint32(endIdx))]
        if st.fwdEdges[fwdEndIdx].len != 1:
          # Reached end of movement in reverse direction
          # look in forward direction
          var fwdWeights : seq[uint32]
          for fwd_node in st.fwdEdges[fwdEndIdx]:
            var weight = 0'u32
            for readId in st.weights[(fwdEndIdx,fwd_node)]:
              if (readId notin excludedReads) and (readId notin fwdEndReads):
                weight += 1'u32
            fwdWeights.add(weight)
          if max(fwdWeights) == 0 and fwdEndReads.len != 0:
            fwdWeights = @[]
            for fwd_node in st.fwdEdges[fwdEndIdx]:
              var weight = 0'u32
              for readId in st.weights[(fwdEndIdx,fwd_node)]:
                if readId notin excludedReads:
                  weight += 1'u32
              fwdWeights.add(weight)
            if max(fwdWeights) == 0:
              break
          nextNodeId = st.fwdEdges[fwdEndIdx][maxIdx(fwdWeights)]
        else:
          nextNodeId = st.fwdEdges[fwdEndIdx][0]

        for i in 0..<reads.len:
          let readIdx = readIdxs[i][1]
          if readIdx == -1:
            continue
          if readIdx + 1 != len(reads[i].correctedPath):
            if reads[i].correctedPath[readIdx + 1] == nextNodeId:
              readIdxs[i][1] += 1
              if (reads[i].correctedPath.len - readIdxs[i][1]) < ends_delta:
                fwdEndReads.incl(reads[i].name)
                readIdxs[i][1] = -1
            else:
              readIdxs[i][1] = -1
              # if true:
              if reads[i].name notin fwdEndReads:
                excludedReads.incl(reads[i].name)
          else:
            readIdxs[i][1] = -1
        endIdx = int(st.nodeIndexes[('b',nextNodeId)])
        # echo "back node", endIdx
        # echo "end - ", endIdx
        if st.fwdEdges[nextNodeId].len == 0:
          endIdx = -1
      
      elif endIdx == -1:
        # Reached end of movement in forward direction
        # look in reverse direction
        let fwdStartIdx = st.nodeIndexes[('f',uint32(startIdx))]
        if st.revEdges[fwdStartIdx].len != 1:
          var revWeights : seq[uint32]
          for rev_node in st.revEdges[fwdStartIdx]:
            var weight = 0'u32
            for readId in st.weights[(rev_node,fwdStartIdx)]:
              if (readId notin excludedReads) and (readId notin revEndReads):
                weight += 1'u32
            revWeights.add(weight)
          if max(revWeights) == 0 and revEndReads.len != 0:
            revWeights = @[]
            for rev_node in st.revEdges[fwdStartIdx]:
              var weight = 0'u32
              for readId in st.weights[(rev_node,fwdStartIdx)]:
                if (readId notin excludedReads):
                  weight += 1'u32
              revWeights.add(weight)
            if max(revWeights) == 0:
              break
          nextNodeId = st.revEdges[fwdStartIdx][maxIdx(revWeights)]
        else:
          nextNodeId = st.revEdges[fwdStartIdx][0]
        for i in 0..<reads.len:
          let readIdx = readIdxs[i][0]
          if readIdx == -1:
            continue
          if readIdx != 0:
            if reads[i].correctedPath[readIdx - 1] == nextNodeId:
              readIdxs[i][0] = readIdx - 1
              if readIdxs[i][0] < ends_delta:
                revEndReads.incl(reads[i].name)
                readIdxs[i][0] = -1
            else:
              readIdxs[i][0] = -1
              # if true:
              if reads[i].name notin revEndReads:
                excludedReads.incl(reads[i].name)
          else:
            readIdxs[i][0] = -1
        startIdx = int(st.nodeIndexes[('b',nextNodeId)])
        # echo "start - ", startIdx
        if st.revEdges[nextNodeId].len == 0:
          startIdx = -1
      
      else:
        # Decide whether to move forward of reverse based on max weight in 
        # next step on the path
        let fwdStartIdx = st.nodeIndexes[('f',uint32(startIdx))]
        let fwdEndIdx = st.nodeIndexes[('f',uint32(endIdx))]
        var fwdFlag = false
        if st.revEdges[fwdStartIdx].len == 1:
          nextNodeId = st.revEdges[fwdStartIdx][0]
        elif st.fwdEdges[fwdEndIdx].len == 1:
          nextNodeId = st.fwdEdges[fwdEndIdx][0]
          fwdFlag = true
        else:
          var revWeights,fwdWeights : seq[uint32]
          for rev_node in st.revEdges[fwdStartIdx]:
            var weight = 0'u32
            for readId in st.weights[(rev_node,fwdStartIdx)]:
              if (readId notin excludedReads) and (readId notin revEndReads):
                weight += 1'u32
            revWeights.add(weight)
          
          for fwd_node in st.fwdEdges[fwdEndIdx]:
            var weight = 0'u32
            for readId in st.weights[(fwdEndIdx,fwd_node)]:
              if (readId notin excludedReads) and (readId notin fwdEndReads):
                weight += 1'u32
            fwdWeights.add(weight)
          
          var maxRevWeights = max(revWeights)
          if maxRevWeights == 0'u32 and revEndReads.len != 0:
            revWeights = @[]
            for rev_node in st.revEdges[fwdStartIdx]:
              var weight = 0'u32
              for readId in st.weights[(rev_node,fwdStartIdx)]:
                if readId notin excludedReads:
                  weight += 1'u32
              revWeights.add(weight)
            maxRevWeights = max(revWeights)

          var maxFwdWeights = max(fwdWeights)
          if maxFwdWeights == 0'u32 and fwdEndReads.len != 0:
            fwdWeights = @[]
            for fwd_node in st.fwdEdges[fwdEndIdx]:
              var weight = 0'u32
              for readId in st.weights[(fwdEndIdx,fwd_node)]:
                if readId notin excludedReads:
                  weight += 1'u32
              fwdWeights.add(weight)
            maxFwdWeights = max(fwdWeights)
          if maxRevWeights < maxFwdWeights:
            # Move in fwd direction
            if maxFwdWeights == 0:
              break
            nextNodeId = st.fwdEdges[fwdEndIdx][maxIdx(fwdWeights)]
            fwdFlag = true
          else:
            # Move in rev direction
            if maxRevWeights == 0:
              break
            nextNodeId = st.revEdges[fwdStartIdx][maxIdx(revWeights)]
        if fwdFlag:
          for i in 0..<reads.len:
            let readIdx = readIdxs[i][1]
            if readIdx == -1:
              continue
            if readIdx + 1 != reads[i].correctedPath.len:
              if reads[i].correctedPath[readIdx + 1] == nextNodeId:
                readIdxs[i][1] = readIdx + 1
                if (reads[i].correctedPath.len - readIdxs[i][1]) < ends_delta:
                  readIdxs[i][1] = -1
                  fwdEndReads.incl(reads[i].name)
              else:
                readIdxs[i][1] = -1
                # if true:
                if reads[i].name notin fwdEndReads:
                  excludedReads.incl(reads[i].name)
            else:
              readIdxs[i][1] = -1
          endIdx = int(st.nodeIndexes[('b',nextNodeId)])
          # echo "end - ", endIdx
          if st.fwdEdges[nextNodeId].len == 0:
            endIdx = -1
        else:
          for i in 0..<reads.len:
            let readIdx = readIdxs[i][0]
            if readIdx == -1:
              continue
            if readIdx != 0:
              if reads[i].correctedPath[readIdx - 1] == nextNodeId:
                readIdxs[i][0] = readIdx - 1
                if readIdxs[i][0] < ends_delta:
                  readIdxs[i][0] = -1
                  revEndReads.incl(reads[i].name)
              else:
                readIdxs[i][0] = -1
                # if true:
                if reads[i].name notin revEndReads:
                  excludedReads.incl(reads[i].name)
            else:
              readIdxs[i][0] = -1
          startIdx = int(st.nodeIndexes[('b',nextNodeId)])
          # echo "start - ", startIdx
          if st.revEdges[nextNodeId].len == 0:
            startIdx = -1
      representativePath.add(nextNodeId)
    # total_walk_time += (cpuTime() - time1)
    let accountedForReads = allReads - excludedReads
    # time1 = cpuTime()
    representativePaths.add(sorted(representativePath))
    # total_sort_time += (cpuTime() - time1)
    var toDeleteReads : seq[int]
    removedReads = @[]
    for readId in accountedForReads:
      toDeleteReads.add(readIdToIdx[readId])
      removedReads.add(po[].reads[remainingReads[readIdToIdx[readId]]])
      # echo readIdToIdx[readId]
      # echo "here"
    var support = 0'u32
    for read in removedReads:
      support += read.support
    readSupports.add(support)
    for idx in sorted(toDeleteReads,order = SortOrder.Descending):
      remainingReads.delete(idx)
  # echo "ST update time - ", total_st_time
  # echo "Search time - ",total_search_time
  # echo "Walk time - ", total_walk_time
  # echo "sort time - ", total_sort_time

  return (representativePaths,readSupports)

proc updateGraph( rep_po : ptr TrimmedPOGraph,
                  newPath,old_path : seq[uint32],
                  weight:uint32 = 1'u32) =
  ## NEW PATH
  # var time1 = cpuTime()
  rep_po.nodes[rep_po.nodeIndexes[('b',newPath[0])]].nanoporeSupport += 1'u32
  for i in 1..<newPath.len:
    let u = newPath[i-1]
    let v = newPath[i]
    if (u,v) in rep_po[].weights:
      rep_po[].weights[(u,v)] += weight
    else:
      rep_po[].weights[(u,v)] = weight
      rep_po[].edges[u].add(v)
    if (u,v) in rep_po.nanoporeCounts:
      rep_po.nanoporeCounts[(u,v)] += 1'u32
    else:
      rep_po.nanoporeCounts[(u,v)] = 1'u32
    rep_po.nodes[rep_po.nodeIndexes[('b',v)]].nanoporeSupport += 1'u32
  # echo "NEW PATH RESOLUTION - ", cpuTime() - time1
  ## OLD PATH
  # time1 = cpuTime()
  rep_po.nodes[rep_po.nodeIndexes[('b',
                                   old_path[0])]].nanoporeSupport -= 1'u32
  if rep_po.nodes[rep_po.nodeIndexes[('b',
                                      old_path[0])]].nanoporeSupport == 0'u32:
    rep_po.deletedNodes.incl(old_path[0])
  for i in 1..<old_path.len:
    let u = old_path[i-1]
    let v = old_path[i]
    rep_po[].weights[(u,v)] -= weight
    if rep_po[].weights[(u,v)] == 0'u32:
      rep_po[].weights.del((u,v))
      for i,j in rep_po[].edges[u]:
        if j == v:
          rep_po[].edges[u].delete(i,i)
          break
    rep_po.nanoporeCounts[(u,v)] -= 1'u32
    if rep_po.nanoporeCounts[(u,v)] == 0'u32:
      rep_po.nanoporeCounts.del((u,v))
      for j,v2 in rep_po.edges[u]:
        if v == v2:
          rep_po.edges[u].del(j)
          break
    rep_po.nodes[rep_po.nodeIndexes[('b',v)]].nanoporeSupport -= 1'u32
    if rep_po.nodes[rep_po.nodeIndexes[('b',v)]].nanoporeSupport == 0'u32:
      rep_po.deletedNodes.incl(v)
  # echo "OLD PATH RESOLUTION - ", cpuTime() - time1


proc walkGreedyPaths(po : ptr TrimmedPOGraph,
                     psi : uint16 = 15'u16) : seq[seq[uint32]] = 
  # var time1 = cpuTime()
  for i in 0..<po[].nodes.len:
    po[].nodes[i].visited = false
  var sourceNodes : seq[uint32]
  for read in po[].reads:
    sourceNodes.add(read.correctedPath[0])
  var greedyWalks : seq[seq[uint32]]
  for i,fwdNodeIdx in sorted(sourceNodes):
    var path = @[fwdNodeIdx]
    var bckNodeIdx = po[].nodeIndexes[('b',fwdNodeIdx)]
    if po[].nodes[bckNodeIdx].visited:
      continue
    po[].nodes[bckNodeIdx].visited = true
    var fwdNodeIdx2 = fwdNodeIdx
    while po[].edges[fwdNodeIdx2].len != 0:
      if po[].edges[fwdNodeIdx2].len == 1:
        fwdNodeIdx2 = po[].edges[fwdNodeIdx2][0]
      else:
        var weightList : seq[uint32]
        for j in po[].edges[fwdNodeIdx2]:
          weightList.add(uint32(po[].weights[(fwdNodeIdx2,uint32(j))]))
        fwdNodeIdx2 = po[].edges[fwdNodeIdx2][maxIdx(weightList)]
      path.add(fwdNodeIdx2)
      let bckNodeIdx2 = po[].nodeIndexes[('b',fwdNodeIdx2)]
      if po[].nodes[bckNodeIdx2].visited:
        for _ in 0'u16..<psi:
          if po[].edges[fwdNodeIdx2].len > 1:
            var weightList : seq[uint32]
            for j in po[].edges[fwdNodeIdx2]:
              weightList.add(uint32(po[].weights[(fwdNodeIdx2,uint32(j))]))
            # if fwdNodeIdx2 == 16939'u32:
            #   echo weightList
            #   echo maxIdx(weightList)
            #   echo po[].edges[fwdNodeIdx2]
            #   echo po[].edges[fwdNodeIdx2][maxIdx(weightList)]
            fwdNodeIdx2 = po[].edges[fwdNodeIdx2][maxIdx(weightList)]
          elif po[].edges[fwdNodeIdx2].len == 1:
            fwdNodeIdx2 = po[].edges[fwdNodeIdx2][0]
          else:
            break
          path.add(fwdNodeIdx2)
        break
      po[].nodes[bckNodeIdx2].visited = true
    greedyWalks.add(path)
    # echo "greedy1 - ", path
  for i in 0..<po[].nodes.len:
    po[].nodes[i].visited = false
  # echo "walk greedy paths - ", cpuTime() - time1
  return greedyWalks


proc detectMinipathDeltas( po : ptr TrimmedPOGraph,
                           minipath,greedy_path : seq[uint32],
                           psi : uint16 = 35) : bool =
  var i = 0
  var j = 0
  var alignment : seq[uint8]
  var alignIdxs : seq[uint32]
  var startFlag = false
  var deltaFlag = false
  # let new_minipath = searchForGreedyPath(po,minipath,greedy_path,psi)
  if minipath[0] in po.deletedNodes or minipath[^1] in po.deletedNodes:
    for i in 1..<minipath.len:
      let u = minipath[i-1]
      let v = minipath[i]
      #TODO: Test if all these comparisons are faster than just having two
      #TODO: separate for loops
      if i != minipath.len - 1: 
        po[].nodes[po.nodeIndexes[('b',v)]].illuminaSupport -= 1'u32
        if po[].nodes[po.nodeIndexes[('b',v)]].illuminaSupport == 0'u32 and
           po[].nodes[po.nodeIndexes[('b',v)]].nanoporeSupport == 0'u32:
          po[].deletedNodes.incl(v)
      po[].illuminaCounts[(u,v)] -= 1'u32
      if (u,v) in po[].nanoporeCounts:
        if po[].illuminaCounts[(u,v)] == 0'u32 and
           po[].nanoporeCounts[(u,v)] == 0'u32:
          po[].illuminaCounts.del((u,v))
          po[].nanoporeCounts.del((u,v))
          po[].weights.del((u,v))
          for j,v2 in po[].edges[u]:
            if v2 == v:
              po[].edges[u].del(j)
              break
      else:
        if po[].illuminaCounts[(u,v)] == 0'u32:
          po[].illuminaCounts.del((u,v))
          po[].weights.del((u,v))
          for j,v2 in po[].edges[u]:
            if v2 == v:
              po[].edges[u].del(j)
              break
    return true
  while true:
    if i == minipath.len:
      break
    if j == greedy_path.len:
      break
    if minipath[i] == greedy_path[j]:
      startFlag = true
      alignment.add(0'u8)
      alignIdxs.add(uint32(j))
      i += 1
      j += 1
    elif minipath[i] < greedy_path[j]:
      i += 1
      if startFlag:
        alignment.add(1'u8)
        deltaFlag = true
    else: # x[i] > y[j]
      j += 1
      if startFlag:
        alignment.add(2'u8)
        deltaFlag = true
  # if 242'u32 in minipath:
  #   echo "mini - ", minipath
  #   assert minipath[0] in po.deletedNodes
  #   # echo "new_mini - ", new_minipath
  #   echo "greedy - ", greedy_path
  #   echo alignment
  # if alignment.len == 0:
  #   return false
  if alignment.len < minipath.len:
    return false

  if deltaFlag:
    # if searchForGreedyPath(po,minipath,greedy_path,psi)
    if (alignment[0] == 0'u8 and alignment[^1] == 0'u8):
      for i in 1..<minipath.len:
        let u = minipath[i-1]
        let v = minipath[i]
        #TODO: Test if all these comparisons are faster than just having two
        #TODOL separate for loops
        if i != minipath.len - 1:
          po[].nodes[po.nodeIndexes[('b',v)]].illuminaSupport -= 1'u32
          if po[].nodes[po.nodeIndexes[('b',v)]].illuminaSupport == 0'u32 and
             po[].nodes[po.nodeIndexes[('b',v)]].nanoporeSupport == 0'u32:
            po[].deletedNodes.incl(v)
        po[].illuminaCounts[(u,v)] -= 1'u32
        if (u,v) in po[].nanoporeCounts:
          if po[].illuminaCounts[(u,v)] == 0'u32 and
             po[].nanoporeCounts[(u,v)] == 0'u32:
            po[].illuminaCounts.del((u,v))
            po[].nanoporeCounts.del((u,v))
            po[].weights.del((u,v))
            for j,v2 in po[].edges[u]:
              if v2 == v:
                po[].edges[u].del(j)
                break
        else:
          if po[].illuminaCounts[(u,v)] == 0'u32:
            po[].illuminaCounts.del((u,v))
            po[].weights.del((u,v))
            for j,v2 in po[].edges[u]:
              if v2 == v:
                po[].edges[u].del(j)
                break
      return true
    # return true
  else:
    # echo "AHH1"
    # if 17096'u32 in minipath:
    #   echo "perfect_match - ", greedy_path
    return true
  # # echo "AHH4"
  # ########################################################################
  # ###---    Didn't hit the greedy path in one or both directions    ---###
  # ###---    Try extending out by at max psi to see if we can        ---###
  # ###---    make it to greedy path. If not, return false            ---###
  # ########################################################################
  return false

proc trimMinipaths( rep_po : ptr TrimmedPOGraph,
                    greedy_path : seq[uint32],
                    psi : uint16 = 35) = 
  # echo "before branches - ", rep_po.illuminaBranches.len
  var toDelete : seq[seq[uint32]] 
  for minipath in rep_po.illuminaBranches.keys:
    if detectMinipathDeltas(rep_po,minipath,greedy_path,psi):
      toDelete.add(minipath)
  for minipath in toDelete:
    rep_po.illuminaBranches.del(minipath)
  # echo "after branches - ", rep_po.illuminaBranches.len

proc getRepresentativePaths3*(rep_po : ptr TrimmedPOGraph,
                              psi : uint16 = 10,
                              ends_delta : uint16 = 10,
                              illumina_weight:uint32 = 10) :
                                (seq[seq[uint32]],seq[uint32]) =
  ##############################################################################
  ###-------------- Collect greedy walks from each source node --------------###
  ##############################################################################
  # echo rep_po.reads.len
  let greedyWalks = walkGreedyPaths(rep_po)
  # echo greedyWalks
  # var outfile3 : File
  # discard open(outfile3,&"testing_repo2.html",fmWrite)
  # htmlOutput(rep_po[],outfile3,greedyWalks[0])
  ##############################################################################
  ###--------------- Correct reads based on these greedy walks --------------###
  ##############################################################################
  # var time1 = cpuTime()
  var differedReads : seq[int]
  var beforePaths : seq[seq[uint32]]
  for i in 0..<rep_po[].reads.len:
    var flag = true
    beforePaths.add(rep_po[].reads[i].correctedPath)
    for j in 0..<greedyWalks.len:
      # var beforePath = rep_po[].reads[i].correctedPath
      var differed = false
      # echo ""
      # echo "GREEDY: "
      # echo greedyWalks[j]
      # echo ""
      (differed,
       rep_po[].reads[i].correctedPath) = correctPaths(
                                            rep_po[].reads[i].correctedPath,
                                            greedyWalks[j],
                                            addr rep_po[].nodes,
                                            rep_po[].nodeIndexes,
                                            psi=psi,
                                            correct_ends = false,
                                            detect_ends=true)
      # updateGraph(rep_po, rep_po[].reads[i].correctedPath, beforePath)
      flag = differed and flag
    if flag:
    # if true:
      differedReads.add(i)
  for i in 0..<rep_po[].reads.len:
    # echo rep_po.reads[i].correctedPath
    # echo beforePaths[i]
    updateGraph(rep_po,rep_po[].reads[i].correctedPath,beforePaths[i])
  
  # echo "correct paths - ", cpuTime() - time1

  ##############################################################################
  ###---------- Trim mini Illumina only paths based on greedy walks ---------###
  ##############################################################################
  # time1 = cpuTime()
  for j in 0..<greedyWalks.len:
    trimMinipaths(rep_po,greedyWalks[j])
  # updateNodesAndEdges(rep_po)
  # echo "trim minipaths - ", cpuTime() - time1
  # for read in rep_po.reads:
  #   echo read.correctedPath
  
  ##############################################################################
  ###---        Update log likelihoods based on illumina weights          ---###
  ##############################################################################
  # time1 = cpuTime()
  for u in rep_po[].edges.keys():
    var total : float64
    for v in rep_po[].edges[u]:
      total += float64(rep_po[].weights[(u,v)])
    for v in rep_po[].edges[u]:
      rep_po[].logLikelihoods[(u,v)] =
        - ln(float(rep_po[].weights[(u,v)]) / total)
  # echo "update ll - ", cpuTime() - time1

  ##############################################################################
  ###---  Rewalk through the graph for each source node, constructing a   ---###
  ###---  p. queue of branches not walked.                                ---###
  ##############################################################################
  # time1 = cpuTime()
  var altPaths : HeapQueue[(int, float64, uint32, uint32)]
  # var altPaths : HeapQueue[(int, uint32, uint32)]
  var altPathSet : HashSet[(uint32,uint32)]
  # for i in differedReads:
    # let greedy_walk = greedyWalks[i]
  for j,greedy_walk in greedyWalks:
    for i in 1..<greedy_walk.len:
      let u = greedy_walk[i-1]
      let bckNodeIdx = rep_po.nodeIndexes[('b',u)]
      if rep_po[].nodes[bckNodeIdx].visited:
        break
      var firstNode = 0
      if i - 1 > int(psi):
        firstNode = int(i) - 1 - int(psi)
      rep_po[].nodes[bckNodeIdx].path = greedy_walk[firstNode..i-1]
      rep_po[].nodes[bckNodeIdx].visited = true
      let v1 = greedy_walk[i]
      altPathSet.incl((u,v1))
      for v2 in rep_po.edges[u]:
        if (u,v2) notin altPathSet:
          altPaths.push((-int(rep_po.weights[(u,v2)]),
                         -rep_po.logLikelihoods[(u,v2)],
                         v2,
                         u))
          altPathSet.incl((u,v2))
    rep_po[].nodes[rep_po.nodeIndexes[('b',greedy_walk[^1])]].visited = true

    # var outfile3 : File
    # discard open(outfile3,&"testing_repo2_{j}.html",fmWrite)
    # htmlOutput(rep_po[],outfile3,greedy_walk)
    # assert rep_po.edges[greedy_walk[^1]].len == 0
  # echo "Construct P. queue - ", cpuTime() - time1
  # echo altPaths
  # var outfile3 : File
  # discard open(outfile3,&"testing_repo2.html",fmWrite)
  # htmlOutput(rep_po[],outfile3,greedyWalks[0])
  ##############################################################################
  ###--- While the p. queue is not empty:                                 ---###
  ###--- If the node has been deleted (meaning it was not part of a       ---###
  ###---   delta > psi in len), ignore it                                 ---###
  ###--- Similarly, if the transition from the preceeding node to the     ---###
  ###---   node have been deleted, ignore                                 ---###
  ###--- If the neither node nor the transition have been deleted:        ---###
  ###---     i - build a path from -1 of the node to the next node        ---###
  ###---           previously seen in the greedy walk correct reads       ---###
  ###---           based on this short path                               ---###
  ###---    ii - store new alternative paths to the priority queue        ---###
  ###---   iii - reconstruct the graph with the corrected paths, keeping  ---###
  ###---           track of the deleted nodes                             ---###
  ##############################################################################
  # time1 = cpuTime()
  var iterations = 0
  while altPaths.len != 0:
    iterations += 1 
    var (_,_,fwdNodeIdx,node1) = altPaths.pop()
    # echo "NEW PATH TO TRAVEL:"
    # echo node1, " --> ", fwdNodeIdx
    # echo ""
    if fwdNodeIdx in rep_po.deletedNodes or
      (node1,fwdNodeIdx) notin rep_po.weights:
      continue
    # echo rep_po.weights[(node1,fwdNodeIdx)]
    var path = rep_po[].nodes[rep_po[].nodeIndexes[('b',node1)]].path
    var bckNodeIdx = rep_po.nodeIndexes[('b',fwdNodeIdx)]
    rep_po.nodes[bckNodeIdx].path = path[^min(int(psi),path.len)..^1] &
      @[fwdNodeIdx]
    path = path[^min(int(psi),path.len)..^1]  & @[fwdNodeIdx]
    # echo path
    while rep_po[].edges[fwdNodeIdx].len != 0 and
      not rep_po[].nodes[bckNodeIdx].visited:
      var weightList : seq[uint32]
      for j in rep_po[].edges[fwdNodeIdx]:
        weightList.add(uint32(rep_po[].weights[(fwdNodeIdx,j)]))
      let lastFwdNodeIdx = fwdNodeIdx
      fwdNodeIdx = rep_po.edges[fwdNodeIdx][maxIdx(weightList)]
      for i,nextNode in rep_po[].edges[lastFwdNodeIdx]:
        if nextNode == fwdNodeIdx:
          altPathSet.incl((nextNode,
                             lastFwdNodeIdx))
          continue
        if (nextNode, lastFwdNodeIdx) notin altPathSet:
          altPaths.push((-int(rep_po[].weights[(lastFwdNodeIdx,nextNode)]),
                         -rep_po[].logLikelihoods[(lastFwdNodeIdx,nextNode)],
                         nextNode,
                         lastFwdNodeIdx))
          altPathSet.incl((nextNode, lastFwdNodeIdx))
      rep_po.nodes[bckNodeIdx].visited = true
      bckNodeIdx = rep_po[].nodeIndexes[('b',fwdNodeIdx)]
      rep_po.nodes[bckNodeIdx].path = path[^min(int(psi),path.len)..^1] &
        @[fwdNodeIdx]
      path.add(fwdNodeIdx)
    if rep_po[].edges[fwdNodeIdx].len == 0:
      rep_po.nodes[bckNodeIdx].visited = true
    if rep_po[].nodes[bckNodeIdx].visited:
      for _ in 0..<int(psi):
        if rep_po[].edges[fwdNodeIdx].len == 0:
          break
        var weightList : seq[uint32]
        for j in rep_po[].edges[fwdNodeIdx]:
          weightList.add(uint32(rep_po[].weights[(fwdNodeIdx,j)]))
        fwdNodeIdx = rep_po[].edges[fwdNodeIdx][maxIdx(weightList)]
        path.add(fwdNodeIdx)
    # var reads : seq[Read]
    # echo "path - ", path
    # echo differedReads
    # var total1 = 0.0
    # var total2 = 0.0
    for i in differedReads:
      # var time3 = cpuTime()
      var j : bool
      var beforePath = rep_po.reads[i].correctedPath
      (j,
       rep_po[].reads[i].correctedPath) = 
         correctPaths(rep_po[].reads[i].correctedPath,
                      path,
                      addr rep_po[].nodes,
                      rep_po[].nodeIndexes,
                      psi,
                      false)
      # total1 += cpuTime() - time3
      # time3 = cpuTime()
      # echo rep_po.reads[i].name
      # echo rep_po.reads[i].name, " - ", beforePath
      # echo rep_po.reads[i].name, " - ", rep_po[].reads[i].correctedPath
      updateGraph(rep_po,
                  rep_po[].reads[i].correctedPath,
                  beforePath)
      # total2 += cpuTime() - time3
    # echo "Update1.1 - ", total1
    # echo "Update1.2 - ", total2
      # if iterations == 68:
      #   if 17096'u32 in rep_po[].reads[i].correctedPath and
      #      17088'u32 in rep_po[].reads[i].correctedPath:
      #     var outfile4 : File
      #     discard open(outfile4,
      #                  &"testing_repo2_{iterations}_{i}_reads.html",fmWrite)
      #     htmlOutput(rep_po[],outfile4,rep_po[].reads[i].correctedPath,path)
    # if iterations == 68:
    # # if true:
    #   var j = 0
    #   for minipath in rep_po.illuminaBranches.keys:
    #     if 17096'u32 in minipath:
    #       var outfile4 : File
    #       discard open(outfile4,
    #                    &"testing_repo2_{iterations}_{j}_minipaths.html",
    #                    fmWrite)
    #       htmlOutput(rep_po[],outfile4,minipath,path)
    #     j+=1
    # var time2 = cpuTime()
    trimMinipaths(rep_po,path)
    # echo "Update2 - ", cpuTime() - time2
    
    # time2 = cpuTime()
    # updateNodesAndEdges(rep_po)
    # echo "Update2 - ", cpuTime() - time2
  # echo "Navigate P. queue - ", cpuTime() - time1
  # var outfile4 : File
  # discard open(outfile4,"testing_repo3.html",fmWrite)
  # htmlOutput(rep_po[],outfile4)
  ##############################################################################
  ###--- Run approach similar to Stringtie approach for resolving splice  ---###
  ###--- graph; i.e. start at highest coverage node, walk heaviest path   ---###
  ###--- based on reads compatible with the walk up to that point. Remove ---###
  ###--- reads compatible with the full walk, repeat.                     ---###
  ##############################################################################
  
  # var outfile3 : File
  # discard open(outfile3,"p_queue_trimmed.html",fmWrite)
  # htmlOutput(rep_po, outfile3)

  # time1 = cpuTime()
  # let my_paths = walkHeaviestPaths( addr rep_po )
  # echo "Walk paths - ", cpuTime() - time1
  # return my_paths
  # for i,node in rep_po.nodes:
  #   # echo rep_po.nodeIndexes[('f',uint32(i))]
  #   assert node.visited or
  #          rep_po.nodeIndexes[('f',uint32(i))] in rep_po.deletedNodes
  return walkHeaviestPaths(rep_po,
                           ends_delta = int(ends_delta))


proc mismatchBasesToGraph( po : ptr TrimmedPOGraph,
                           pos : uint32,
                           base : string) : int = 
  var newNodeIndex : uint32
  var nodeIndex = po[].nodes.len
  var nodeNumber = po[].ogNodes
  let newNode = Node(nts : base,
                      alignRingPartner : int32(-1),
                      visited : false,
                      indegree : 1'u16,
                      logLikelihood : 0.0,
                      startNodeFlag : false,
                      endNodeFlag : false)

  po[].nodes.add(newNode)
  newNodeIndex = nodeNumber
  # if newNodeIndex == 1219'u32:
  #   echo "here1"
  #   echo base
  #   echo pos
  po[].nodeIndexes[('f',uint32(nodeIndex))] = uint32(nodeNumber)
  po[].nodeIndexes[('b',uint32(nodeNumber))] = uint32(nodeIndex)
  po[].edges[newNodeIndex] = @[]
  nodeIndex += 1
  nodeNumber += 1

  if po[].nodes[po[].nodeIndexes[('b', pos)]].alignRingPartner == -1'i32:
    po[].nodes[po[].nodeIndexes[('b', pos)]].alignRingPartner =
      int32(newNodeIndex)
    po[].nodes[po[].nodeIndexes[('b',newNodeIndex)]].alignRingPartner =
      int32(pos)
  else:
    var nextNode = po[].nodeIndexes[('b',
      uint32(po[].nodes[po.nodeIndexes[('b',pos)]].alignRingPartner))]
    while po[].nodes[nextNode].alignRingPartner != int32(pos):
      nextNode =
        po[].nodeIndexes[('b', uint32(po[].nodes[nextNode].alignRingPartner))]
    po[].nodes[nextNode].alignRingPartner = int32(newNodeIndex)
    po[].nodes[po[].nodeIndexes[('b',newNodeIndex)]].alignRingPartner =
      int32(pos)

  po[].ogNodes += 1'u32

  return int(newNodeIndex)

proc insertBasesToGraph( po : ptr TrimmedPOGraph,
                         pos : uint32,
                         insert : string) : seq[uint32] = 
  var newNodeIndexes : seq[uint32]
  var nodeIndex = po[].nodes.len
  var nodeNumber = po[].ogNodes
  for base in insert:
    let newNode = Node(nts : $base,
                        alignRingPartner : int32(-1),
                        visited : false,
                        indegree : 1'u16,
                        logLikelihood : 0.0,
                        startNodeFlag : false,
                        endNodeFlag : false)

    po[].nodes.add(newNode)
    newNodeIndexes.add(nodeNumber)
    # if nodeNumber == 1219'u32:
    #   echo "here2"
    #   echo insert
    po[].nodeIndexes[('f',uint32(nodeIndex))] = uint32(nodeNumber)
    po[].nodeIndexes[('b',uint32(nodeNumber))] = uint32(nodeIndex)
    po[].edges[nodeNumber] = @[]
    nodeIndex += 1
    nodeNumber += 1
  po[].ogNodes += uint32(insert.len)

  return newNodeIndexes

proc topologicalSort( po : ptr TrimmedPOGraph) =
  var nodesWithNoIncomingEdges : seq[uint32]
  var topoOrdering : seq[uint32]
  var indegrees : seq[uint16]
  for _ in po[].nodes:
    indegrees.add(0'u16)
  for i,_ in po[].nodes:
    # echo po.nodeIndexes[('f',uint32(i))]
    for fwd_neighbor in po[].edges[po.nodeIndexes[('f',uint32(i))]]:
      let bckNeighbor = po[].nodeIndexes[('b',fwd_neighbor)]
      indegrees[bckNeighbor] += 1'u16
  # echo indegrees
  # echo po.edges
  for i in 0..<po[].nodes.len:
    if indegrees[i] == 0'u16:
      nodesWithNoIncomingEdges.add(uint32(i))
  while nodesWithNoIncomingEdges.len > 0:
    let node = nodesWithNoIncomingEdges.pop()
    topoOrdering.add(node)
    # echo "node\t", node
    for fwd_neighbor in po[].edges[po[].nodeIndexes[('f',node)]]:
      # echo "fwd\t", fwd_neighbor
      let bckNeighbor = po[].nodeIndexes[('b',fwd_neighbor)]
      # echo "bck\t", bckNeighbor
      indegrees[bckNeighbor] -= 1'u16
      if indegrees[bckNeighbor] == 0'u16:
        nodesWithNoIncomingEdges.add(bckNeighbor)
  # for i,indegree in indegrees:
  #   if indegree != 0'u16:
  #     echo &"Improper Node: {i}"
  #     echo indegree
  assert topoOrdering.len == po[].nodes.len

  var newNodeIndexes : Table[(char,uint32),uint32]
  for i,j in topoOrdering:
    let bckIndex = j
    # let old_fwd_index = po[].nodeIndexes[('f',bckIndex)]
    let newFwdIndex = uint32(i)
    # if ((3000'u32 < old_fwd_index) and
    #    (old_fwd_index < 3030'u32)) or
    #    (old_fwd_index == 6456'u32):
    #   echo old_fwd_index, "\t", newFwdIndex
    newNodeIndexes[('b',newFwdIndex)] = bckIndex
    newNodeIndexes[('f',bckIndex)] = newFwdIndex
  var newEdges : Table[uint32,seq[uint32]]
  var newWeights : Table[(uint32,uint32),uint32]
  for old_u in po[].edges.keys:
    let newU = newNodeIndexes[('f',po.nodeIndexes[('b',old_u)])]
    newEdges[newU] = @[]
    for old_v in po[].edges[old_u]:
      let newV = newNodeIndexes[('f',po[].nodeIndexes[('b',old_v)])]
      newEdges[newU].add(newV)
      # assert newU < newV
      # echo old_u, " - ", old_v
      newWeights[(newU,newV)] = po[].weights[(old_u,old_v)]
  var newNanoporeCounts : Table[(uint32,uint32),uint32]
  for (old_u,old_v) in po[].nanoporeCounts.keys:
    let newU = newNodeIndexes[('f',po[].nodeIndexes[('b',old_u)])]
    let newV = newNodeIndexes[('f',po[].nodeIndexes[('b',old_v)])]
    newNanoporeCounts[(newU,newV)] = po[].nanoporeCounts[(old_u,old_v)]
  for i in 0..<po[].nodes.len:
    if po[].nodes[i].alignRingPartner != -1:
      po[].nodes[i].alignRingPartner = 
        int32(newNodeIndexes[('f',
          po[].nodeIndexes[('b',
                            uint32(po[].nodes[i].alignRingPartner))])])
  for i in 0..<po.reads.len:
    var newPath : seq[uint32]
    for old_node in po[].reads[i].correctedPath:
      let newNode = newNodeIndexes[('f',po[].nodeIndexes[('b',old_node)])]
      newPath.add(newNode)
    po[].reads[i].correctedPath = newPath
  var newIlluminaBranches : Table[seq[uint32],seq[uint32]]
  for s in po[].illuminaBranches.keys():
    var newPath : seq[uint32]
    for u in s:
      newPath.add(newNodeIndexes[('f',po[].nodeIndexes[('b',u)])])
    for i in 1..<newPath.len:
      # echo s[i-1], " - ", s[i]
      # if s[i] == 470'u32:
        # assert 1219'u32 notin po[].edges[470'u32]
      # echo newPath[i-1], " - ", newPath[i]
      assert newPath[i-1] < newPath[i]
    newIlluminaBranches[newPath] = po[].illuminaBranches[s]
  
  var newIlluminaSupported : HashSet[(uint32,uint32)]
  for (old_u,old_v) in po[].illuminaSupported.items:
    let newU = newNodeIndexes[('f',po[].nodeIndexes[('b',old_u)])]
    let newV = newNodeIndexes[('f',po[].nodeIndexes[('b',old_v)])]
    newIlluminaSupported.incl((newU,newV))

  po[].nodeIndexes = newNodeIndexes
  po[].edges = newEdges
  po[].weights = newWeights
  # po[].illumina_reads = new_illumina_reads
  po[].illuminaBranches = newIlluminaBranches
  po[].nanoporeCounts = newNanoporeCounts
  po[].illuminaSupported = newIlluminaSupported

proc illuminaPolishPOGraph*( po : ptr TrimmedPOGraph,
                             bam : Bam,
                             illumina_weight : uint32 = 10,
                             debug = true) =
  po[].nanoporeCounts = po[].weights
  var readIdToIdx : Table[string,uint32]
  var existingEdges : HashSet[(uint32,uint32)]
  for i in 0..<po.nodes.len:
    po[].nodes[i].alignRingPartner = -1'i32
  for i,read in po.reads:
    readIdToIdx[read.name] = uint32(i)
    for j in 1..<read.correctedPath.len:
      existingEdges.incl((read.correctedPath[j-1], read.correctedPath[j]))
  
  var inserts : Table[(uint32,string),seq[uint32]]
  # var illumina_reads : CountTable[seq[uint32]]
  var illuminaBranches : Table[seq[uint32],seq[uint32]]

  # var total_first_mapping_time = 0.0
  # var total_secondary_mapping_time = 0.0
  for k,read in po.reads:
    # echo read.name
    for record in bam.query(read.name):
      # let time1 = cpuTime()
      if record.flag.unmapped:
        continue
      # if record.chrom == "":
      #   continue
      # echo ""
      # echo record.cigar
      # if not (record.chrom in readIdToIdx):
      #   continue
      # if record.mapping_quality == 0:
      #   continue
      # let path = read_paths[readIdToIdx[record.chrom]]
      let path = po[].reads[readIdToIdx[record.chrom]].correctedPath
      let cigar = record.cigar
      var refIndex = record.start
      var queIndex = 0
      var traveledNodes : seq[uint32]
      var newEdges : seq[(uint32,uint32)]
      var illuminaWalks : seq[(uint16,uint16)]

      var newEdgeFlag = false
      var deleteFlag = false
      var ambiguousFlag = false
      var newBlockFlag = false
      var startFlag = true
      var newBlockStart = -1
      var newBlockOffset = 0
      var s : string
      var recordInserts : seq[(uint32,string)]
      record.sequence(s)
      # echo  s
      # echo ""
      var lastMatchOp = cigar.len
      var opNum = 0
      for op in cigar:
        if int(op.op) == 7:
          lastMatchOp = opNum
        opNum += 1
      opNum = 0
      for op in cigar:
        case int(op.op):
          of 1: #Insert
            if not startFlag:
              if opNum > lastMatchOp:
                # echo "BREAKING1"
                break
              if not newBlockFlag:
                newBlockStart = queIndex - 1 - newBlockOffset
                newBlockFlag = true
              var sequence : seq[string]
              for i in 0..<op.len:
                if record.base_at(queIndex) == 'N':
                  ambiguousFlag = true
                  break
                sequence.add($record.base_at(queIndex))
                queIndex += 1
              if ambiguousFlag:
                break
              var insertWalk : seq[uint32]
              if not ((traveledNodes[^1],sequence.join()) in inserts):
                insertWalk = insertBasesToGraph(po,
                                                traveledNodes[^1],
                                                sequence.join())
                inserts[(traveledNodes[^1],sequence.join())] = insertWalk
                recordInserts.add((traveledNodes[^1],sequence.join()))
                if not deleteFlag:
                  if not existingEdges.containsOrIncl((traveledNodes[^1],
                                                       insertWalk[0])):
                    newEdges.add((traveledNodes[^1],insertWalk[0]))
                    # echo "insert1 adding - ",
                    #      traveledNodes[^1], " ", insertWalk[0]
                for i in 1..<insertWalk.len:
                  if not existingEdges.containsOrIncl((insertWalk[i-1],
                                                       insertWalk[i])):
                    newEdges.add((insertWalk[i-1], insertWalk[i]))
                    # echo "insert1 adding - ",
                    #      insertWalk[i-1], " ", insertWalk[i]
                # newEdgeFlag = true
              else:
                insertWalk = inserts[(traveledNodes[^1],sequence.join())]
              if deleteFlag:
                if not existingEdges.containsOrIncl((traveledNodes[^1],
                                                     insertWalk[0])):
                  newEdges.add((traveledNodes[^1],insertWalk[0]))
                  # echo "insert2 adding - ",
                  #      traveledNodes[^1], " ", insertWalk[0]
                deleteFlag = false
              for node in insertWalk:
                traveledNodes.add(node)
                # illumina_walk.add(true)
              newEdgeFlag = true
            else:
              queIndex += op.len
              newBlockOffset += op.len
          of 2: #Delete
            if not startFlag:
              if opNum > lastMatchOp:
                # echo "BREAKING2"
                break
              if not newBlockFlag:
                newBlockFlag = true
                newBlockStart = queIndex - 1 - newBlockOffset
              # let start_path = traveledNodes[^1]
              refIndex += op.len
              deleteFlag = true
              newEdgeFlag = false              
              # if newEdgeFlag:
              #   newEdges.add((traveledNodes[^1],path[refIndex]))
              #   newEdgeFlag = false

              # let delete_path = @[start_path,end_path]
              # if delete_path notin deletes:
              #   newEdgeFlag = true
              #   deletes.incl(delete_path)
            else:
              refIndex += op.len

          
          # of 3: #Intron # Shouldn't be possible for us

          of 4: #Soft-clip
            # echo "Soft-clip"
            queIndex += op.len
            newBlockOffset += op.len
          
          # of 5: #Hard-clip # Shouldn't be possible for us

          # of 6: #Pad # Shouldn't be possible for us
          
          of 7: #Equal
            # if debug:
            #   echo "Equal"
            startFlag = false
            if newBlockFlag:
              illuminaWalks.add((uint16(newBlockStart),
                                 uint16(queIndex - newBlockOffset)))
              newBlockStart = -1
              newBlockFlag = false
            if deleteFlag:
              if not existingEdges.containsOrIncl((traveledNodes[^1],
                                                   path[refIndex])):
                newEdges.add((traveledNodes[^1],path[refIndex]))
                # echo "match1 adding - ",
                #      traveledNodes[^1], " ", path[refIndex]
              newEdgeFlag = false
              deleteFlag = false
            elif newEdgeFlag:
              # if (traveledNodes[^1],path[refIndex]) == ()
              if not existingEdges.containsOrIncl((traveledNodes[^1],
                                                   path[refIndex])):
                newEdges.add((traveledNodes[^1], path[refIndex]))
                # echo "match2 adding - ",
                #      traveledNodes[^1], " ", path[refIndex]
              newEdgeFlag = false
            for i in 0..<op.len:
              if debug:
                # echo "Match!"
                # echo po.nodes[po[].nodeIndexes[('b',path[refIndex])]].nts
                # echo $record.base_at(queIndex)
                # echo refIndex
                assert po[].nodes[po[].nodeIndexes[('b',
                                                    path[refIndex])]].nts == 
                                                      $record.base_at(queIndex)
              traveledNodes.add(path[refIndex])
              refIndex += 1
              queIndex += 1
          of 8: #Diff
            # if debug:
            #   echo "Diff"
            if not startFlag:
              if opNum > lastMatchOp:
                # echo cigar
                # echo opNum, " ", lastMatchOp
                # echo traveledNodes.len
                # echo traveledNodes
                # echo "BREAKING3"
                break
              if not newBlockFlag:
                newBlockStart = queIndex - 1 - newBlockOffset
                newBlockFlag = true
              for i in 0..<op.len:
                if record.base_at(queIndex) == 'N':
                  ambiguousFlag = true
                  break
                # if test.nodes[path[refIndex]].nts !=
                #   $record.base_at(queIndex):
                if debug:
                  # echo "Mismatch!"
                  # echo po.nodes[po.nodeIndexes[('b',path[refIndex])]].nts
                  # echo $record.base_at(queIndex)
                  # if path[refIndex] == 3017:
                  #   echo "here7"
                  #   echo po.nodes[po.nodeIndexes[('b',path[refIndex])]].nts
                  #   echo $record.base_at(queIndex)
                  assert po[].nodes[po[].nodeIndexes[('b',
                                                      path[refIndex])]].nts !=
                                                      $record.base_at(queIndex)
                  # echo "Diff - Mismatch! - ",
                  #      po.nodes[po.nodeIndexes[('b',path[refIndex])]].nts,
                  #      $record.base_at(queIndex)
                  # test.nodes[path[refIndex]].nts = $record.base_at(queIndex)
                  # po.nodes[po.nodeIndexes[('b',path[refIndex])]].nts =
                  #   $record.base_at(queIndex)
                var nextNode = po[].nodeIndexes[('b',path[refIndex])]
                var correctNode = -1
                for i in 0..<3:
                  if po[].nodes[nextNode].alignRingPartner == -1'i32:
                    break
                  if po[].nodes[nextNode].alignRingPartner ==
                    int32(path[refIndex]):
                    break
                  if po[].nodes[
                       po[].nodeIndexes[('b',
                         uint32(
                           po[].nodes[nextNode].alignRingPartner))]].nts ==
                              $record.base_at(queIndex):
                    correctNode = po[].nodes[nextNode].alignRingPartner
                    var flag = false
                    for v in po[].edges[traveledNodes[^1]]:
                      flag = flag or v == uint32(correctNode)
                    if not flag:
                      newEdgeFlag = true
                    break
                  else:
                    nextNode =
                      po[].nodeIndexes[('b',
                        uint32(po[].nodes[nextNode].alignRingPartner))]
                if correctNode == -1:
                  correctNode = mismatchBasesToGraph(po,
                                                     path[refIndex],
                                                     $record.base_at(queIndex))
                  # if correctNode == 11354:
                  #   echo "WHY?"
                  #   echo traveledNodes.len
                  #   echo deleteFlag
                  #   echo (traveledNodes[^1],uint32(correctNode))
                  #   echo (traveledNodes[^1],
                  #         uint32(correctNode)) in existingEdges
                  if traveledNodes.len > 0 and not deleteFlag:
                    if not existingEdges.containsOrIncl((traveledNodes[^1],
                                                         uint32(correctNode))):
                      newEdges.add((traveledNodes[^1],uint32(correctNode)))
                      # if correctNode == 11354:
                      #   echo "WHAT THE ACTUAL FUCK"
                      #   echo (traveledNodes[^1],uint32(correctNode))
                      #   echo newEdges
                      # echo "mismatch1 adding - ",
                      #      traveledNodes[^1], " ", correctNode
                  # newEdgeFlag = true
                elif newEdgeFlag and not deleteFlag:
                  if not existingEdges.containsOrIncl((traveledNodes[^1],
                                                       uint32(correctNode))):
                    newEdges.add((traveledNodes[^1],uint32(correctNode)))
                    # echo "mismatch2 adding - ",
                    #      traveledNodes[^1], " ", correctNode
                if deleteFlag:
                  if not existingEdges.containsOrIncl((traveledNodes[^1],
                                                       uint32(correctNode))):
                    newEdges.add((traveledNodes[^1],uint32(correctNode)))
                    # echo "mismatch3 adding - ",
                    #      traveledNodes[^1], " ", correctNode
                  deleteFlag = false
                traveledNodes.add(uint32(correctNode))
                refIndex += 1
                queIndex += 1
                newEdgeFlag = true
              if ambiguousFlag:
                break
            else:
              queIndex += op.len
              refIndex += op.len
              newBlockOffset += op.len
          # of 9: #Back # Shouldn't be possible for us
          else:
            echo "Something went wrong"
        opNum += 1
      if ambiguousFlag:
        # echo "here!"
        # echo newEdges
        for edge in newEdges:
          existingEdges.excl(edge)
        for insert in recordInserts:
          inserts.del(insert)
        continue
      # echo record.cigar
      # echo traveledNodes
      for (s,e) in illuminaWalks:
        # echo s, " ", e
        # echo "traveled - ",traveledNodes[s..e]
        # illuminaBranches.incl(traveledNodes[s..e])
        if traveledNodes[s..e] in illuminaBranches:
          if illuminaBranches[traveledNodes[s..e]][^1] != uint32(k):
            illuminaBranches[traveledNodes[s..e]].add(uint32(k))
        else:
          illuminaBranches[traveledNodes[s..e]] = @[uint32(k)]

        # if 21829'u32 in traveledNodes[s..e]:
        #   echo cigar
        #   echo traveledNodes
        #   echo newEdges
      # illumina_reads.inc(traveledNodes)
      for (u,v) in newEdges:
        # echo "adding - ", u, " ", v
        # assert u in po[].edges
        # if u in po[].edges:
        po[].edges[u].add(v)
      po[].nodes[po.nodeIndexes[('b',
                                 traveledNodes[0])]].illuminaSupported = true
      for i in 1..<traveledNodes.len:
        let u = traveledNodes[i-1]
        let v = traveledNodes[i]
        if (u,v) in po[].weights:
          po[].weights[(u,v)] += illumina_weight
        else:
          po[].weights[(u,v)] = illumina_weight
        po[].nodes[po.nodeIndexes[('b',v)]].illuminaSupported = true
        po[].illuminaSupported.incl((u,v))
      # po = topologicalSort(po)
      # if debug:
      #   assert queIndex == s.len
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
      #     # echo traveledNodes
      #     for i in 1..<traveledNodes.len:
      #       let node1 = traveledNodes[i-1]
      #       let node2 = traveledNodes[i]
      #       echo node1, "\t", po.nodes[po.nodeIndexes[('b',node1)]].nts
      #       echo po.edges[node1]
      #       for v in po.edges[node1]:
      #         echo node1,
      #              "\t", v, "\t",po.nodes[po.nodeIndexes[('b',v)]].nts,
      #              "\t", po.weights[(node1,v)] 
      #     echo traveledNodes[^1], "\t",
      #          po.nodes[po.nodeIndexes[('b',traveledNodes[^1])]].nts
      #   if report_flag:
      #     echo "here4"
      #     echo record.cigar
      #     # echo traveledNodes
      #     for i in 1..<traveledNodes.len:
      #       let node1 = traveledNodes[i-1]
      #       let node2 = traveledNodes[i]
      #       echo node1, "\t", po.nodes[po.nodeIndexes[('b',node1)]].nts
      #       echo po.edges[node1]
      #       for v in po.edges[node1]:
      #         echo node1, "\t", v, "\t",
      #              po.nodes[po.nodeIndexes[('b',v)]].nts, "\t",
      #              po.weights[(node1,v)] 
      #     echo traveledNodes[^1], "\t",
      #          po.nodes[po.nodeIndexes[('b',traveledNodes[^1])]].nts
      # if record.mapping_quality == 0:
      #   total_secondary_mapping_time += (cpuTime() - time1)
      # else:
      #   total_first_mapping_time += (cpuTime() - time1)

  # echo "Primary alignments - ", total_first_mapping_time
  # echo "Secondary alignments - ", total_secondary_mapping_time
  # var outfile3 : File
  # discard open(outfile3,"illumina_weighted1.html",fmWrite)
  # htmlOutput(po[], outfile3)
  # echo po.weights
  # po.illumina_reads = illumina_reads
  po.illuminaBranches = illuminaBranches
  # echo "EDGES:"
  # for u in sorted(toSeq(po.edges.keys())):
  #   echo u, " - ", po.edges[u]
  # echo "WEIGHTS:"
  # for (u,v) in sorted(toSeq(po.weights.keys())):
  #   echo u, " - ", v
  
  # for u in sorted(toSeq(po.edges.keys())):
  #   po.edges[u] = sorted(po.edges[u])
  #   for i in 1..<po.edges[u].len:
  #     assert po.edges[u][i-1] != po.edges[u][i]
  # for branch in po.illuminaBranches:
  #   for i in 1..<branch.len:
  #     # echo branch[i-1], " - ", branch[i]
  #     assert branch[i] in po.edges[branch[i-1]]
  # if 470'u32 in po.edges:
  #   echo "hmm, ", po.edges[470'u32]
  topologicalSort(po)
  # for read in po.reads:
  #   if 15813 in read.correctedPath:
  #     echo read.correctedPath
  # echo po.weights
  # assert po[].weights != po[].nanoporeCounts
  
  var illuminaCounts : Table[(uint32,uint32),uint32]
  for branch in po.illuminaBranches.keys:
    # if branch == @[21987'u32, 21990'u32, 21998'u32]:
      # echo "got here"
    for i in 1..<branch.len:
      let u = branch[i-1]
      let v = branch[i]
      if (u,v) in illuminaCounts:
        illuminaCounts[(u,v)] += 1'u32
      else:
        illuminaCounts[(u,v)] = 1'u32
      if i != branch.len - 1:
        po.nodes[po.nodeIndexes[('b',v)]].illuminaSupport += 1'u32
  po.illuminaCounts = illuminaCounts


  # var outfile4 : File
  # discard open(outfile4,"illumina_weighted2.html",fmWrite)
  # htmlOutput(po[], outfile4)

proc convertPOGraphtoTrimmedPOGraph*( po : var POGraph) : TrimmedPOGraph = 
  var logLikelihoods : Table[(uint32,uint32),float64]
  var nodeIndexes : Table[(char,uint32), uint32]
  var sourceNodes : seq[uint32]
  var endNodes : seq[uint32]
  var deletedNodes : HashSet[uint32]
  for i in 0'u32..<uint32(po.nodes.len):
    nodeIndexes[('f',i)] = i
    nodeIndexes[('b',i)] = i
    # if i == 5590'u32:
    #   echo "HMMM"
  for i in 0..<po.reads.len:
    po.reads[i].correctedPath = po.reads[i].path
  return TrimmedPOGraph( nodes : po.nodes,
                         reads : po.reads, 
                         edges : po.edges,
                         weights : po.weights,
                         ogNodes : po.ogNodes,
                         logLikelihoods : logLikelihoods,
                         nodeIndexes : nodeIndexes,
                         sourceNodes : sourceNodes,
                         endNodes : endNodes,
                         deletedNodes : deletedNodes,
                         nanoporeCounts : po.weights)

proc buildConsensusPO*( po : ptr POGraph,
                        paths : seq[seq[uint32]],
                        supports : seq[uint32],
                        read_prefix : string) : TrimmedPOGraph = 
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  var uSet,vSet : HashSet[uint32]
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
        uSet.incl(path[i])
        vSet.incl(path[i+1])
      sequence.add(po[].nodes[path[i]].nts)
    sequence.add(po[].nodes[path[^1]].nts)
    reads.add(Read(name : &"{read_prefix}_{j}",
                   length : uint32(path.len),
                   path : path,
                   support : supports[j],
                   correctedPath : path,
                   sequence : sequence.join()))
    # for i in 0..<read.correctedPath.len:
    #   echo read.correctedPath[i]
  var sourceNodes1,sourceNodes2 : HashSet[uint32]
  for i,node in po[].nodes:
    if node.startNodeFlag:
      sourceNodes1.incl(uint32(i))
  let sinkNodes = vSet - uSet
  for sink_node in sinkNodes:
    edges[sink_node] = @[]
  var nodes : seq[Node]
  var nodeIndexes = sorted(toSeq(uSet + vSet))
  var nodeIndexes2 : Table[(char,uint32),uint32]
  var logLikelihoods : Table[(uint32,uint32),float64]
  for i,j in nodeIndexes:
    nodes.add(po[].nodes[j])
    if j in sourceNodes1:
      sourceNodes2.incl(j)
    nodeIndexes2[('f',uint32(i))] = j
    nodeIndexes2[('b',j)] = uint32(i)
    if j in edges:
      var totalReadCount = 0'u32
      for k in edges[j]:
        totalReadCount += weights[(j,k)]
      let tmp = float(totalReadCount)
      for k in edges[j]:
        logLikelihoods[(j,k)] = - ln(float(weights[(j,k)]) / tmp)
  var endNodes : seq[uint32]
  for i in 0..<nodes.len:
    nodes[i].indegree = 0'u16
    if nodes[i].endNodeFlag:
      endNodes.add(nodeIndexes2[('f',uint32(i))])
    if nodes[i].alignRingPartner != -1:
      if uint32(nodes[i].alignRingPartner) notin nodeIndexes:
        nodes[i].alignRingPartner = -1
  for j in edges.keys:
    for k in edges[j]:
      nodes[nodeIndexes2[('b',k)]].indegree += 1'u16
  var trim = TrimmedPOGraph(nodes : nodes,
                            edges : edges,
                            weights : weights,
                            reads : reads,
                            ogNodes : po[].ogNodes,
                            logLikelihoods : logLikelihoods,
                            sourceNodes : toSeq(sourceNodes2),
                            endNodes : endNodes,
                            nodeIndexes : nodeIndexes2)
  topologicalSort(addr trim)
  return trim

proc getTrimmedGraphFromFastaRecord*(record : FastaRecord) : TrimmedPOGraph = 
  var nodes : seq[Node]
  var path : seq[uint32]
  var edges : Table[uint32,seq[uint32]]
  var weights : Table[(uint32,uint32),uint32]
  var nodeIndexes : Table[(char,uint32), uint32]
  for i,base in record.sequence:
    nodes.add(Node( nts : $base,
                    # supportingReads : @[0'u32],
                    alignRingPartner : -1'i32,
                    nanoporeSupport : 1'u32))
    let v = uint32(i)
    path.add(v)
    nodeIndexes[('f',v)] = v
    nodeIndexes[('b',v)] = v
    if i != 0:
      let u = uint32(i - 1)
      edges[u] = @[v]
      weights[(u,v)] = 1'u32
  edges[uint32(record.sequence.len - 1)] = @[]
  # assert record.sequence.len == nodes.len
  # echo edges
  let reads = @[Read(name : record.readId, 
                     length : uint32(record.sequence.len),
                     path : path,
                     correctedPath : path,
                     sequence : record.sequence)]
  return TrimmedPOGraph(nodes : nodes,
                        reads : reads,
                        edges : edges,
                        weights : weights,
                        ogNodes : uint32(nodes.len),
                        nodeIndexes : nodeIndexes,
                        nanoporeCounts : weights)

proc getSequenceFromPath*(po : TrimmedPOGraph, path : seq[uint32]) : string = 
  var sequence1 : seq[string]
  for idx in path:
    sequence1.add(po.nodes[po.nodeIndexes[('b',idx)]].nts)
  return sequence1.join()

proc getFastaRecordsFromTrimmedPOGraph*(po : ptr TrimmedPOGraph,
                                        paths : seq[seq[uint32]],
                                        supports : seq[uint32],
                                        label : string) : seq[FastaRecord] = 
  var records : seq[FastaRecord]
  for i,path in paths:
    var sequence : seq[string]
    for fwd_u in path:
      let bckU = po[].nodeIndexes[('b',fwd_u)]
      sequence.add(po[].nodes[bckU].nts)
    records.add(FastaRecord(readId : &"{label}_{i}_{supports[i]}",
                            sequence : sequence.join()))
  return records
