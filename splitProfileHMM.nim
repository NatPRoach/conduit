proc initSplitPHMM( po : TrimmedPOGraph,representative_paths : seq[seq[uint32]]) : SplitProfileHMM = 
  # let nt_to_idx = {'A': 0'u8, 'a' : 0'u8,
  #                  'G': 1'u8, 'g' : 1'u8,
  #                  'C': 2'u8, 'c' : 2'u8,
  #                  'T': 3'u8, 't' : 3'u8, 'U': 3'u8, 'u' : 3'u8}.toTable()
  # let idx_to_nt = { 0'u8 : 'A',
  #                   1'u8 : 'G',
  #                   2'u8 : 'C',
  #                   3'u8 : 'T'}.toTable()
  var new_edges : Table[uint32,seq[uint32]]
  var u_set,v_set,source_nodes,sink_nodes : HashSet[uint32]
  var source_nodes2,sink_nodes2 : seq[uint32]
  # echo "source_nodes1: ", po.source_nodes
  # echo "sink_nodes1: ", po.end_nodes
  # echo po.end_nodes
  for i in po.source_nodes:
    source_nodes.incl(i)
  for i in po.end_nodes:
    sink_nodes.incl(i)
  for u in po.edges.keys:
    let u_back = po.node_indexes[('b',u)]
    if u in source_nodes:
      source_nodes2.add(u_back)
    u_set.incl(u_back)
    if po.edges[u].len > 0:
      for v in po.edges[u]:
        let v_back = po.node_indexes[('b',v)]
        if v in sink_nodes:
          sink_nodes2.add(v_back)
        v_set.incl(v_back)
        if u_back in new_edges:
          new_edges[u_back].add(v_back)
        else:
          new_edges[u_back] = @[v_back]
    else:
      new_edges[u_back] = @[]
  # echo "source_nodes2: ", source_nodes2
  # echo "sink_nodes2: ", sink_nodes2
  var end_nodes : HashSet[uint32]
  var node_lengths : seq[uint32]
  var node_starts = @[0'u32]
  for node in po.nodes:
    node_lengths.add(uint32(node.nts.len))
    node_starts.add(node_starts[^1] + uint32(node.nts.len))
  node_starts.del(node_starts.len - 1)
  let num_match_states = sum(node_lengths)
  let num_delete_states = num_match_states
  let num_insert_states = num_match_states + uint32(source_nodes2.len)
  let num_phmm_nodes = 3'u32*num_match_states + uint32(source_nodes2.len) + 2'u32
  let start_node_idx = num_phmm_nodes - 2
  let end_node_idx = num_phmm_nodes - 1
  var transitions = newSeqWith(int(num_phmm_nodes),newSeq[float64](int(num_phmm_nodes)))
  var emissions = newSeqWith(int(num_phmm_nodes),@[1.0'f64,1.0'f64,1.0'f64,1.0'f64]) # A G C T emission counts for each node
  ########################################################
  ### INITIALIZE START NODE TO SOURCE NODE TRANSITIONS ###
  ########################################################
  var source_node_order : Table[uint32,uint32]
  for i,j in source_nodes2:
    source_node_order[j] = uint32(i)
    let i1 = 3'u32*num_match_states + uint32(i)
    let m2 = node_starts[j]
    let d2 = m2 + num_match_states
    transitions[start_node_idx][m2] += 1.0 #S -> M2
    transitions[start_node_idx][d2] += 1.0 #S -> D2
    transitions[start_node_idx][i1] += 1.0 #S -> I1
    transitions[i1][i1] += 1.0             #I1 -> I1
    transitions[i1][m2] += 1.0             #I1 -> M2
    transitions[i1][d2] += 1.0             #I1 -> D2
  
  for i, node_start in node_starts:
    ###############################################
    ### INITIALIZE LINEAR BLOCKS OF TRANSITIONS ###
    ###############################################
    for j in node_start..<node_start + node_lengths[i] - 1:
      let m1 = j
      let m2 = m1 + 1'u32
      let d1 = m1 + num_match_states
      let d2 = d1 + 1'u32
      let i1 = d1 + num_match_states
      transitions[m1][m2] += 1.0'f64 #M1 -> M2
      transitions[m1][d2] += 1.0'f64 #M1 -> D2
      transitions[m1][i1] += 1.0'f64 #M1 -> I1
      transitions[d1][m2] += 1.0'f64 #D1 -> M2
      transitions[d1][d2] += 1.0'f64 #D1 -> D2
      transitions[d1][i1] += 1.0'f64 #D1 -> I2
      transitions[i1][m2] += 1.0'f64 #I1 -> M2
      transitions[i1][d2] += 1.0'f64 #I1 -> D2
      transitions[i1][i1] += 1.0'f64 #I1 -> I1
    let m1 = node_start + node_lengths[i] - 1'u32
    let d1 = m1 + num_match_states
    let i1 = d1 + num_match_states
    transitions[m1][i1] += 1.0'f64 #M1 -> I1
    transitions[d1][i1] += 1.0'f64 #D1 -> I1
    transitions[i1][i1] += 1.0'f64 #I1 -> I1
    ###############################################
    ###   INITIALIZE NODE TO NODE TRANSITIONS   ###
    ###############################################
    for j in new_edges[uint32(i)]:
      let m2 = node_starts[j]
      let d2 = m2 + num_match_states
      transitions[m1][m2] += 1.0'f64 #M1 -> M2
      transitions[m1][d2] += 1.0'f64 #M1 -> D2
      transitions[d1][m2] += 1.0'f64 #D1 -> M2
      transitions[d1][d2] += 1.0'f64 #D1 -> D2
      transitions[i1][m2] += 1.0'f64 #I1 -> M2
      transitions[i1][d2] += 1.0'f64 #I1 -> D2
    ###############################################
    ### INITIALIZE NODE TO END NODE TRANSITIONS ###
    ###############################################
    if uint32(i) in sink_nodes2:
      transitions[m1][end_node_idx] += 1.0'f64
      transitions[d1][end_node_idx] += 1.0'f64
      transitions[i1][end_node_idx] += 1.0'f64
      end_nodes.incl(m1)
      end_nodes.incl(d1)
      end_nodes.incl(i1) 
  ###############################################
  ###     FIND TRANSITIONS INTO EACH NODE     ###
  ###############################################
  # TODO: can likely do this in the loop above (and should in optimization steps)
  # wait("allocate in_nodes")
  var in_nodes : Table[uint32,seq[uint32]]
  for i in 0..<num_phmm_nodes:
    for j in 0..<num_phmm_nodes:
      if transitions[i][j] != 0.0'f64:
        if j in in_nodes:
          in_nodes[j].add(i)
        else:
          in_nodes[j] = @[i]
  # wait("check in_nodes")
  # wait("allocate in_node_array")
  # var in_node_array : seq[uint32]
  let in_node_array_base_len = 3'u32 * num_phmm_nodes + 100'u32
  var in_node_array = newSeq[uint32](in_node_array_base_len)
  # var in_node_array_bounds : seq[uint32]
  var in_node_array_bounds = newSeq[uint32](num_phmm_nodes)
  var counter = 0'u32
  var counter2 = 0'u32
  var counter3 = 0'u32
  for i in 0..<num_phmm_nodes:
    var row_idxs : seq[uint32]
    var counter4 = 0'u32
    for j in 0..<num_phmm_nodes:
      if transitions[j][i] != 0.0'f64:
        # row_idxs.add(j)
        if counter2 < in_node_array_base_len:
          in_node_array[counter2] = j
        else:
          in_node_array.add(j)
        counter2 += 1'u32
        counter4 += 1'u32
    # in_node_array = in_node_array & row_idxs
    in_node_array_bounds[counter3]= counter
    counter += counter4
    counter3 += 1'u32
  # wait("check in_node_array")
  let in_node_array2 = vector(in_node_array)
  let in_node_array_bounds2 = vector(in_node_array_bounds)
  
  
  ###################################################################
  ###  ADD SELF TRANSITION IN END NODE TO AVOID DIVIDING BY ZERO  ###
  ###################################################################
  transitions[end_node_idx][end_node_idx] += 1.0'f64
  ###################################################################
  ### FIND BEST PATHS THROUGH THE PHMM FOR EACH READ, UPDATE PHMM ###
  ###################################################################
  var phmm_paths : seq[seq[uint32]]
  for i,path in representative_paths:
    let (path,nts) = reducePath(path, po.nodes, po.node_indexes)
    var start_time = cpuTime()
    let alignment = alignSeqs(po.reads[i].sequence,nts)
    var end_time = cpuTime()
    echo "Alignment: ---", end_time - start_time, "--- seconds"
    start_time = cpuTime()
    var seq_idx,align_idx = 0
    var phmm_path = @[start_node_idx]
    var last_insert_node = 3'u32 * num_match_states + source_node_order[path[0]]
    for node_idx in path:
      var node_char_count = 0'u32
      let node = po.nodes[node_idx]
      while node_char_count < uint32(node.nts.len):
        let m_idx = node_starts[node_idx] + node_char_count
        let d_idx = m_idx + num_match_states
        let i_idx = d_idx + num_match_states
        if alignment[align_idx] == 0'u8: # Match
          phmm_path.add(m_idx)
          emissions[m_idx][NT_TO_IDX[int(po.reads[i].sequence[seq_idx])]] += 1.0'f64
          last_insert_node = i_idx
          seq_idx += 1
          node_char_count += 1
        elif alignment[align_idx] == 1'u8: # Insertion
          phmm_path.add(last_insert_node)
          emissions[last_insert_node][NT_TO_IDX[int(po.reads[i].sequence[seq_idx])]] += 1.0'f64
          seq_idx += 1
        else: # alignment[align_idx] == 2; Deletion
          phmm_path.add(d_idx)
          last_insert_node = i_idx
          node_char_count += 1
        align_idx += 1
    phmm_path.add(end_node_idx)
    for i in 0..<phmm_path.len-1:
      assert transitions[phmm_path[i]][phmm_path[i+1]] != 0.0'f64
      transitions[phmm_path[i]][phmm_path[i+1]] += 1.0'f64
    phmm_paths.add(phmm_path)
    end_time = cpuTime()
    echo "Update:    ---", end_time - start_time, "--- seconds"
  # var non_zeros : seq[float64]
  for i in 0..<num_phmm_nodes:
    var row_sum = sum(transitions[i])
    for j in 0..<num_phmm_nodes:
      # if transitions[i][j] != 0.0'f64:
      #   non_zeros.add(transitions[i][j] / row_sum)
      transitions[i][j] = ln(transitions[i][j] / row_sum)
    row_sum = sum(emissions[i])
    for j in 0..3:
      emissions[i][j] = ln(emissions[i][j] / row_sum / 0.25'f64)
  # for i,read in po.reads:
  #   echo ">", read.name
  #   echo correctRead(addr emissions,read.sequence,)
  return SplitProfileHMM( nodes : po.nodes,
                          num_source_nodes : uint32(source_nodes2.len),
                          source_node_order : source_node_order,
                          num_match_states : num_match_states,
                          num_delete_states : num_delete_states,
                          num_insert_states : num_insert_states,
                          num_phmm_nodes : num_phmm_nodes,
                          sink_nodes : sink_nodes2,
                          start_node : start_node_idx,
                          end_node : end_node_idx,
                          end_nodes : end_nodes,
                          in_nodes : in_nodes,
                          in_node_array : in_node_array2,
                          in_node_array_bounds : in_node_array_bounds2,
                          transitions : matrix(transitions),
                          emissions : matrix(emissions))

proc viterbiScore( phmm : ptr SplitProfileHMM, read : string, max_threshold : float = Inf) : float64 = 
  # Initialize the dynamic programming matrix
  var time1 = cpuTime()
  var dp_matrix = newSeqWith(read.len + 1,newSeq[float64](int(phmm[].num_phmm_nodes)))
  var time2 = cpuTime()
  echo "Initialization: ", time2 - time1
  #Transition from start node to initial insert states
  for j in 3'u32 * phmm[].num_match_states..<3'u32 * phmm[].num_match_states + phmm[].num_source_nodes:
    for i in 1..read.len: # range(1,read.len + 1)
      dp_matrix[i][j] = phmm[].emissions[int(j),int(NT_TO_IDX[int(read[i-1])])] + dp_matrix[i-1][j] + phmm[].transitions[int(j),int(j)]
  for j in 0..<int(phmm[].num_match_states):
    let match_idx = j
    let delete_idx = j + int(phmm[].num_match_states)
    let insert_idx = j + 2 * int(phmm[].num_match_states)
    for i in 1..read.len:
      # Match state
      var max_ll = -Inf
      # for k in phmm[].in_nodes[uint32(match_idx)]:
      # for k in phmm[].in_node_array[phmm[].in_node_array_bounds[match_idx]..<phmm[].in_node_array_bounds[match_idx+1]]:
      for l in int(phmm[].in_node_array_bounds[match_idx])..<int(phmm[].in_node_array_bounds[match_idx+1]):
        let k = phmm[].in_node_array[l]
        let ll = dp_matrix[i-1][k] + phmm[].transitions[int(k),match_idx]
        if ll > max_ll:
          max_ll = ll
      dp_matrix[i][match_idx] = phmm[].emissions[match_idx,int(NT_TO_IDX[int(read[i-1])])] + max_ll
      if dp_matrix[i][match_idx] > max_threshold:
        echo dp_matrix[i][match_idx]
        return dp_matrix[i][match_idx]
      # Delete state
      max_ll = -Inf
      # for k in phmm[].in_nodes[uint32(delete_idx)]:
      # for k in phmm[].in_node_array[phmm[].in_node_array_bounds[delete_idx]..<phmm[].in_node_array_bounds[delete_idx+1]]:
      for l in int(phmm[].in_node_array_bounds[delete_idx])..<int(phmm[].in_node_array_bounds[delete_idx+1]):
        let k = phmm[].in_node_array[l]
        let ll = dp_matrix[i][k] + phmm[].transitions[int(k),delete_idx]
        if ll > max_ll:
          max_ll = ll
      dp_matrix[i][delete_idx] = max_ll
      # Insert state
      max_ll = -Inf
      # for k in phmm[].in_nodes[uint32(insert_idx)]:
      # for k in phmm[].in_node_array[phmm[].in_node_array_bounds[insert_idx]..<phmm[].in_node_array_bounds[insert_idx+1]]:
      for l in int(phmm[].in_node_array_bounds[insert_idx])..<int(phmm[].in_node_array_bounds[insert_idx+1]):
        let k = phmm[].in_node_array[l]
        let ll = dp_matrix[i-1][k] + phmm[].transitions[int(k),insert_idx]
        if ll > max_ll:
          max_ll = ll
      dp_matrix[i][insert_idx] = phmm[].emissions[insert_idx,int(NT_TO_IDX[int(read[i-1])])] + max_ll
      if dp_matrix[i][insert_idx] > max_threshold:
        echo dp_matrix[i][insert_idx]
        return dp_matrix[i][insert_idx]
  let i = read.len
  var max_ll = -Inf
  # for k in phmm[].in_nodes[phmm[].end_node]:
  # for k in phmm[].in_node_array[phmm[].in_node_array_bounds[phmm[].end_node]..<phmm[].in_node_array.len]:
  for l in int(phmm[].in_node_array_bounds[int(phmm[].end_node)])..<phmm[].in_node_array.len:
    let k = phmm[].in_node_array[l]
    let ll = dp_matrix[i][k] + phmm[].transitions[int(k),int(phmm[].end_node)]
    if ll > max_ll:
      max_ll = ll
  echo max_ll
  return max_ll

proc getWindowIterator( phmm : ptr SplitProfileHMM,j,read_len : uint32, window : uint32 = 10) : seq[uint32] =
  var s : HashSet[uint32]
  for w in phmm[].windows[j]:
    for i in w - window..w + window:
      if 1'u32 <= i and i <= read_len:
        s.incl(i)
  # let s_p = toSeq(s)
  # echo s_p
  # return s_p
  return sorted(toSeq(s))

proc findWindows( phmm : ptr SplitProfileHMM, node_idx : var uint32, read_idx : var uint32) =
  var multipath : seq[uint32]
  while node_idx in phmm[].in_nodes and read_idx != 0 and multipath.len == 0:
    if node_idx in phmm[].windows:
      phmm[].windows[node_idx].incl(read_idx)
    else:
      phmm[].windows[node_idx] = toHashSet([read_idx])
    read_idx -= 1'u32 
    multipath = @[]
    for idx in phmm[].in_nodes[node_idx]:
      if idx < phmm[].num_match_states:
        multipath.add(idx)
    if multipath.len == 1:
      node_idx = multipath[0]
      multipath = @[]
    elif multipath.len == 0:
      break
  if (multipath.len > 0) and (read_idx != 0) and (node_idx in phmm[].in_nodes):
    for path in multipath:
      var path = path
      var read_idx2 = read_idx
      findWindows(phmm,path,read_idx2)

proc windowedViterbiScore( phmm : ptr SplitProfileHMM, read : string, max_threshold : float = Inf,window : uint32 = 10): float64 = 
  phmm[].windows = initTable[uint32,HashSet[uint32]]()
  var start_time = cpuTime()
  for end_node in phmm[].end_nodes:
    if end_node < phmm[].num_match_states:
      var end_node = end_node
      var read_idx = uint32(read.len)
      findWindows(phmm,end_node,read_idx)
  var end_time = cpuTime()
  echo "Finding Windows: ---", end_time - start_time, "--- seconds"
  # echo phmm[].windows
  # for key in sorted(toSeq(phmm[].windows.keys)):
  #   echo key, sorted(toSeq(phmm[].windows[key]))

  # Initialize the dynamic programming matrix
  var dp_matrix = newSeqWith(read.len + 1,newSeq[float64](int(phmm[].num_phmm_nodes)))
  #Transition from start node to initial insert states
  for j in 3'u32 * phmm[].num_match_states..<3'u32 * phmm[].num_match_states + phmm[].num_source_nodes:
    for i in 1..read.len: # range(1,read.len + 1)
      dp_matrix[i][j] = phmm[].emissions[int(j),int(NT_TO_IDX[int(read[i-1])])] + dp_matrix[i-1][j] + phmm[].transitions[int(j),int(j)]
  for j in 0..<int(phmm[].num_match_states):
    let match_idx = j
    let delete_idx = j + int(phmm[].num_match_states)
    let insert_idx = j + 2 * int(phmm[].num_match_states)
    # for i in 1..read.len:
    for i in getWindowIterator(phmm,uint32(j),uint32(read.len),window):
      # Match state
      var max_ll = -Inf
      for k in phmm[].in_nodes[uint32(match_idx)]:
        let ll = dp_matrix[i-1][k] + phmm[].transitions[int(k),match_idx]
        if ll > max_ll:
          max_ll = ll
      dp_matrix[i][match_idx] = phmm[].emissions[match_idx,int(NT_TO_IDX[int(read[i-1])])] + max_ll
      if dp_matrix[i][match_idx] > max_threshold:
        echo dp_matrix[i][match_idx]
        return dp_matrix[i][match_idx]
      # Delete state
      max_ll = -Inf
      for k in phmm[].in_nodes[uint32(delete_idx)]:
        let ll = dp_matrix[i][k] + phmm[].transitions[int(k),delete_idx]
        if ll > max_ll:
          max_ll = ll
      dp_matrix[i][delete_idx] = max_ll
      # Insert state
      max_ll = -Inf
      for k in phmm[].in_nodes[uint32(insert_idx)]:
        let ll = dp_matrix[i-1][k] + phmm[].transitions[int(k),insert_idx]
        if ll > max_ll:
          max_ll = ll
      dp_matrix[i][insert_idx] = phmm[].emissions[insert_idx,int(NT_TO_IDX[int(read[i-1])])] + max_ll
      if dp_matrix[i][insert_idx] > max_threshold:
        echo dp_matrix[i][insert_idx]
        return dp_matrix[i][insert_idx]
  let i = read.len
  var max_ll = -Inf
  for k in phmm[].in_nodes[phmm[].end_node]:
    let ll = dp_matrix[i][k] + phmm[].transitions[int(k),int(phmm[].end_node)]
    if ll > max_ll:
      max_ll = ll
  echo max_ll
  return max_ll

proc alignSeqs( x,y:string,a:int32 = 2,b:int32 = 4,o:int32=4,e:int32=2) : seq[uint8] = 

  # let time1 = cpuTime()
  var
    score_matrix = newSeqWith(x.len + 1, newSeq[int32](y.len + 1))
    trace_matrix = newSeqWith(x.len + 1, newSeq[uint8](y.len + 1))
  for i in 1'i32..int32(x.len):
    score_matrix[i][0] = -o - i*e
  for j in 1'i32..int32(y.len):
    score_matrix[0][j] = -o - j*e
  
  var col_score = @[0'i32]
  for j in 1'i32..int32(y.len + 1):
    col_score.add(-2'i32 * o - j*e)
  
  # let time2 = cpuTime()
  for i in 1'i32..int32(x.len):
    var row_score = -2'i32 * o - i*e
    for j in 1'i32..int32(y.len):
      var no_gap_score : int32
      if x[i-1'i32] == y[j-1'i32]:
        no_gap_score = score_matrix[i-1'i32][j-1'i32] + a
      else:
        no_gap_score = score_matrix[i-1'i32][j-1'i32] - b
      row_score = max(score_matrix[i][j-1'i32] - o, row_score) - e
      col_score[j] = max(score_matrix[i-1'i32][j] - o,col_score[j]) - e
      let best_score = max(max(no_gap_score,col_score[j]),row_score)
      score_matrix[i][j] = best_score
      
      var possible_origins = 0'u8
      if best_score == no_gap_score:
        possible_origins += 1
      if best_score == col_score[j]:
        possible_origins += 2
      if best_score == row_score:
        possible_origins += 4
      trace_matrix[i][j] = possible_origins
  # let time3 = cpuTime()
  var i = x.len
  var j = y.len
  var walkback : seq[uint8]
  while i != 0 and j != 0:
    if trace_matrix[i][j] == 1'u8:
      i -= 1
      j -= 1
      walkback.add(0'u8)
    elif trace_matrix[i][j] == 2'u8:
      i -= 1
      walkback.add(1'u8)
    elif trace_matrix[i][j] == 4'u8:
      j -= 1
      walkback.add(2'u8)
    elif trace_matrix[i][j] == 3'u8: #Ambiguous alignment
      i -= 1 #Arbitrarily decide to use match
      j -= 1
      walkback.add(0'u8)
    elif trace_matrix[i][j] == 5'u8: #Ambiguous alignment
      i -= 1 #Arbitrarily decide to use match
      j -= 1
      walkback.add(0'u8)
    elif trace_matrix[i][j] == 6'u8: #Ambiguous alignment
      j -= 1 #Arbitrarily decide to use deletion
      walkback.add(2'u8)
    elif trace_matrix[i][j] == 7'u8: #Ambiguous alignment
      i -= 1 #Arbitrarily decide to use match
      j -= 1
      walkback.add(0'u8)
  for _ in 0..<i:
    walkback.add(1'u8)
  for _ in 0..<j:
    walkback.add(2'u8)
  walkback = reversed(walkback)
  # let time4 = cpuTime()
  # echo "Init:     ", time2 - time1
  # echo "Fill-in:  ", time3 - time2
  # echo "Walkback: ", time4 - time3
  # echo "Total:    ", time4 - time1
  return walkback


proc getSplitPHMM(po : var POGraph, psi:uint16 = 10) : SplitProfileHMM = 
  let paths = getRepresentativePaths(po, psi = psi)
  let closest = assignReadsToPaths(po.reads,paths)
  var selected_paths : seq[seq[uint32]]
  for c in closest:
    selected_paths.add(paths[c])
    # for idx in paths[c]:
    #   echo idx
  
  var rep_po = constructGraphFromPaths(addr po, selected_paths)
  # for idx in sorted(toSeq(rep_po.node_indexes.keys)):
  #   echo "idx, ", idx, " ", rep_po.node_indexes[idx]
  for read in rep_po.reads:
    rep_po.nodes[rep_po.node_indexes[('b',read.corrected_path[0])]].start_node_flag = true
    rep_po.nodes[rep_po.node_indexes[('b',read.corrected_path[^1])]].end_node_flag = true
    
  rep_po = collapseLinearStretches(rep_po)
  # for idx in sorted(toSeq(rep_po.node_indexes.keys)):
  #   echo "idx, ", idx, " ", rep_po.node_indexes[idx]
  return initSplitPHMM(rep_po,selected_paths)

var infile : File
discard open(infile,paramStr(1))
var time1 = cpuTime()
var test =  initPOGraph(infile)
var time2 = cpuTime()
echo "POGraph Initial Assembly:    ", time2 - time1

time1 = cpuTime()
var (trim,selected_paths) = trimAndCollapsePOGraph(test,psi = 15)
time2 = cpuTime()
echo "POGraph Trimming & Collapsing:    ", time2 - time1

var outfile1 : File
discard open(outfile1,paramStr(2),fmWrite)
writeCorrectedReads(trim,outfile1)

time1 = cpuTime()
var bam:Bam
discard open(bam,"test_illumina_to_nanopore.bam",index=true)
var trim2 = illuminaPolishPOGraph(trim,bam,update_weights=true,debug=false)
time2 = cpuTime()
trim2 = getRepresentativePaths2(trim2,psi = 15)
echo "POGraph Illumina Polishing:    ", time2 - time1
var outfile2 : File
discard open(outfile2,paramStr(3),fmWrite)
writeCorrectedReads(trim2,outfile2)



time1 = cpuTime()
# var phmm = getSplitPHMM(test)
var phmm = initSplitPHMM(trim, selected_paths)
time2 = cpuTime()
echo "Split pHMM Assembly: ", time2 - time1

echo type(phmm)
echo sizeof(phmm)
echo sizeof(phmm.in_nodes)

var size = 0
for key in phmm.in_nodes.keys:
  size += sizeof(phmm.in_nodes[key])

for i,read in test.reads:
  time1 = cpuTime()
  discard viterbiScore(addr phmm, read.sequence)
  time2 = cpuTime()
  echo "Viterbi walk, read ", i, " ", time2 - time1, " seconds"
  time1 = cpuTime()
  discard windowedViterbiScore(addr phmm,read.sequence,100)
  time2 = cpuTime()
  echo "Windowed Viterbi walk, read ", i, " ", time2 - time1, " seconds"
