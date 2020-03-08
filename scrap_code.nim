
type
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


  # for i in 0..<greedy_walks.len:
  #   var to_remove : seq[seq[uint32]]
  #   for minipath in rep_po[].illumina_branches:
  #     var print_flag = false
  #     # if 17278 in minipath:
  #     #   echo "greedy!!! - ", greedy_walks[i]
  #     #   echo ""
  #     #   echo "mini!!! - ", minipath
  #     #   echo ""
  #     #   print_flag = true
  #     ###TODO: reduce the number of copying of minipath that needs to happen, redundant and slows us down.
  #     var minipath2 = rep_po.nodes[rep_po.node_indexes[('b',minipath[0])]].path[^min(rep_po.nodes[rep_po.node_indexes[('b',minipath[0])]].path.len,int(psi))..^1] & minipath[1..^1]
  #     var fwd_node_idx2 = minipath[^1]
  #     for _ in 0'u16..<psi:
  #       if rep_po[].edges[fwd_node_idx2].len > 1:
  #         var weight_list : seq[uint32]
  #         for j in rep_po[].edges[fwd_node_idx2]:
  #           weight_list.add(uint32(rep_po[].weights[(fwd_node_idx2,uint32(j))]))
  #         fwd_node_idx2 = rep_po[].edges[fwd_node_idx2][maxIdx(weight_list)]
  #       elif rep_po[].edges[fwd_node_idx2].len == 1:
  #         fwd_node_idx2 = rep_po[].edges[fwd_node_idx2][0]
  #       else:
  #         break
  #       minipath2.add(fwd_node_idx2)
  #     var minipath3 = minipath2
  #     var (_,after_path) = correctPaths(minipath2, greedy_walks[i],addr rep_po[].nodes, rep_po[].node_indexes,psi = psi,correct_ends = false)
  #     let align = alignPaths(minipath3, after_path)
  #     var delete_flag = false
  #     var idx = 0
  #     var diff_flag = false
  #     for op in align:
  #       case op:
  #         of 0'u8:
  #           if diff_flag:
  #             if (minipath3[idx-1],minipath3[idx]) == (17197'u32,17204'u32):
  #               echo "THIS IS BAD1"
  #             if print_flag:
  #               echo "deleting1 - ", minipath[idx-1], " -> ", minipath[idx]
  #             rep_po[].deleted_nodes.incl(minipath3[idx])
  #             rep_po[].weights.del((minipath3[idx-1],minipath3[idx]))
  #             for l,k in rep_po[].edges[minipath3[idx-1]]:
  #               if k == minipath3[idx]:
  #                 rep_po[].edges[minipath3[idx-1]].del(l)
  #                 break
  #           diff_flag = false
  #           idx += 1
  #         of 1'u8: #ins in OG path; del in after
  #           diff_flag = true
  #           rep_po[].deleted_nodes.incl(minipath3[idx])
  #           rep_po[].weights.del((minipath3[idx-1],minipath3[idx]))
  #           if (minipath3[idx-1],minipath3[idx]) == (17197'u32,17204'u32):
  #             echo "THIS IS BAD2"
  #           for l,k in rep_po[].edges[minipath3[idx-1]]:
  #             if k == minipath3[idx]:
  #               rep_po[].edges[minipath3[idx-1]].del(l)
  #               if print_flag:
  #                 echo "deleting2 - ", minipath[idx-1], " -> ", minipath[idx]
  #               break
  #           idx += 1
  #           delete_flag = true
  #         else: #2'u8 del in OG path; ins in after
  #           diff_flag = true
  #           delete_flag = true
  #     if delete_flag:
  #       to_remove.add(minipath)
  #   for s in to_remove:
  #     rep_po[].illumina_branches.excl(s)


    # var to_remove : seq[seq[uint32]]
    # for minipath in rep_po.illumina_branches:
    #   ###TODO: reduce the number of copying of minipath that needs to happen, redundant and slows us down.
    #   # var minipath2 = minipath
    #   var minipath2 = rep_po.nodes[rep_po.node_indexes[('b',minipath[0])]].path[^min(rep_po.nodes[rep_po.node_indexes[('b',minipath[0])]].path.len,int(psi))..^1] & minipath[1..^1]
    #   var fwd_node_idx2 = minipath[^1]
    #   for _ in 0'u16..<psi:
    #     if rep_po[].edges[fwd_node_idx2].len > 1:
    #       var weight_list : seq[uint32]
    #       for j in rep_po[].edges[fwd_node_idx2]:
    #         weight_list.add(uint32(rep_po[].weights[(fwd_node_idx2,uint32(j))]))
    #       fwd_node_idx2 = rep_po[].edges[fwd_node_idx2][maxIdx(weight_list)]
    #     elif rep_po[].edges[fwd_node_idx2].len == 1:
    #       fwd_node_idx2 = rep_po[].edges[fwd_node_idx2][0]
    #     else:
    #       break
    #     minipath2.add(fwd_node_idx2)
    #   var minipath3 = minipath2
    #   var (_,after_path) = correctPaths(minipath2, path,addr rep_po[].nodes, rep_po[].node_indexes,psi = psi,correct_ends = false)
    #   let align = alignPaths(minipath3, after_path)
    #   var delete_flag = false
    #   var idx = 0
    #   var diff_flag = false
    #   var print_flag = false
    #   # if 17278 in minipath:
    #   #   echo "path - ", path
    #   #   echo ""
    #   #   echo "mini - ", minipath
    #   #   echo ""
    #   #   print_flag = true

    #   for op in align:
    #     case op:
    #       of 0'u8:
    #         if diff_flag:
    #           rep_po[].deleted_nodes.incl(minipath3[idx])
    #           if (minipath3[idx-1],minipath3[idx]) == (17197'u32,17204'u32):
    #             echo "THIS IS BAD3"
    #           rep_po[].weights.del((minipath3[idx-1],minipath3[idx]))
    #           if print_flag:
    #             echo "deleting1 - ", minipath3[idx-1], " -> ", minipath3[idx]
    #           for l,k in rep_po[].edges[minipath3[idx-1]]:
    #             if k == minipath3[idx]:
    #               rep_po[].edges[minipath3[idx-1]].del(l)
    #               break
    #         diff_flag = false
    #         idx += 1
    #       of 1'u8: #ins in OG path; del in after
    #         diff_flag = true
    #         rep_po[].deleted_nodes.incl(minipath3[idx])
    #         rep_po[].weights.del((minipath3[idx-1],minipath3[idx]))
    #         if (minipath3[idx-1],minipath3[idx]) == (17197'u32,17204'u32):
    #           echo "THIS IS BAD4"
    #           print_flag = true
    #         if print_flag:
    #           echo "deleting2 - ", minipath3[idx-1], " -> ", minipath3[idx]
    #         for l,k in rep_po[].edges[minipath3[idx-1]]:
    #           if k == minipath3[idx]:
    #             rep_po[].edges[minipath3[idx-1]].del(l)
    #             break
    #         idx += 1
    #         delete_flag = true
    #       else: #2'u8 del in OG path; ins in after
    #         diff_flag = true
    #         delete_flag = true
    #   if delete_flag:
    #     to_remove.add(minipath)

const IDX_TO_NT = ['A','G','C','T']

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


proc assignReadsToPaths( reads : seq[Read], paths : seq[seq[uint32]]) : seq[int] = 
  var closest_paths : seq[int]
  for read in reads:
    var path_scores : seq[int]
    for path in paths:
      path_scores.add(scorePaths(read.path,path))
    closest_paths.add(maxIdx(path_scores))
  return closest_paths


proc walkHeaviestPaths( po : ptr POGraph,psi = 15) : seq[seq[uint32]] = 
  var representative_paths : seq[seq[uint32]]
  var remaining_reads = toSeq(0..<po[].reads.len)
  var total_st_time = 0.0
  var total_sort_time = 0.0
  var total_walk_time = 0.0
  while remaining_reads.len != 0:
    # echo "here"
    var reads : seq[Read] = @[]
    for i in remaining_reads:
      reads.add(po[].reads[i])
    var time1 = cpuTime()
    let st = getStringtieLikeDataStructures(po,reads)
    total_st_time += (cpuTime() - time1)
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
    time1 = cpuTime()
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
    total_walk_time += (cpuTime() - time1)
    let accounted_for_reads = all_reads - excluded_reads
    time1 = cpuTime()
    representative_paths.add(sorted(representative_path))
    total_sort_time += (cpuTime() - time1)
    var to_delete_reads : seq[int]
    for read_id in accounted_for_reads:
      to_delete_reads.add(read_id_to_idx[read_id])
    for idx in sorted(to_delete_reads,order = SortOrder.Descending):
      remaining_reads.delete(idx)
  echo "ST time - ", total_st_time
  echo "Walk time - ", total_walk_time
  echo "sort time - ", total_sort_time
  return representative_paths