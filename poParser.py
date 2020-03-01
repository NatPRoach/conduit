#!/usr/bin/env python
import sys
import copy
import textwrap
import heapq as pq
import numpy as np
import time
import numba
#
# sys.setrecursionlimit(1500)
class Read():
    def __init__(self, name, length,path = []):
        self.name = name
        self.length = length
        self.path = []
        self.seq = []

class Node():
    def __init__(self,char,align_ring_partner = None,indegree = 0,distance = -np.inf,likelihood = np.inf,start_node_flag = False,end_node_flag=False):
        self.char = char
        self.supporting_reads = []
        self.align_ring_partner = align_ring_partner
        self.visited = False
        self.indegree = 0
        self.distance = distance
        self.log_likelihood = likelihood
        self.log_likelihood_scaled = likelihood
        self.path = []
        self.corrected_path = []
        self.start_node_flag = start_node_flag
        self.end_node_flag = end_node_flag
# class GraphNode():
#     def __init__(self,name,parent,list_):
#         self.name = name
#         self.visited = False
#         self.parent = parent
#         self.list = list_
        

@numba.jit()
def alignSeqs(x,y,a=2,b=4,o=4,e=2,anchor="both"):
    #x - first seq (query for the purposes of defining insertion / deletion)
    #y - second seq (reference for the purposes of defining insertion / deletion)
    #a - match score
    #b - mismatch penalty
    #o - gap open penalty
    #e - gap extend penalty
    inf = np.inf
    # x = x.upper()
    # y = y.upper()
    # print(len(x),len(y))
    
    # time1 = time.time()
    ## Create score and traceback matrices in shape x (down) by y (across)
    
    score_matrix = np.zeros([len(x) + 1, len(y) + 1])
    trace_matrix = np.zeros([len(x) + 1, len(y) + 1])
    if anchor == "both": ## Working as intended
        ## Both end gaps are penalized
        for i in range(1,len(x)+1): 
            score_matrix[i,0] = -o - i*e
        for j in range(1,len(y)+1):
            score_matrix[0,j] = -o - j*e
    col_score = [0]
    for j in range(1,len(y) + 1):
        col_score.append(-2*o - j*e)
    # time2 = time.time()
    # print("Initialization: %s seconds" %( time2 - time1))
    ## Fill in the score matrices:
    total_time = 0
    for i in range(1,len(x)+1):
        row_score = -2 * o - i*e
        for j in range(1,len(y)+1):
            # time3 = time.time()
            if x[i-1] == y[j-1]:
                no_gap_score = score_matrix[i-1,j-1] + a
            else:
                no_gap_score = score_matrix[i-1,j-1] - b
            # no_gap_score = score_matrix[i-1,j-1] + self.s(x[i-1],y[j-1],a,b)
            # m[i,j] = v[i-1,j-1] + self.s(x[i-1],y[j-1],a,b)
            # row_open = score_matrix[i,j-1] - o
            # row_extend = row_score - e
            # row_score = max(row_open,row_extend)
            row_score = max(score_matrix[i,j-1] - o, row_score) - e
            # col_open = score_matrix[i-1,j] - o
            # col_extend = col_score[j] - e
            # col_score[j] = max(col_open,col_extend)
            col_score[j] = max(score_matrix[i-1,j] - o,col_score[j]) - e
            best_score = max(no_gap_score,col_score[j],row_score)
            score_matrix[i,j] = best_score
            # time4 = time.time()
            # total_time += time4 - time3
            
            possible_origins = 0
            if no_gap_score == best_score:
                possible_origins += 1 #Align seqA with seqB
            if col_score[j] == best_score:
                possible_origins += 2
            if row_score == best_score:
                possible_origins += 4
            trace_matrix[i,j] = possible_origins
    # print(score_matrix)
    # time1 = time.time()
    # print("max logic: %s seconds" %(total_time))
    # print("Fill-in: %s seconds" %(time1 - time2))
    ## Walkback
    if anchor == "both":
        i = len(x)
        j = len(y)
        walkback = []

    while i != 0 and j != 0:
        direction = trace_matrix[i,j]
        if direction == 1:
            i = i-1
            j = j-1
            walkback.append(0)
        elif direction == 2:
            i= i-1
            walkback.append(1)
        elif direction == 4:
            j=j-1
            walkback.append(2)
        elif direction == 3: #Ambiguous alignment
            i = i-1 #Arbitrarily decide to use match
            j = j-1
            walkback.append(0)
        elif direction == 5: #Ambiguous alignment
            i = i-1 #Arbitrarily decide to use match
            j = j-1
            walkback.append(0)
        elif direction == 6: #Ambiguous alignment
            j=j-1 #Arbitrarily decide to use deletion
            walkback.append(2)
        elif direction == 7: #Ambiguous alignment
            i = i-1 #Arbitrarily decide to use match
            j = j-1
            walkback.append(0)
    if i > 0:
        walkback =  walkback + [1 for x in range(i)] # 1 is insertion
    elif j > 0:
        walkback =  walkback + [2 for x in range(j)] # 2 is deletion
    walkback.reverse()
    # time2 = time.time()
    # print("Walkback: %s seconds" %(time2 - time1))
    return walkback

class splitProfileHMM(object):
    def __init__(self,nodes,node_indexes,edges,representative_paths,read_seqs,read_ids,source_nodes,sink_nodes,outfile):
        self.nodes = copy.deepcopy(nodes)
        self.nt_idxs = {'A': 0, 'a' : 0,
                        'G': 1, 'g' : 1,
                        'C': 2, 'c' : 2,
                        'T': 3, 't' : 3, 'U': 3, 'u' : 3,
                        0 : 'A',
                        1 : 'G',
                        2 : 'C',
                        3 : 'T'}
        new_edges = {}
        u_set = set()
        v_set = set()
        print("source_nodes1:",source_nodes)
        print("sink_nodes1:",sink_nodes)
        source_nodes2 = []
        sink_nodes2 = []
        for u in edges:
            u_back = node_indexes['b',u]
            if u in source_nodes:
                source_nodes2.append(u_back)
            u_set.add(u_back)
            if len(edges[u]) > 0:
                for v in edges[u]:
                    v_back = node_indexes['b',v]
                    if v in sink_nodes:
                        sink_nodes2.append(v_back)
                    v_set.add(v_back)
                    if u_back in new_edges:
                        new_edges[u_back].append(v_back)
                    else:
                        new_edges[u_back] = [v_back]
            else:
                new_edges[u_back] = []
        # source_nodes = list(u_set - v_set)
        source_nodes = source_nodes2
        print("source_nodes2:",source_nodes)
        # sink_nodes = list(v_set - u_set)
        sink_nodes = sink_nodes2
        print("sink_nodes2:",sink_nodes)
        self.end_nodes = set()
        node_lengths = []
        node_starts = [0]
        for node in self.nodes:
            node_lengths.append(len(node.char))
            node_starts.append(node_starts[-1] + len(node.char))
        del node_starts[-1]
        num_match_states = sum(node_lengths)
        self.num_match_states = num_match_states
        num_delete_states = num_match_states
        self.num_delete_states = num_delete_states
        num_insert_states = num_match_states + len(source_nodes)
        self.num_insert_states = num_insert_states
        num_phmm_nodes = 3*num_match_states + len(source_nodes) + 2 # match, delete, insert for each char; + insert nodes before starts; + start and end nodes
        self.num_phmm_nodes = num_phmm_nodes
        start = num_phmm_nodes - 2
        self.start = start
        end = num_phmm_nodes - 1
        self.end = end
        self.num_source_nodes = len(source_nodes)
        transitions = np.zeros((num_phmm_nodes,num_phmm_nodes))
        # emissions = np.zeros((num_phmm_nodes,4)) # A G C T emission counts for each node
        emissions = np.ones((num_phmm_nodes,4))
        ########################################################
        ### INITIALIZE START NODE TO SOURCE NODE TRANSITIONS ###
        ########################################################
        self.source_node_order = {}
        for i,j in enumerate(source_nodes):
            self.source_node_order[j] = i
            i1 = 3*num_match_states + i
            m2 = node_starts[j]
            d2 = m2 + num_match_states
            transitions[start,m2] += 1 #S  -> M2
            transitions[start,d2] += 1 #S  -> D2
            transitions[start,i1] += 1 #S  -> I1
            transitions[i1,i1] += 1    #I1 -> I1
            transitions[i1,m2] += 1    #I1 -> M2
            transitions[i1,d2] += 1    #I1 -> D2
        
        for i, start in enumerate(node_starts):
            ###############################################
            ### INITIALIZE LINEAR BLOCKS OF TRANSITIONS ###
            ###############################################
            for j in range(start,start + node_lengths[i] - 1):
                m1 = j
                m2 = m1+1 
                d1 = m1+num_match_states
                d2 = d1+1
                i1 = d1+num_match_states
                transitions[m1,m2] += 1 #M1 -> M2
                transitions[m1,d2] += 1 #M1 -> D2
                transitions[m1,i1] += 1 #M1 -> I1
                transitions[d1,m2] += 1 #D1 -> M2
                transitions[d1,d2] += 1 #D1 -> D2
                transitions[d1,i1] += 1 #D1 -> I1
                transitions[i1,m2] += 1 #I1 -> M2
                transitions[i1,d2] += 1 #I1 -> D2
                transitions[i1,i1] += 1 #I1 -> I1
            m1 = start + node_lengths[i] - 1
            d1 = m1 + num_match_states
            i1 = d1 + num_match_states
            transitions[m1,i1] += 1 #M1 -> I1
            transitions[d1,i1] += 1 #D1 -> I1
            transitions[i1,i1] += 1 #I1 -> I1
            
            ###############################################
            ###   INITIALIZE NODE TO NODE TRANSITIONS   ###
            ###############################################
            for j in new_edges[i]:
                m2 = node_starts[j]
                d2 = m2 + num_match_states
                transitions[m1,m2] += 1 #M1 -> M2
                transitions[m1,d2] += 1 #M1 -> D2
                transitions[d1,m2] += 1 #D1 -> M2
                transitions[d1,d2] += 1 #D1 -> D2
                transitions[i1,m2] += 1 #I1 -> M2
                transitions[i1,d2] += 1 #I1 -> D2
            ###############################################
            ### INITIALIZE NODE TO END NODE TRANSITIONS ###
            ###############################################
            # if len(new_edges[i]) == 0:
            if i in sink_nodes:
                transitions[m1,end] += 1
                transitions[d1,end] += 1
                transitions[i1,end] += 1
                self.end_nodes.add(m1)
                self.end_nodes.add(d1)
                self.end_nodes.add(i1)
        # print(transitions)
        ###################################################################
        ### ADD PSUEDO COUNT TO EACH NT FOR EACH INSERT AND MATCH STATE ###
        ###################################################################
        ### Switched initialization to np.ones
        # for i in range(num_match_states):
        #     for j in range(4):
        #         emissions[i,j] += 1
        #
        # # ### For now, add psuedo counts to delete states, and start and end states, for the purposes of not dividing by zero. Wont actually matter. We can change structure of pHMM later.
        # for i in range(num_match_states,num_match_states + num_delete_states):
        #     for j in range(4):
        #         emissions[i,j] += 1
        # for i in range(num_match_states + num_delete_states + num_insert_states, num_match_states + num_delete_states + num_insert_states + 2 ):
        #     for j in range(4):
        #         emissions[i,j] += 1
        #
        # for i in range(num_match_states + num_delete_states, num_match_states + num_delete_states + num_insert_states):
        #     for j in range(4):
        #         emissions[i,j] += 1
        
        ###################################################################
        ###  ADD SELF TRANSITION IN END NODE TO AVOID DIVIDING BY ZERO  ###
        ###################################################################
        ### Find the transitions in to each node:
        in_nodes = {}
        for i in range(num_phmm_nodes):
            for j in range(num_phmm_nodes):
                if transitions[i,j] != 0:
                    if j in in_nodes:
                        in_nodes[j].append(i)
                    else:
                        in_nodes[j] = [i]
        
        self.in_nodes = in_nodes
        transitions[end,end] += 1
        
        self.transitions = transitions
        self.emissions = emissions
        reduced_paths = []
        path_seqs = []
        phmm_paths = []
        for i,path in enumerate(representative_paths):
            path,seq = self.reducePath(path,node_indexes)
            # reduced_paths.append(path)
            # path_seqs.append(seq)
            start_time = time.time()
            # alignment = self.alignSeqs(read_seqs[i], seq)
            # print(read_seqs[i])
            # print(seq)
            alignment = alignSeqs(read_seqs[i], seq)
            # print(seq)
            # print(read_seqs[i])
            print("Alignment: ---%s--- seconds" %(time.time() - start_time))
            # try:
            #     alignment2 = alignSeqs(read_seqs[i], seq)
            #     assert alignment == alignment2
            # except:
            #     raise
                # print(alignment)
                # print(alignment2)
            start_time = time.time()
            phmm_paths.append(self.updateHMM(path,alignment,read_seqs[i],self.source_node_order,num_match_states,node_starts,node_indexes))
            print("Update: ---%s--- seconds" %(time.time() - start_time))
        # print(self.transitions)
        # print("\n\n")
        # print(self.emissions)
        self.transitions = (self.transitions.T / np.sum(self.transitions,axis = 1)).T #Dividing by zero in end node where we're not transitioning anywhere.
        self.emissions = (self.emissions.T / np.sum(self.emissions,  axis = 1)).T #Dividing by zero in delete nodes
        # print(self.transitions)
        # print("\n\n")
        # print(self.emissions)
        
        corrected_seqs = []
        for i in range(len(phmm_paths)):
            corrected_seqs.append(self.correctReads(phmm_paths[i],read_seqs[i],num_match_states,len(source_nodes)))
        self.outputReads(read_ids,corrected_seqs,outfile)
    
    def outputReads(self,read_ids,corrected_seqs,outfile):
        for i in range(len(read_ids)):
            outfile.write('>' + read_ids[i] + '\n')
            outfile.write(corrected_seqs[i] + '\n')
        outfile.close()
    
    def reducePath(self,path,node_indexes):
        new_path = []
        i = 0
        seq_list = []
        # print(path)
        # print(len(path))
        while i < len(path):
        # while i != len(path):
            # print(i)
            # print(path[i])
            node = self.nodes[node_indexes[('b',path[i])]]
            new_path.append(node_indexes[('b',path[i])])
            seq_list.append(node.char)
            i += len(node.char)
        print(new_path)
        return new_path, "".join(seq_list)
    
    def updateHMM(self,path,alignment,read_seq,source_node_order,num_match_states,node_starts,node_indexes):
        # path_index = 0
        seq_idx = 0
        align_idx = 0
        start_idx = 3 * num_match_states + len(source_node_order)
        end_idx = start_idx + 1
        phmm_path = [start_idx]
        last_insert_node = 3 * num_match_states + source_node_order[path[0]]
        debug_dict = {}
        for node_idx in path:
            node_char_count = 0
            node = self.nodes[node_idx]
            print("align_idx:", align_idx)
            while node_char_count < len(node.char):
                m_idx = node_starts[node_idx] + node_char_count
                d_idx = m_idx + num_match_states
                i_idx = d_idx + num_match_states
                if alignment[align_idx] == 0: #Match
                    phmm_path.append(m_idx)
                    debug_dict[m_idx] = ("Match",node_idx,node_indexes[('f',node_idx)],node_char_count)
                    self.emissions[m_idx,self.nt_idxs[read_seq[seq_idx]]] += 1
                    last_insert_node = i_idx
                    seq_idx += 1
                    node_char_count += 1
                elif alignment[align_idx] == 1: # Insertion
                    phmm_path.append(last_insert_node)
                    debug_dict[last_insert_node] = ("Insert",node_idx,node_indexes[('f',node_idx)],node_char_count)
                    self.emissions[last_insert_node,self.nt_idxs[read_seq[seq_idx]]] += 1
                    seq_idx += 1
                elif alignment[align_idx] == 2: # Deletion
                    phmm_path.append(d_idx)
                    debug_dict[d_idx] = ("Delete",node_idx,node_indexes[('f',node_idx)],node_char_count)
                    last_insert_node = i_idx
                    node_char_count += 1
                align_idx += 1
        phmm_path.append(end_idx)
        for i in range(len(phmm_path) - 1):
            # if i != 0:
            #     print(debug_dict[phmm_path[i]])
            try:
                # pass
                assert self.transitions[phmm_path[i],phmm_path[i+1]] != 0
            except:
                print(start_idx)
                print(i,phmm_path[i],sep='\t')
                print(i+1,phmm_path[i+1],sep='\t')
                print(end_idx)
                raise
            self.transitions[phmm_path[i],phmm_path[i+1]] += 1
        return phmm_path
    
    def correctReads(self,phmm_path,seq,num_match_states,num_source_nodes,psi = 10):
        new_seq = []
        seq_idx = 0
        cont_ins = []
        cont_del = []
        
        for node_idx in phmm_path:
            if node_idx < num_match_states: # Match state
                if len(cont_ins) >= psi:
                    for base in cont_ins:
                        new_seq.append(base)
                elif len(cont_del) < psi:
                    for base in cont_del:
                        new_seq.append(base)
                new_seq.append(self.nt_idxs[max(range(4), key = self.emissions[node_idx,:].__getitem__)])
                seq_idx += 1
                cont_ins = []
                cont_del = []
            elif num_match_states <= node_idx < 2*num_match_states: #Delete state
                cont_del.append(self.nt_idxs[max(range(4), key = self.emissions[node_idx - num_match_states,:].__getitem__)])
                if len(cont_ins) >= psi:
                    for base in cont_ins:
                        new_seq.append(base)
                cont_ins = []
            elif 2*num_match_states <= node_idx < 3*num_match_states + num_source_nodes: #Insert state
                cont_ins.append(seq[seq_idx])
                if len(cont_del) < psi:
                    for base in cont_del:
                        new_seq.append(base)
                cont_del = []
                seq_idx += 1
            elif node_idx == 3*num_match_states + num_source_nodes: # Start state
                pass
            else: #End state
                if len(cont_ins) >= psi:
                    for base in cont_ins:
                        new_seq.append(base)
                elif len(cont_del) < psi:
                    for base in cont_del:
                        new_seq.append(base)
                return ''.join(new_seq)
    
    # def alignSeqs(self,x,y,a=2,b=4,o=4,e=2,anchor="both"):
#         #x - first seq (query for the purposes of defining insertion / deletion)
#         #y - second seq (reference for the purposes of defining insertion / deletion)
#         #a - match score
#         #b - mismatch penalty
#         #o - gap open penalty
#         #e - gap extend penalty
#         inf = np.inf
#         # x = x.upper()
#         # y = y.upper()
#         print(len(x),len(y))
#
#         time1 = time.time()
#         ## Create score and traceback matrices in shape x (down) by y (across)
#
#         score_matrix = np.zeros([len(x) + 1, len(y) + 1])
#         trace_matrix = np.zeros([len(x) + 1, len(y) + 1])
#         if anchor == "both": ## Working as intended
#             ## Both end gaps are penalized
#             for i in range(1,len(x)+1):
#                 score_matrix[i,0] = -o - i*e
#             for j in range(1,len(y)+1):
#                 score_matrix[0,j] = -o - j*e
#         col_score = [0]
#         for j in range(1,len(y) + 1):
#             col_score.append(-o - j*e)
#         time2 = time.time()
#         print("Initialization: %s seconds" %( time2 - time1))
#         ## Fill in the score matrices:
#         total_time = 0
#         for i in range(1,len(x)+1):
#             row_score = -o - i*e
#             for j in range(1,len(y)+1):
#                 time3 = time.time()
#                 if x[i-1] == y[j-1]:
#                     no_gap_score = score_matrix[i-1,j-1] + a
#                 else:
#                     no_gap_score = score_matrix[i-1,j-1] - b
#                 # no_gap_score = score_matrix[i-1,j-1] + self.s(x[i-1],y[j-1],a,b)
#                 # m[i,j] = v[i-1,j-1] + self.s(x[i-1],y[j-1],a,b)
#                 # row_open = score_matrix[i,j-1] - o
#                 # row_extend = row_score - e
#                 # row_score = max(row_open,row_extend)
#                 row_score = max(score_matrix[i,j-1] - o, row_score - e)
#                 # col_open = score_matrix[i-1,j] - o
#                 # col_extend = col_score[j] - e
#                 # col_score[j] = max(col_open,col_extend)
#                 col_score[j] = max(score_matrix[i-1,j] - o,col_score[j] - e)
#                 best_score = max(no_gap_score,col_score[j],row_score)
#                 score_matrix[i,j] = best_score
#                 time4 = time.time()
#                 total_time += time4 - time3
#
#                 possible_origins = 0
#                 if no_gap_score == best_score:
#                     possible_origins += 1 #Align seqA with seqB
#                 if row_score == best_score:
#                     possible_origins += 2
#                 if col_score[j] == best_score:
#                     possible_origins += 4
#                 trace_matrix[i,j] = possible_origins
#         time1 = time.time()
#         print("max logic: %s seconds" %(total_time))
#         print("Fill-in: %s seconds" %(time1 - time2))
#         ## Walkback
#         if anchor is None or anchor == "left":
#             i = len(x)
#             #print i
#             optimal_score = max(v[i,:])
#             #print optimal_score
#             for k in range(len(y)+1):
#                 if v[i,k] == optimal_score:
#                     j = k
#                     break
#             walkback = [ 2 for x in range(len(y)- k)]
#         elif anchor == "both" or anchor == "right":
#             i = len(x)
#             j = len(y)
#             walkback = []
#
#         while i != 0 and j != 0:
#             direction = trace_matrix[i,j]
#             if direction == 1:
#                 i = i-1
#                 j = j-1
#                 walkback.append(0)
#             elif direction == 2:
#                 i= i-1
#                 walkback.append(1)
#             elif direction == 4:
#                 j=j-1
#                 walkback.append(2)
#             elif direction == 3: #Ambiguous alignment
#                 i = i-1 #Arbitrarily decide to use match
#                 j = j-1
#                 walkback.append(0)
#             elif direction == 5: #Ambiguous alignment
#                 i = i-1 #Arbitrarily decide to use match
#                 j = j-1
#                 walkback.append(0)
#             elif direction == 6: #Ambiguous alignment
#                 j=j-1 #Arbitrarily decide to use deletion
#                 walkback.append(2)
#             elif direction == 7: #Ambiguous alignment
#                 i = i-1 #Arbitrarily decide to use match
#                 j = j-1
#                 walkback.append(0)
#         if i > 0:
#             walkback =  walkback + [1 for x in range(i)] # 1 is insertion
#         elif j > 0:
#             walkback =  walkback + [2 for x in range(j)] # 2 is deletion
#         walkback.reverse()
#         time2 = time.time()
#         print("Walkback: %s seconds" %(time2 - time1))
#         return walkback

    def alignSeqs(self,x,y,a=2,b=4,o=4,e=2,anchor="both"):
        #x - first seq (query for the purposes of defining insertion / deletion)
        #y - second seq (reference for the purposes of defining insertion / deletion)
        #a - match score
        #b - mismatch penalty
        #o - gap open penalty
        #e - gap extend penalty
        inf = np.inf
        # x = x.upper()
        # y = y.upper()
        print(len(x),len(y))
        time1 = time.time()
        v = np.zeros([len(x) + 1, len(y) + 1])
        m = np.zeros([len(x) + 1, len(y) + 1])
        x_gap = np.zeros([len(x) + 1, len(y) + 1])
        y_gap = np.zeros([len(x) + 1, len(y) + 1])
        walkback_matrix = np.zeros([len(x) + 1, len(y) + 1])

        ### Initialize matrix values
        if anchor is None: #Semi-global alignment, allow gaps at the begining and end of y
            ## Both end gaps are free. (for y only)
            for i in range(1,len(x)+1):
                v[i,0] = -o - i*e
                #v[i,0] = 0. #Redundant, but do it for now
                #x_gap[i,0] = -o - i*e
                y_gap[i,0] = -o - i*e
                #y_gap[i,0] = -inf
                #x_gap[i,0] = -inf
            for j in range(1,len(y)+1):
                #v[0,j] = -o - j*e
                x_gap[0,j] = -o - j*e
                #y_gap[0,j] = -o - j*e
                v[0,j] = 0. #Redundant, but do it for now
                #x_gap[0,j] = -inf
                #y_gap[0,j] = -inf
        elif anchor == "both": ## Working as intended
            ## Both end gaps are penalized
            for i in range(1,len(x)+1):
                v[i,0] = -o - i*e
                x_gap[i,0] = -o - i*e
                y_gap[i,0] = -inf
            for j in range(1,len(y)+1):
                v[0,j] = -o - j*e
                y_gap[0,j] = -o - j*e
                x_gap[0,j] = -inf
        elif anchor == "left" or anchor == "l" or anchor == "L":
            ## Gaps on the left side are penalized
            for i in range(1,len(x)+1):
                v[i,0] = -o - i*e
                #v[i,0] = 0. #Redundant, but do it for now
                x_gap[i,0] = -o - i*e
                #y_gap[i,0] = -o - i*e
                y_gap[i,0] = -inf
                #x_gap[i,0] = -inf
            for j in range(1,len(y)+1):
                v[0,j] = -o - j*e
                #x_gap[0,j] = -o - j*e
                y_gap[0,j] = -o - j*e
                #v[0,j] = 0. #Redundant, but do it for now
                x_gap[0,j] = -inf
                #y_gap[0,j] = -inf
        elif anchor == "right" or anchor == "r" or anchor == "R":
            ## Gaps on the right side are penalized
            for i in range(1,len(x)+1):
                v[i,0] = -o - i*e
                #v[i,0] = 0. #Redundant, but do it for now
                x_gap[i,0] = -o - i*e
                #y_gap[i,0] = -o - i*e
                y_gap[i,0] = -inf
                #x_gap[i,0] = -inf
            for j in range(1,len(y)+1):
                #v[0,j] = -o - j*e
                #x_gap[0,j] = -o - j*e
                y_gap[0,j] = -o - j*e
                v[0,j] = 0. #Redundant, but do it for now
                x_gap[0,j] = -inf
                #y_gap[0,j] = -inf
        time2 = time.time()
        print("Initialization: %s seconds" %( time2 - time1))
        ## Fill in the score matrices:
        total_time = 0
        for i in range(1,len(x)+1):
            for j in range(1,len(y)+1):
                time3 = time.time()
                if x[i-1] == y[j-1]:
                    m[i,j] = v[i-1,j-1] + a
                else:
                    m[i,j] = v[i-1,j-1] - b
                # m[i,j] = v[i-1,j-1] + self.s(x[i-1],y[j-1],a,b)
                time4 = time.time()
                x_gap[i,j] = max( x_gap[i-1,j], v[i-1,j] - o) - e
                y_gap[i,j] = max( y_gap[i,j-1], v[i,j-1] - o) - e
                v[i,j] = max([m[i,j],x_gap[i,j],y_gap[i,j]])
                total_time += time4 - time3
                possible_origins = 0
                if v[i,j] == m[i,j]:
                    possible_origins += 1
                if v[i,j] == x_gap[i,j]:
                    possible_origins += 2
                if v[i,j] == y_gap[i,j]:
                    possible_origins += 4
                walkback_matrix[i,j] = possible_origins
        time1 = time.time()
        # print(v)
        print("max logic: %s seconds" %(total_time))
        print("Fill-in: %s seconds" %(time1 - time2))
        ## Walkback
        if anchor is None or anchor == "left":
            i = len(x)
            #print i
            optimal_score = max(v[i,:])
            #print optimal_score
            for k in range(len(y)+1):
                if v[i,k] == optimal_score:
                    j = k
                    break
            walkback = [ 2 for x in range(len(y)- k)]
        elif anchor == "both" or anchor == "right":
            i = len(x)
            j = len(y)
            walkback = []

        while i != 0 and j != 0:
            direction = walkback_matrix[i,j]
            if direction == 1:
                i = i-1
                j = j-1
                walkback.append(0)
            elif direction == 2:
                i= i-1
                walkback.append(1)
            elif direction == 4:
                j=j-1
                walkback.append(2)
            elif direction == 3: #Ambiguous alignment
                i = i-1 #Arbitrarily decide to use match
                j = j-1
                walkback.append(0)
            elif direction == 5: #Ambiguous alignment
                i = i-1 #Arbitrarily decide to use match
                j = j-1
                walkback.append(0)
            elif direction == 6: #Ambiguous alignment
                j=j-1 #Arbitrarily decide to use deletion
                walkback.append(2)
            elif direction == 7: #Ambiguous alignment
                i = i-1 #Arbitrarily decide to use match
                j = j-1
                walkback.append(0)
        if i > 0:
            walkback =  walkback + [1 for x in range(i)] # 1 is insertion
        elif j > 0:
            walkback =  walkback + [2 for x in range(j)] # 2 is deletion
        walkback.reverse()
        time2 = time.time()
        print("Walkback: %s seconds" %(time2 - time1))
        # print(walkback)
        return walkback

    def s(self,x,y,a,b):
        if x == y:
            return a
        else:
            return -b
    
    def find_windows(self, node_idx, read_idx):
        num_match_states = self.num_match_states
        num_delete_states = self.num_delete_states
        num_insert_states = self.num_insert_states
        num_phmm_nodes = self.num_phmm_nodes
        start = self.start
        end = self.end
        in_nodes = self.in_nodes
        multipath = []
        while node_idx in in_nodes and read_idx != 0 and not multipath:
            if node_idx in self.windows:
                self.windows[node_idx].add(read_idx)
            else:
                self.windows[node_idx] = set([read_idx])
            read_idx -= 1
            multipath = []
            for idx in in_nodes[node_idx]:
                if idx < num_match_states:
                    multipath.append(idx)
            if len(multipath) == 1:
                node_idx = multipath[0]
                multipath = []
            elif len(multipath) == 0:
                break
            # else:
            #     print(node_idx,read_idx,multipath)
        if multipath and node_idx in in_nodes and read_idx != 0:
            for path in multipath:
                self.find_windows(path,read_idx)
    
    def get_window_iterator(self,j,read_len,window = 10):
        s = set()
        for w in self.windows[j]:
            for i in range(w - window, w + window + 1):
                if 1 <= i <= read_len:
                    s.add(i)
        return s
    
    def viterbi_score2(self,read,max_threshold = None,min_threshold=None,window = 10):
        num_match_states = self.num_match_states
        num_delete_states = self.num_delete_states
        num_insert_states = self.num_insert_states
        num_phmm_nodes = self.num_phmm_nodes
        start = self.start
        end = self.end
        num_source_nodes = self.num_source_nodes
        # transitions = self.transitions
        transitions = np.log(self.transitions)
        # emissions = self.emissions
        emissions = np.log(self.emissions / 0.25)
        # for i in range(num_phmm_nodes):
        #     for j in range(num_phmm_nodes):
        #         if transitions[i,j] == 0:
        #             continue
        #         transitions[i,j] = np.log(transitions[i,j])
        # emissions[:num_match_states,:] = np.log(emissions[:num_match_states,:])
        # emissions[2*num_match_states:2*num_match_states + num_insert_states,:] = np.log(emissions[2*num_match_states:2*num_match_states + num_insert_states,:])
        # for i in range(num_match_states):
        #
        # for i in range(2*num_match_states,2*num_match_states + num_insert_states):
        end_nodes = self.end_nodes
        self.windows = {}
        start_time = time.time()
        print(end_nodes)
        for end_node in end_nodes:
            if end_node < num_match_states:
                self.find_windows(end_node,len(read))
        print("Finding Windows: ---%s--- seconds" %(time.time() - start_time))
        # print(self.windows)
        # for key in sorted(list(self.windows)):
        #     print(key,sorted(list(self.windows[key])))
        nt_idxs = self.nt_idxs
        in_nodes = self.in_nodes
        ### Structure the dynamic programming matrix
        dp_matrix = np.zeros((len(read) + 1,self.num_phmm_nodes))
        ### Structure this in the same order as the transitioning matrix
        ### Find diagonals:
        #Transition from start node to initial insert states, 
        for j in range(3 * num_match_states, 3 * num_match_states + num_source_nodes):
            for i in range(1,len(read) + 1):
                #Self-transitions
                # dp_matrix[i,j] = np.log(emissions[j,nt_idxs[read[i-1]]] / 0.25) + dp_matrix[i-1,j] + np.log(transitions[j,j])
                # dp_matrix[i,j] = np.log(emissions[j,nt_idxs[read[i-1]]] / 0.25) + dp_matrix[i-1,j] + transitions[j,j]
                dp_matrix[i,j] = emissions[j,nt_idxs[read[i-1]]] + dp_matrix[i-1,j] + transitions[j,j]
                # if min_threshold is not None:
                #     if dp_matrix[i,j] < min_threshold:
                #         break
        #Iterate through the other nodes, which should all be semi topologically sorted
        for j in range(num_match_states):
            match_idx = j
            delete_idx = num_match_states + j
            insert_idx = 2 * num_match_states + j
            for i in self.get_window_iterator(j,len(read),window):
            # for i in range(1, len(read) + 1):
                # flag = True
                # for w in self.windows[j]:
                #     if abs(i - w) < window:
                #         flag = False
                #         break
                # if flag:
                #     continue
                ### Match state
                # emission_ll = np.log(emissions[match_idx,nt_idxs[read[i - 1]]] / 0.25)
                # emission_ll = emissions[match_idx,nt_idxs[read[i - 1]]]
                lls = []
                max_ll = -np.inf
                for k in in_nodes[match_idx]:
                    # transition_ll = np.log(transitions[k,match_idx])
                    # transition_ll = transitions[k,match_idx]
                    # previous_ll = dp_matrix[i-1, k]
                    # lls.append( previous_ll + transition_ll)
                    # lls.append( dp_matrix[i-1, k] + transitions[k,match_idx])
                    ll = dp_matrix[i-1, k] + transitions[k,match_idx]
                    if ll > max_ll:
                        max_ll = ll
                # dp_matrix[i,match_idx] = emission_ll + max(lls)
                # dp_matrix[i,match_idx] = emissions[match_idx,nt_idxs[read[i - 1]]] + max(lls)
                dp_matrix[i,match_idx] = emissions[match_idx,nt_idxs[read[i - 1]]] + max_ll
                if max_threshold is not None:
                    if dp_matrix[i,match_idx] > max_threshold:
                        print("Ending early!")
                        return dp_matrix
                ### Delete state
                lls = []
                max_ll = -np.inf
                for k in in_nodes[delete_idx]:
                    # transition_ll = np.log(transitions[k,delete_idx])
                    # transition_ll = transitions[k,delete_idx]
                    # previous_ll = dp_matrix[i,k]
                    # lls.append(previous_ll + transition_ll)
                    # lls.append(dp_matrix[i,k] + transitions[k,delete_idx])
                    ll = dp_matrix[i,k] + transitions[k,delete_idx]
                    if ll > max_ll:
                        max_ll = ll
                # dp_matrix[i,delete_idx] = max(lls)
                dp_matrix[i,delete_idx] = max_ll
    
                ### Insert state
                # emission_ll = np.log(emissions[insert_idx,nt_idxs[read[i - 1]]] / 0.25)
                # emission_ll = emissions[insert_idx,nt_idxs[read[i - 1]]]
                lls = []
                max_ll = -np.inf
                for k in in_nodes[insert_idx]:
                    # transition_p = np.log(transitions[k,insert_idx])
                    # transition_p = transitions[k,insert_idx]
                    # previous_ll = dp_matrix[i-1,k]
                    # lls.append(previous_ll + transition_ll)
                    # lls.append(dp_matrix[i-1,k] + transitions[k,insert_idx])
                    ll = dp_matrix[i-1,k] + transitions[k,insert_idx]
                    if ll > max_ll:
                        max_ll = ll
                # dp_matrix[i,insert_idx] = emission_ll + max(lls)
                # dp_matrix[i,insert_idx] = emissions[insert_idx,nt_idxs[read[i - 1]]] + max(lls)
                dp_matrix[i,insert_idx] = emissions[insert_idx,nt_idxs[read[i - 1]]] + max_ll
                if max_threshold is not None:
                    if dp_matrix[i,insert_idx] > max_threshold:
                        print("Ending early!")
                        return dp_matrix
        
        #handle end state:
        i = len(read) #Want to fully consume the read
        lls = []
        max_ll = -np.inf
        for k in in_nodes[end]:
            # transition_ll = np.log(transitions[k,end])
            # transition_ll = transitions[k,end]
            # previous_ll = dp_matrix[i,k]
            # lls.append(previous_ll + transition_ll)
            # lls.append(dp_matrix[i,k] + transitions[k,end])
            ll = dp_matrix[i,k] + transitions[k,end]
            if ll > max_ll:
                max_ll = ll
        dp_matrix[i,end] = max_ll
        print(dp_matrix)
        print(dp_matrix[i,end])
        return dp_matrix
    
    def viterbi_score(self,read,max_threshold = None,min_threshold=None):
        num_match_states = self.num_match_states
        num_delete_states = self.num_delete_states
        num_insert_states = self.num_insert_states
        num_phmm_nodes = self.num_phmm_nodes
        start = self.start
        end = self.end
        num_source_nodes = self.num_source_nodes
        # transitions = self.transitions
        transitions = np.log(self.transitions)
        emissions = np.log(self.emissions / 0.25)
        nt_idxs = self.nt_idxs
        in_nodes = self.in_nodes
        ### Structure the dynamic programming matrix
        dp_matrix = np.zeros((len(read) + 1,self.num_phmm_nodes))
        ### Structure this in the same order as the transitioning matrix
        
        #Transition from start node to initial insert states, 
        for j in range(3 * num_match_states, 3 * num_match_states + num_source_nodes):
            for i in range(1,len(read) + 1):
                #Self-transitions
                # dp_matrix[i,j] = np.log(emissions[j,nt_idxs[read[i-1]]] / 0.25) + dp_matrix[i-1,j] + np.log(transitions[j,j])
                # dp_matrix[i,j] = np.log(emissions[j,nt_idxs[read[i-1]]] / 0.25) + dp_matrix[i-1,j] + transitions[j,j]
                dp_matrix[i,j] = emissions[j,nt_idxs[read[i-1]]] + dp_matrix[i-1,j] + transitions[j,j]
                # if min_threshold is not None:
                #     if dp_matrix[i,j] < min_threshold:
                #         break
        #Iterate through the other nodes, which should all be semi topologically sorted
        
        for j in range(num_match_states):
            match_idx = j
            delete_idx = num_match_states + j
            insert_idx = 2 * num_match_states + j
            for i in range(1, len(read) + 1):
                ### Match state
                # emission_ll = np.log(emissions[match_idx,nt_idxs[read[i - 1]]] / 0.25)
                # emission_ll = emissions[match_idx,nt_idxs[read[i - 1]]]
                lls = []
                max_ll = -np.inf
                for k in in_nodes[match_idx]:
                    # transition_ll = np.log(transitions[k,match_idx])
                    # transition_ll = transitions[k,match_idx]
                    # previous_ll = dp_matrix[i-1, k]
                    # lls.append( previous_ll + transition_ll)
                    # lls.append( dp_matrix[i-1, k] + transitions[k,match_idx])
                    ll = dp_matrix[i-1, k] + transitions[k,match_idx]
                    if ll > max_ll:
                        max_ll = ll
                # dp_matrix[i,match_idx] = emission_ll + max(lls)
                # dp_matrix[i,match_idx] = emissions[match_idx,nt_idxs[read[i - 1]]] + max(lls)
                dp_matrix[i,match_idx] = emissions[match_idx,nt_idxs[read[i - 1]]] + max_ll
                if max_threshold is not None:
                    if dp_matrix[i,match_idx] > max_threshold:
                        print("Ending early!")
                        return dp_matrix
                ### Delete state
                lls = []
                max_ll = -np.inf
                for k in in_nodes[delete_idx]:
                    # transition_ll = np.log(transitions[k,delete_idx])
                    # transition_ll = transitions[k,delete_idx]
                    # previous_ll = dp_matrix[i,k]
                    # lls.append(previous_ll + transition_ll)
                    # lls.append(dp_matrix[i,k] + transitions[k,delete_idx])
                    ll = dp_matrix[i,k] + transitions[k,delete_idx]
                    if ll > max_ll:
                        max_ll = ll
                # dp_matrix[i,delete_idx] = max(lls)
                dp_matrix[i,delete_idx] = max_ll
                
                ### Insert state
                # emission_ll = np.log(emissions[insert_idx,nt_idxs[read[i - 1]]] / 0.25)
                # emission_ll = emissions[insert_idx,nt_idxs[read[i - 1]]]
                lls = []
                max_ll = -np.inf
                for k in in_nodes[insert_idx]:
                    # transition_p = np.log(transitions[k,insert_idx])
                    # transition_p = transitions[k,insert_idx]
                    # previous_ll = dp_matrix[i-1,k]
                    # lls.append(previous_ll + transition_ll)
                    # lls.append(dp_matrix[i-1,k] + transitions[k,insert_idx])
                    ll = dp_matrix[i-1,k] + transitions[k,insert_idx]
                    if ll > max_ll:
                        max_ll = ll
                # dp_matrix[i,insert_idx] = emission_ll + max(lls)
                # dp_matrix[i,insert_idx] = emissions[insert_idx,nt_idxs[read[i - 1]]] + max(lls)
                dp_matrix[i,insert_idx] = emissions[insert_idx,nt_idxs[read[i - 1]]] + max_ll
                if max_threshold is not None:
                    if dp_matrix[i,insert_idx] > max_threshold:
                        print("Ending early!")
                        return dp_matrix
                    
        #handle end state:
        i = len(read) #Want to fully consume the read
        lls = []
        max_ll = -np.inf
        for k in in_nodes[end]:
            # transition_ll = np.log(transitions[k,end])
            # transition_ll = transitions[k,end]
            # previous_ll = dp_matrix[i,k]
            # lls.append(previous_ll + transition_ll)
            # lls.append(dp_matrix[i,k] + transitions[k,end])
            ll = dp_matrix[i,k] + transitions[k,end]
            if ll > max_ll:
                max_ll = ll
        
        dp_matrix[i,end] = max_ll
        print(dp_matrix)
        print(dp_matrix[i,end])
        return dp_matrix
    
class poGraph():
    def __init__(self):
        self.nodes = []
        self.reads = []
        self.visited = []
        self.edges = {}
        self.weights = {}
        self.source_node_indices = []
        self.end_node_indices = []
        self.vertex_list = []
        self.read_paths = []
        self.trimmed_edges = {}
        # self.trimmed_weights = {}
        self.log_likelihoods = {}
    
    def alignSeqs(self,x,y,a=2,b=4,o=4,e=2,anchor="both"):
        #x - first seq (query for the purposes of defining insertion / deletion)
        #y - second seq (reference for the purposes of defining insertion / deletion)
        #a - match score
        #b - mismatch penalty
        #o - gap open penalty
        #e - gap extend penalty
        inf = np.inf
        # x = x.upper()
        # y = y.upper()
        v = np.zeros([len(x) + 1, len(y) + 1])
        m = np.zeros([len(x) + 1, len(y) + 1])
        x_gap = np.zeros([len(x) + 1, len(y) + 1])
        y_gap = np.zeros([len(x) + 1, len(y) + 1])
        walkback_matrix = np.zeros([len(x) + 1, len(y) + 1])
    
        ### Initialize matrix values
        if anchor is None: #Semi-global alignment, allow gaps at the begining and end of y
            ## Both end gaps are free. (for y only)
            for i in range(1,len(x)+1): 
                v[i,0] = -o - i*e
                #v[i,0] = 0. #Redundant, but do it for now
                #x_gap[i,0] = -o - i*e
                y_gap[i,0] = -o - i*e
                #y_gap[i,0] = -inf
                #x_gap[i,0] = -inf
            for j in range(1,len(y)+1): 
                #v[0,j] = -o - j*e
                x_gap[0,j] = -o - j*e
                #y_gap[0,j] = -o - j*e
                v[0,j] = 0. #Redundant, but do it for now
                #x_gap[0,j] = -inf
                #y_gap[0,j] = -inf
        elif anchor == "both": ## Working as intended
            ## Both end gaps are penalized
            for i in range(1,len(x)+1): 
                v[i,0] = -o - i*e
                x_gap[i,0] = -o - i*e
                y_gap[i,0] = -inf
            for j in range(1,len(y)+1):
                v[0,j] = -o - j*e
                y_gap[0,j] = -o - j*e
                x_gap[0,j] = -inf
        elif anchor == "left" or anchor == "l" or anchor == "L":
            ## Gaps on the left side are penalized
            for i in range(1,len(x)+1): 
                v[i,0] = -o - i*e
                #v[i,0] = 0. #Redundant, but do it for now
                x_gap[i,0] = -o - i*e
                #y_gap[i,0] = -o - i*e
                y_gap[i,0] = -inf
                #x_gap[i,0] = -inf
            for j in range(1,len(y)+1): 
                v[0,j] = -o - j*e
                #x_gap[0,j] = -o - j*e
                y_gap[0,j] = -o - j*e
                #v[0,j] = 0. #Redundant, but do it for now
                x_gap[0,j] = -inf
                #y_gap[0,j] = -inf
        elif anchor == "right" or anchor == "r" or anchor == "R":
            ## Gaps on the right side are penalized
            for i in range(1,len(x)+1): 
                v[i,0] = -o - i*e
                #v[i,0] = 0. #Redundant, but do it for now
                x_gap[i,0] = -o - i*e
                #y_gap[i,0] = -o - i*e
                y_gap[i,0] = -inf
                #x_gap[i,0] = -inf
            for j in range(1,len(y)+1): 
                #v[0,j] = -o - j*e
                #x_gap[0,j] = -o - j*e
                y_gap[0,j] = -o - j*e
                v[0,j] = 0. #Redundant, but do it for now
                x_gap[0,j] = -inf
                #y_gap[0,j] = -inf
        ## Fill in the score matrices:
        for i in range(1,len(x)+1):
            for j in range(1,len(y)+1):
                m[i,j] = v[i-1,j-1] + self.s(x[i-1],y[j-1],a,b)
                x_gap[i,j] = max( x_gap[i-1,j], v[i-1,j] - o) - e
                y_gap[i,j] = max( y_gap[i,j-1], v[i,j-1] - o) - e
                v[i,j] = max([m[i,j],x_gap[i,j],y_gap[i,j]])
                possible_origins = 0
                if v[i,j] == m[i,j]:
                    possible_origins += 1
                if v[i,j] == x_gap[i,j]:
                    possible_origins += 2
                if v[i,j] == y_gap[i,j]:
                    possible_origins += 4
                walkback_matrix[i,j] = possible_origins
        ## Walkback
        if anchor is None or anchor == "left":
            i = len(x)
            #print i
            optimal_score = max(v[i,:])
            #print optimal_score
            for k in range(len(y)+1):
                if v[i,k] == optimal_score:
                    j = k
                    break
            #print "j"
            #print j
            walkback = [ 2 for x in range(len(y)- k)] 
        elif anchor == "both" or anchor == "right":
            i = len(x)
            j = len(y)
            walkback = []
    
        while i != 0 and j != 0:
            direction = walkback_matrix[i,j]
            if direction == 1:
                i = i-1
                j = j-1
                walkback.append(0)
            elif direction == 2:
                i= i-1
                walkback.append(1)
            elif direction == 4:
                j=j-1
                walkback.append(2)
            elif direction == 3: #Ambiguous alignment
                i = i-1 #Arbitrarily decide to use match
                j = j-1
                walkback.append(0)
            elif direction == 5: #Ambiguous alignment
                i = i-1 #Arbitrarily decide to use match
                j = j-1
                walkback.append(0)
            elif direction == 6: #Ambiguous alignment
                j=j-1 #Arbitrarily decide to use deletion
                walkback.append(2)
            elif direction == 7: #Ambiguous alignment
                i = i-1 #Arbitrarily decide to use match
                j = j-1
                walkback.append(0)
        if i > 0:
            walkback =  walkback + [1 for x in range(i)] # 1 is insertion
        elif j > 0:
            walkback =  walkback + [2 for x in range(j)] # 2 is deletion
        walkback.reverse()
        # print(walkback)
        return walkback
    
    def makeGraph(self,file_object):
        paths = 1
        for _ in range(3):
            _ = file_object.readline()
        num_nodes = int(file_object.readline().strip().split('=')[1])
        num_reads = int(file_object.readline().strip().split('=')[1])
        # print(num_nodes)
        # print(num_reads)
        for _ in range(num_reads):
            name = file_object.readline().strip().split('=')[1]
            length = int(file_object.readline().strip().split('=')[1].split(' ')[0])
            self.reads.append(Read(name,length))
        
        for i in range(num_nodes):
            char,info = file_object.readline().strip().split(':')
            a_pos = info.find('A')
            if a_pos != -1:
                a_pair = int(info[a_pos+1:])
                # print(i,a_pair)
                info = info[:a_pos]
            else:
                a_pair = None
            self.nodes.append(Node(char,align_ring_partner = a_pair))
            s_split = info.split('S')
            info = s_split[0]
            seq_idxs = [int(j) for j in s_split[1:]]
            for seq_idx in seq_idxs:
                if len(self.reads[seq_idx].path) == 0:
                    # print("here")
                    self.nodes[i].start_node_flag = True
                    self.source_node_indices.append(i)
                self.reads[seq_idx].path.append(i)
                self.reads[seq_idx].seq.append(char)
                self.nodes[i].supporting_reads.append(seq_idx)
            if info == '':
                # self.source_node_indices.append(i)
                pass
            else:
                l_split = info.split('L')
                in_nodes = [int(j) for j in l_split[1:]]
                paths = paths * len(in_nodes)
                for in_node in in_nodes:
                    if in_node in self.edges:
                        self.edges[in_node].append(i)
                        # self.edges[in_node][i] = 0
                    else:
                        self.edges[in_node]  = [i]
                        # self.edges[in_node] = {i:0}
                    self.weights[(in_node,i)] = 0
                    # self.visited.append(False)
        for read in self.reads:
            self.nodes[read.path[-1]].end_node_flag = True
            self.end_node_indices.append(read.path[-1])
        # print(self.end_node_indices)
        for i in range(len(self.nodes)):
            if i not in self.edges:
                # self.end_node_indices.append(i)
                self.edges[i] = []
            else:
                for j in self.edges[i]:
                    self.nodes[j].indegree += 1
        for read in self.reads:
            # print(read.path)
            for i in range(len(read.path) - 1):
                # self.edges[read.path[i]][read.path[i+1]] += 1
                self.weights[(read.path[i],read.path[i+1])] += 1
        # print(paths)
    
    def getPossiblePaths(self):
        possible_paths = []
        for source in self.source_node_indices:
            source_paths = self.dfs_paths(source,self.end_node_indices)
            # print(list(source_paths))
            for path in source_paths:
                possible_paths.append(path)
        return possible_paths
    
    def dfs_paths(self, start, goals):
        stack = [(start, [start])]
        counter = 0
        while stack:
            (vertex, path) = stack.pop()
            # print(vertex,path)
            for next_ in set(self.edges[vertex]) - set(path):
                if next_ in goals:
                    yield path + [next_]
                    counter += 1
                    print(counter,end='\r')
                else:
                    stack.append((next_, path + [next_]))
    
    def dfs_paths_one_return(self, start, goals):
        stack = [(start, [start])]
        counter = 0
        while stack:
            (vertex, path) = stack.pop()
            for next_ in set(self.edges[vertex]) - set(path):
                if next_ in goals:
                    yield path + [next_]
                    counter += 1
                    print(counter,end='\r')
                else:
                    stack.append((next_, path + [next_]))
        return []
    
    def getLongestPath(self):
        possible_paths = []
        for source in self.source_node_indices:
            self.nodes[source].distance = 0
            self.nodes[source].path = [source]
        for i,node in enumerate(self.nodes):
            if node.distance is not -np.inf:
                for j in self.edges[i]:
                    node2 = self.nodes[j]
                    if node2.distance < node.distance + self.weights[(i,j)]:
                        node2.distance = node.distance + self.weights[(i,j)]
                        node2.path = node.path + [j]
        longest_path_tuples = []
        for node in self.nodes:
            longest_path_tuples.append((node.distance,node.path))
            # print(node.distance)
            # print(node.path)
        return sorted(longest_path_tuples,key=lambda x : x[0],reverse=True)[0][1]
    
    def getLongestPath2(self,nodes,node_indexes,source_nodes,edges,weights):
        for source in source_nodes:
            nodes[source].distance = 0
            nodes[source].path = [source]
        for i,node in enumerate(nodes):
            k = node_indexes[('f',i)]
            if node.distance is not -np.inf:
                for j in edges[k]:
                    node2 = nodes[node_indexes[('b',j)]]
                    if node2.distance < node.distance + weights[(k,j)]:
                        node2.distance = node.distance + weights[(k,j)]
                        node2.path = node.path + [j]
        longest_path_tuples = []
        for node in nodes:
            longest_path_tuples.append((node.distance,node.path))
            # print(node.distance)
            # print(node.path)
        return sorted(longest_path_tuples,key=lambda x : x[0],reverse=True)[0][1]
    
    def getMostLikelyPath(self,nodes,node_indexes,source_nodes,end_nodes,edges,log_likelihoods):
        print(source_nodes)
        for node in nodes:
            node.log_likelihood = np.inf
            node.log_likelihood_scaled = np.inf
            node.path = []
        for source in source_nodes:
            # print(node_indexes[('b',source)])
            nodes[node_indexes[('b',source)]].log_likelihood = 0
            nodes[node_indexes[('b',source)]].log_likelihood_scaled = 0
            nodes[node_indexes[('b',source)]].path = [source]
        for i,node in enumerate(nodes):
            k = node_indexes[('f',i)]
            if node.log_likelihood is not np.inf:
                for j in edges[k]:
                    node2 = nodes[node_indexes[('b',j)]]
                    if node2.log_likelihood_scaled > (node.log_likelihood + log_likelihoods[(k,j)]) - np.log(len(node.path) + 1):
                    # if node2.log_likelihood > node.log_likelihood + log_likelihoods[(k,j)]:
                        node2.log_likelihood = node.log_likelihood + log_likelihoods[(k,j)]
                        node2.path = node.path + [j]
                        node2.log_likelihood_scaled = (node.log_likelihood + log_likelihoods[(k,j)]) - np.log(len(node.path))
            # if 50 <= k <= 55:
            #     print(k)
            #     for j in edges[k]:
            #         print((k,j),'-',log_likelihoods[(k,j)])
            #         print("node",k,'-',node.log_likelihood)
            #         print("node",j,'-',nodes[node_indexes[('b',j)]].log_likelihood)
            #     print(node.log_likelihood_scaled)
                # print(node.path)
        most_likely_path_tuples = []
        for node_id in end_nodes:
            node = nodes[node_indexes[('b',node_id)]]
            most_likely_path_tuples.append((node.log_likelihood_scaled,node.path))
        # print(most_likely_path_tuples)
        # print(len(end_nodes))
        # print(len(source_nodes))
        # print(sorted(most_likely_path_tuples,key=lambda x : x[0])[0][1])
        return sorted(most_likely_path_tuples,key=lambda x : x[0])[0][1]
        
    def getMostLikelyPaths(self,nodes,node_indexes,source_nodes,end_nodes,edges,log_likelihoods):
        most_likely_path_tuples = []
        # if len(source_nodes) == 1:
        #     print(edges)
        # print(source_nodes)
        for source in source_nodes:
            for node in nodes:
                node.log_likelihood = np.inf
                node.log_likelihood_scaled = np.inf
                node.path = []
            nodes[node_indexes[('b',source)]].log_likelihood = 0
            nodes[node_indexes[('b',source)]].log_likelihood_scaled = 0
            nodes[node_indexes[('b',source)]].path = [source]
            for i,node in enumerate(nodes):
                k = node_indexes[('f',i)]
                # if k == source:
                #      print("HERE")
                if node.log_likelihood is not np.inf:
                    # if len(source_nodes) == 1:
                        # print(k)
                        # print(edges[k])
                    for j in edges[k]:
                        node2 = nodes[node_indexes[('b',j)]]
                        # if node2.log_likelihood > node.log_likelihood + log_likelihoods[(k,j)] :
                        if node2.log_likelihood_scaled > (node.log_likelihood + log_likelihoods[(k,j)]) - np.log(len(node.path) + 1):
                            node2.log_likelihood = node.log_likelihood + log_likelihoods[(k,j)]
                            node2.path = node.path + [j]
                            node2.log_likelihood_scaled = node2.log_likelihood - np.log(len(node2.path))
            for node_id in end_nodes:
                node = nodes[node_indexes[('b',node_id)]]
                most_likely_path_tuples.append((node.log_likelihood,node.path))
        # print(most_likely_path_tuples)
        # print(len(end_nodes))
        # print(len(source_nodes))
        return [ x[1] for x in sorted(most_likely_path_tuples,key=lambda t : t[0])]
    
    def getMostLikelyPaths2(self,nodes,node_indexes,source_nodes,end_nodes,edges,log_likelihoods,psi=10):
        most_likely_path_tuples = []
        # if len(source_nodes) == 1:
        #     print(edges)
        # print(source_nodes)
        if 5035 in end_nodes:
            print("WHAT THE ACTUAL FUCK")
        for source in source_nodes:
            for node in nodes:
                node.log_likelihood = np.inf
                node.log_likelihood_scaled = np.inf
                node.path = []
            nodes[node_indexes[('b',source)]].log_likelihood = 0
            nodes[node_indexes[('b',source)]].log_likelihood_scaled = 0
            nodes[node_indexes[('b',source)]].path = [source]
            for i,node in enumerate(nodes):
                k = node_indexes[('f',i)]
                # if k == source:
                #      print("HERE")
                if node.log_likelihood is not np.inf:
                    # if len(source_nodes) == 1:
                        # print(k)
                        # print(edges[k])
                    for j in edges[k]:
                        node2 = nodes[node_indexes[('b',j)]]
                        # if node2.log_likelihood > node.log_likelihood + log_likelihoods[(k,j)] :
                        if node2.log_likelihood_scaled > (node.log_likelihood + log_likelihoods[(k,j)]) - np.log(len(node.path) + 1):
                            node2.log_likelihood = node.log_likelihood + log_likelihoods[(k,j)]
                            node2.path = node.path + [j]
                            node2.log_likelihood_scaled = node2.log_likelihood - np.log(len(node2.path))
            for node_id in end_nodes:
                node = nodes[node_indexes[('b',node_id)]]
                most_likely_path_tuples.append((node.log_likelihood,node.path))
        # print(most_likely_path_tuples)
        # print(len(end_nodes))
        # print(len(source_nodes))
        different_paths = []
        most_likely_paths = [ x[1] for x in sorted(most_likely_path_tuples,key=lambda t : t[0])]
        for i, path in enumerate(most_likely_paths):
            if i == 0:
                different_paths.append(path)
            else:
                flag = True
                for path2 in different_paths:
                    flag2,__ = self.correctPathsAndEnds(path,path2,psi=psi)
                    flag = flag and flag2
                if flag:
                    different_paths.append(path)
        return different_paths
    
    def getRepresentativePaths2(self,outfile1=None,outfile2=None,output_every_iteration = False,outfile_prefix=None,psi=10):
        if output_every_iteration:
            if outfile_prefix is not None:
                outfile_prefix = "out"
            else:
                print("No outfile prefix provided, defaulting to out")
        if outfile1 is not None:
            read_paths = []
        for i in range(len(self.reads)):
            self.reads[i].corrected_path = self.reads[i].path.copy()
            if outfile1 is not None:
                read_paths.append(self.reads[i].path)
        # print(len(self.reads))
        if outfile1 is not None:
            nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructGraphFromPaths(read_paths)
            del read_paths
            # self.popBubbles(nodes,node_indexes,edges,weights,log_likelihoods)
            self.htmlOutput(outfile1,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods)
            
        representative_paths = []
        differed_reads = list(range(len(self.reads)))
        index = 0
        flag = False
        last_path_considered_idx = 0
        iteration = 0
        while differed_reads:
            new_differed_reads = []
            # nodes,source_nodes,edges,weights = self.constructNewGraph()
            reads = []
            for i in differed_reads:
                reads.append(copy.deepcopy(self.reads[i]))
            nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructNewGraph(reads) ## This is changing every time because we're using the corrected paths.
            
            if output_every_iteration:
                for i in  differed_reads:
                    outfile = open(outfile_prefix + "_Iteration%d_Read%d_Before.html" %(iteration,i),'w')
                    self.htmlOutput(outfile,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path = self.reads[i].corrected_path)
            
            # if index != 0:
            #     assert (last_nodes,last_node_indexes,last_source_nodes,last_edges,last_weights) == (nodes,node_indexes,source_nodes,edges,weights)
            # last_nodes,last_node_indexes,last_source_nodes,last_edges,last_weights = nodes,node_indexes,source_nodes,edges,weights 
            # new_rep_path = self.getLongestPath2(nodes,node_indexes,source_nodes,edges,weights)
            new_rep_paths = self.getMostLikelyPaths(nodes,node_indexes,source_nodes,end_nodes,edges,log_likelihoods)
            # representative_paths.append(new_rep_path)
            # if debug:
            #     return new_rep_paths,reads[1],reads[10]
            if flag:
            # if index != 0 and new_rep_path == representative_paths[index - 1]:
                #Do what? - For now select a random read's path and use that as new possibility (need to fix this in future with better idea)
                representative_paths.append(self.reads[differed_reads[0]].corrected_path)
                # print("here!")
            else:
                for path in new_rep_paths:
                    representative_paths.append(path)
                # if output_every_iteration:
                #     nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructGraphFromPaths(representative_paths)
                #     for j,read in enumerate(reads):
                #         outfile = open(outfile_prefix + "_Iteration%d_Read%d.html" %(iteration,j),'w')
                #         self.htmlOutput(outfile,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path = read.corrected_path)
                    
                    
            # print(representative_paths[index])
            # print(differed_reads)
            # if len(differed_reads) == 1:
            #     if 1618 in edges:
            #         print(edges[1618])
            #     if ('f',0) in node_indexes:
            #         f = node_indexes[('f',0)]
            #         print(f)
            #     if f in edges:
            #         print(edges[f])
            #     print(representative_paths[-1])
            #     print(self.reads[differed_reads[0]].corrected_path)
            for i in differed_reads:
                flag2 = True
                for j in range(last_path_considered_idx,len(representative_paths)):
                    differed,self.reads[i].corrected_path = self.correctPaths(self.reads[i].corrected_path,representative_paths[j],psi= psi)
                # print(differed)
                # if differed:
                    flag2 = differed and flag2
                print(flag2)
                if flag2:
                    new_differed_reads.append(i)
            if output_every_iteration:
                for i in differed_reads:
                    outfile = open(outfile_prefix + "_Iteration%d_Read%d_After.html" %(iteration,i),'w')
                    self.htmlOutput(outfile,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path = self.reads[i].corrected_path)
            # for k,i in enumerate(differed_reads):
            #     if reads[k].corrected_path == self.reads[i].corrected_path:
            #         print("NOT CORRECTING")
            if differed_reads == new_differed_reads:
                flag = True
            else:
                differed_reads = new_differed_reads
                flag = False
            # index += 1
            print(differed_reads)
            # print(index)
            iteration += 1
            last_path_considered_idx = len(representative_paths)
        if outfile2 is not None:
            nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructGraphFromPaths(representative_paths)
            # print(source_nodes)
            # print(node_indexes)
            # print(edges)
            self.htmlOutput(outfile2,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods)
        return representative_paths
    
    def getRepresentativePaths3(self,outfile1=None,outfile2=None,output_every_iteration = False,outfile_prefix=None,psi=10):
        if output_every_iteration:
            if outfile_prefix is not None:
                outfile_prefix = "out"
            else:
                print("No outfile prefix provided, defaulting to out")
        if outfile1 is not None:
            read_paths = []
        for i in range(len(self.reads)):
            self.reads[i].corrected_path = self.reads[i].path.copy()
            if outfile1 is not None:
                read_paths.append(self.reads[i].path)
        # print(len(self.reads))
        if outfile1 is not None:
            nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructGraphFromPaths(read_paths)
            del read_paths
            # self.popBubbles(nodes,node_indexes,edges,weights,log_likelihoods)
            self.htmlOutput(outfile1,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods)
            
        representative_paths = []
        differed_reads = list(range(len(self.reads)))
        index = 0
        flag = False
        last_path_considered_idx = 0
        iteration = 0
        while differed_reads:
            new_differed_reads = []
            # nodes,source_nodes,edges,weights = self.constructNewGraph()
            reads = []
            for i in differed_reads:
                # print(type(self.reads[i]))
                reads.append(copy.deepcopy(self.reads[i]))
            nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructNewGraph(reads) ## This is changing every time because we're using the corrected paths.
            
            if output_every_iteration:
                for i in  differed_reads:
                    outfile = open(outfile_prefix + "_Iteration%d_Read%d_Before.html" %(iteration,i),'w')
                    self.htmlOutput(outfile,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path = self.reads[i].corrected_path)
            
            # if index != 0:
            #     assert (last_nodes,last_node_indexes,last_source_nodes,last_edges,last_weights) == (nodes,node_indexes,source_nodes,edges,weights)
            # last_nodes,last_node_indexes,last_source_nodes,last_edges,last_weights = nodes,node_indexes,source_nodes,edges,weights 
            # new_rep_path = self.getLongestPath2(nodes,node_indexes,source_nodes,edges,weights)
            if 5035 in end_nodes:
                print("WHY IS THIS HAPPENING")
            else:
                print("Hmm")
            new_rep_paths = self.getMostLikelyPaths2(nodes,node_indexes,source_nodes,end_nodes,edges,log_likelihoods,psi=psi)
            # representative_paths.append(new_rep_path)
            # if debug:
            #     return new_rep_paths,reads[1],reads[10]
            if flag:
            # if index != 0 and new_rep_path == representative_paths[index - 1]:
                #Do what? - For now select a random read's path and use that as new possibility (need to fix this in future with better idea)
                representative_paths.append(self.reads[differed_reads[0]].corrected_path)
                if self.reads[differed_reads[0]].corrected_path[-1] == 5035:
                    print("huh") 
                # print("here!")
            else:
                for path in new_rep_paths:
                    if path[-1] == 5035:
                        print("????")
                    representative_paths.append(path)
                # if output_every_iteration:
                #     nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructGraphFromPaths(representative_paths)
                #     for j,read in enumerate(reads):
                #         outfile = open(outfile_prefix + "_Iteration%d_Read%d.html" %(iteration,j),'w')
                #         self.htmlOutput(outfile,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path = read.corrected_path)
                    
                    
            # print(representative_paths[index])
            # print(differed_reads)
            # if len(differed_reads) == 1:
            #     if 1618 in edges:
            #         print(edges[1618])
            #     if ('f',0) in node_indexes:
            #         f = node_indexes[('f',0)]
            #         print(f)
            #     if f in edges:
            #         print(edges[f])
            #     print(representative_paths[-1])
            #     print(self.reads[differed_reads[0]].corrected_path)
            for i in differed_reads:
                flag2 = True
                for j in range(last_path_considered_idx,len(representative_paths)):
                    differed,self.reads[i].corrected_path = self.correctPaths(self.reads[i].corrected_path,representative_paths[j],psi=psi)
                # print(differed)
                # if differed:
                    flag2 = differed and flag2
                print(flag2)
                if flag2:
                    new_differed_reads.append(i)
            if output_every_iteration:
                for i in differed_reads:
                    outfile = open(outfile_prefix + "_Iteration%d_Read%d_After.html" %(iteration,i),'w')
                    self.htmlOutput(outfile,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path = self.reads[i].corrected_path)
            # for k,i in enumerate(differed_reads):
            #     if reads[k].corrected_path == self.reads[i].corrected_path:
            #         print("NOT CORRECTING")
            if differed_reads == new_differed_reads:
                flag = True
            else:
                differed_reads = new_differed_reads
                flag = False
            # index += 1
            print(differed_reads)
            # print(index)
            iteration += 1
            last_path_considered_idx = len(representative_paths)
        if outfile2 is not None:
            nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructGraphFromPaths(representative_paths)
            # print(source_nodes)
            # print(node_indexes)
            # print(edges)
            self.htmlOutput(outfile2,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods)
        print(len(representative_paths))
        return representative_paths
    
    def greedyPathWalk(self,nodes,node_indexes,source_nodes,edges,weights):
        greedy_walks = []
        for fwd_node_idx in source_nodes:
            path = [fwd_node_idx]
            bck_node_idx = node_indexes[('b',fwd_node_idx)]
            node = nodes[bck_node_idx]
            while len(edges[fwd_node_idx]) != 0:
                weight_list = []
                for j in edges[fwd_node_idx]:
                    weight_list.append(weights[(fwd_node_idx,j)])
                fwd_node_idx = edges[fwd_node_idx][max(range(len(weight_list)),key=weight_list.__getitem__)]
                path.append(fwd_node_idx)
            greedy_walks.append(path)
        return greedy_walks
    
    def getRepresentativePaths4(self,outfile1=None,outfile2=None, psi=10):
        time1 = time.time()
        greedy_walks = []
        for node in self.nodes:
            node.visited = False
            node.path = []
        
        reads = []
        for i in range(len(self.reads)):
            # print(type(self.reads[i]))
            self.reads[i].corrected_path = self.reads[i].path.copy()
            reads.append(copy.deepcopy(self.reads[i]))
        nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructNewGraph(reads)
        # print(source_nodes)
        ### collect greedy walks from each possible source node, build priority queue of possible branches not walked. 
        for i,fwd_node_idx in enumerate(source_nodes):
            # print("fwd_node_idx",fwd_node_idx)
            path = [fwd_node_idx]
            bck_node_idx = node_indexes[('b',fwd_node_idx)]
            node = nodes[bck_node_idx]
            node.visited = True
            while len(edges[fwd_node_idx]) != 0:
                # nodes[bck_node_idx].path = path
                weight_list = []
                for j in edges[fwd_node_idx]:
                    weight_list.append(weights[(fwd_node_idx,j)])
                # last_fwd_node_idx = fwd_node_idx
                fwd_node_idx = edges[fwd_node_idx][max(range(len(weight_list)),key=weight_list.__getitem__)]
                bck_node_idx = node_indexes[('b',fwd_node_idx)]
                nodes[bck_node_idx].visited = True
                path.append(fwd_node_idx)
            # nodes[bck_node_idx].path = path
            greedy_walks.append(path)
            # print(greedy_walks[i])
            # for idx in greedy_walks[i]:
            #     print(idx)
        ### correct reads based on these greedy walks.
        differed_reads = []
        for i in range(len(self.reads)):
            flag2 = True
            for j in range(len(greedy_walks)):
                differed,self.reads[i].corrected_path = self.correctPathsAndEnds3(self.reads[i].corrected_path,greedy_walks[j],nodes,node_indexes,psi=psi)
                # differed,self.reads[i].corrected_path = self.correctPathsAndEnds2(self.reads[i].corrected_path,greedy_walks[j],psi=psi)
                # differed,self.reads[i].corrected_path = self.correctPaths(self.reads[i].corrected_path,greedy_walks[j],psi=psi)
            # print(differed)
            # if differed:
                flag2 = differed and flag2
            # print(flag2)
            # for j in self.reads[i].corrected_path:
            #     print(j)
            if flag2:
                differed_reads.append(i)
        # print(time.time() - time1, " seconds")
        ### reconstruct graph based on the corrected walks, keeping track of deleted nodes
        reads = []
        for i in differed_reads:
            # print(i)
            reads.append(copy.deepcopy(self.reads[i]))
        deletedNodes = set()
        # print("here1")
        nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,deletedNodes = self.constructNewGraph2(reads,deletedNodes)
        # print(weights)
        # print("here2")
        # print("here")
        
        # if outfile2 is not None:
        #     reads2 = []
        #     for i in range(len(self.reads)):
        #         reads2.append(copy.deepcopy(self.reads[i]))
        #     nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,deletedNodes = self.constructNewGraph2(reads,deletedNodes)
        #     self.htmlOutput(outfile2,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path = self.reads[0].corrected_path)
        #     ##Reset
        #     nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,deletedNodes = self.constructNewGraph2(reads,deletedNodes)
        # for key in sorted(weights.keys()):
        #     print(key,weights[key])
        ### Rewalk through the graph for each source node, constructing a priority queue of branches not walked
        ### (Could do this above, but it creates a bunch of branches that end up being deleted after correction)
        ### I think it's faster to rewalk the graph and build it after correction
        alt_paths = []
        alt_path_set = set()
        # print("here0",nodes[node_indexes[('b',5198)]].path)
        for i,fwd_node_idx in enumerate(source_nodes):
            bck_node_idx = node_indexes[('b',fwd_node_idx)]
            node = nodes[bck_node_idx]
            node.path = [fwd_node_idx]
            node.visited = True
            while len(edges[fwd_node_idx]) != 0:
                weight_list = []
                for j in edges[fwd_node_idx]:
                    weight_list.append(weights[(fwd_node_idx,j)])
                last_fwd_node_idx = fwd_node_idx
                fwd_node_idx = edges[fwd_node_idx][max(range(len(weight_list)),key=weight_list.__getitem__)]
                # if last_fwd_node_idx == 5172:
                #     print("WTF")
                #     print(edges[last_fwd_node_idx])
                for i,next_node in enumerate(edges[last_fwd_node_idx]):
                    # if (last_fwd_node_idx,next_node)  == (5172,5207):
                    # if (last_fwd_node_idx,next_node)  == (5172,5198):
                    #     print(next_node,nodes[node_indexes[('b',last_fwd_node_idx)]].path)
                        # print("5172 path",nodes[node_indexes[('b',last_fwd_node_idx)]].path)
                        # print("HERE2, why not including?")
                    nodes[node_indexes[('b',next_node)]].path = nodes[node_indexes[('b',last_fwd_node_idx)]].path + [next_node]
                    if next_node == fwd_node_idx:
                        alt_path_set.add((weights[(last_fwd_node_idx,next_node)],log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx))
                        continue
                    # if next_node == 5172:
                    if (weights[(last_fwd_node_idx,next_node)],log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx) not in alt_path_set:
                        pq.heappush(alt_paths,(-weights[(last_fwd_node_idx,next_node)],-log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx))
                        alt_path_set.add((weights[(last_fwd_node_idx,next_node)],log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx))
                bck_node_idx = node_indexes[('b',fwd_node_idx)]
        
        ### While the p. queue is not empty:
        ###  - If the node has been deleted (meaning it was not part of a delta > psi in length), ignore it.
        ###  - If the node has not been deleted:
        ###       i - build a path from -1 of the node to the next node previously seen in a greedy walk, correct reads based on this short path
        ###      ii - store new alternative paths to priority queue
        ###     iii - reconstruct the graph with the corrected paths, keeping track of deleted nodes.
        
        # print("here1",nodes[node_indexes[('b',5198)]].path)
        # for node in nodes:
        #     for i in range(len(node.path) - 1):
        #         assert node.path[i+1] > node.path[i]
        # ### So paths are valid up to here.
        
        iterations = 0
        # print(len(alt_paths))
        while alt_paths:
            iterations += 1
            # print(nodes[node_indexes['b',5217]].visited)
            test1,test2,fwd_node_idx,node1 = pq.heappop(alt_paths)
            # print(test1,fwd_node_idx,node1)
            if fwd_node_idx in deletedNodes or (node1, fwd_node_idx) not in weights:
                # print("here")
                # print(fwd_node_idx in deletedNodes)
                continue
            # print('iter',iterations)
            # if node1 == 5217:
            #     print("WHY?")
            #     print(fwd_node_idx)
            # else:
            #     print("Why2", fwd_node_idx, node1)
            # if (fwd_node_idx,node1) == (5198,5172):
            #     print("HERE, WHY NOT CORRECTING?")
            # path = [node1,fwd_node_idx]
            # path = copy.deepcopy(nodes[node_indexes['b',fwd_node_idx]].path)
            path = copy.deepcopy(nodes[node_indexes['b',node1]].path)[-psi:] + [fwd_node_idx]
            # print(path)
            # if node1 == 1784:
            #     print("huh...",path)
            # print(iterations)
            # print(len(deletedNodes))
            bck_node_idx = node_indexes[('b',fwd_node_idx)]
            node = nodes[bck_node_idx]
            while len(edges[fwd_node_idx]) != 0 and not nodes[bck_node_idx].visited:
                weight_list = []
                for j in edges[fwd_node_idx]:
                    weight_list.append(weights[(fwd_node_idx,j)])
                last_fwd_node_idx = fwd_node_idx
                fwd_node_idx = edges[fwd_node_idx][max(range(len(weight_list)),key=weight_list.__getitem__)]
                for i,next_node in enumerate(edges[last_fwd_node_idx]):
                    # if last_fwd_node_idx == 5217:
                    #     print("here")
                    if next_node == fwd_node_idx:
                        alt_path_set.add((weights[(last_fwd_node_idx,next_node)],log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx))
                        continue
                    # nodes[node_indexes[('b',next_node)]].path = path + [next_node]
                    if (weights[(last_fwd_node_idx,next_node)],log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx) not in alt_path_set:
                        pq.heappush(alt_paths,(-weights[(last_fwd_node_idx,next_node)],-log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx))
                        alt_path_set.add((weights[(last_fwd_node_idx,next_node)],log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx))
                # nodes[bck_node_idx].path = path
                nodes[bck_node_idx].visited = True
                bck_node_idx = node_indexes[('b',fwd_node_idx)]
                path.append(fwd_node_idx)
            # for i in sorted(list(alt_path_set)):
            #     print(i)
            if nodes[bck_node_idx].visited:
                for _ in range(psi):
                    if len(edges[fwd_node_idx]) == 0:
                        break
                    weight_list = []
                    for j in edges[fwd_node_idx]:
                        weight_list.append(weights[(fwd_node_idx,j)])
                    last_fwd_node_idx = fwd_node_idx
                    fwd_node_idx = edges[fwd_node_idx][max(range(len(weight_list)),key=weight_list.__getitem__)]
                    # for i,next_node in enumerate(edges[last_fwd_node_idx]):
                    #     if next_node == fwd_node_idx:
                    #         alt_path_set.add((weights[(last_fwd_node_idx,next_node)],log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx))
                    #         continue
                    #     nodes[node_indexes[('b',next_node)]].path = path + [next_node]
                    #     if (weights[(last_fwd_node_idx,next_node)],log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx) not in alt_path_set:
                    #         pq.heappush(alt_paths,(weights[(last_fwd_node_idx,next_node)],log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx))
                    #         alt_path_set.add((weights[(last_fwd_node_idx,next_node)],log_likelihoods[(last_fwd_node_idx,next_node)],next_node,last_fwd_node_idx))
                    # nodes[bck_node_idx].path = path
                    bck_node_idx = node_indexes[('b',fwd_node_idx)]
                    # nodes[bck_node_idx].visited = True
                    path.append(fwd_node_idx)
            # Correct reads based on this short path
            # print(path)
            # for idx in path:
            #     print("path",idx)
            reads = []
            for i in differed_reads:
                # _,self.reads[i].corrected_path = self.correctPaths(self.reads[i].corrected_path,path,psi=psi)
                # if path[0:9] == [5172,5200,5201,5202,5203,5204,5205,5206,5207]:
                # if path[0] == 5172:
                debug = False
                # if 5172 in path and 5207 in path and 5172 in self.reads[i].corrected_path and 5207 in self.reads[i].corrected_path and 5200 in self.reads[i].corrected_path:
                    # print("OMG")
                    # print(path)
                    # print(self.reads[i].corrected_path)
                    # debug = True
                # _,self.reads[i].corrected_path = self.correctPathsAndEnds2(self.reads[i].corrected_path,path,psi=psi,debug=debug)
                _,self.reads[i].corrected_path = self.correctPathsAndEnds3(self.reads[i].corrected_path,path,nodes,node_indexes,psi=psi,debug=debug)
                reads.append(copy.deepcopy(self.reads[i]))
            ### Reconstruct the graph with the corrected paths, keeping track of deleted nodes
            # print("here3")
            nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,deletedNodes = self.constructNewGraph2(reads,deletedNodes)
            # print("here4")
            
        # if outfile1 is not None:
        #     # print("Are we reaching here")
        #     reads = []
        #     for i in range(len(self.reads)):
        #         reads.append(copy.deepcopy(self.reads[i]))
        #     nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,deletedNodes = self.constructNewGraph2(reads,deletedNodes)
        #     self.htmlOutput(outfile1,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path = self.reads[0].corrected_path)
        ### Run approach similar to that used by Stringtie
        ### i.e. start at highest coverage node, and walk heaviest path based on the reads compatible with the walk up to that point (will require changing how somethings are stored)
        ###     - Need to include edges in both directions
        ###     - Need to store which reads contribute to each edge (convert weights in to a list of read ids?)
        ###     - Need to store how many reads are contributing to each node.
        ###
        ### Remove reads compatible with the full walk, repeat.
        # print(source_nodes)
        representative_paths = []
        remaining_reads = list(range(len(self.reads)))
        while remaining_reads:
            # print(remaining_reads)
            reads = []
            for i in remaining_reads:
                reads.append(copy.deepcopy(self.reads[i]))
            nodes,node_indexes,fwd_edges,rev_edges,weights,node_support_list = self.getStringtieDataStructures(reads)
            max_node_idx = max(range(len(node_support_list)),key=node_support_list.__getitem__)
            if max_node_idx == 0:
                start_idx = None
            else:
                start_idx = max_node_idx
            if max_node_idx == len(nodes) - 1:
                end_idx = None
            else:
                end_idx = max_node_idx
            fwd_max_node_idx = node_indexes[('f',max_node_idx)]
            # read_idxs = [[0,0] for _ in range(len(reads))]
            # read_idxs = []
            read_idxs = np.zeros((len(reads),2),dtype=np.int64)
            excluded_reads = set()
            all_reads = set()
            read_id_to_idx = {}
            read_idx_to_id = {}
            for i,read in enumerate(reads):
                read_id_to_idx[read.name] = i
                read_idx_to_id[i] = read.name
                all_reads.add(read.name)
                flag = True
                for j,node in enumerate(read.corrected_path): ### TODO: Replace with binary search later.
                    if node == fwd_max_node_idx:
                        read_idxs[i,:] = [j,j]
                        # read_idxs.append([j,j])
                        flag = False
                        break
                if flag:
                    read_idxs[i,:] = [-1,-1]
                    # read_idxs.append([None,None])
                    excluded_reads.add(read.name)
            # print(excluded_reads)
            representative_path = [fwd_max_node_idx]
            while start_idx is not None or end_idx is not None:
                if start_idx is None:
                    # break
                    ## move along end index
                    fwd_end_idx = node_indexes[('f',end_idx)]
                    
                    fwd_weights = []
                    for fwd_node in  fwd_edges[fwd_end_idx]:
                        weight = 0
                        for read_id in weights[(fwd_end_idx,fwd_node)]:
                            if read_id in excluded_reads:
                                continue
                            weight += 1
                        fwd_weights.append(weight)
                    
                    if max(fwd_weights) == 0:
                        break
                    next_end_idx = fwd_edges[fwd_end_idx][max(range(len(fwd_weights)),key=fwd_weights.__getitem__)]
                    next_node_id = next_end_idx
                    for i in range(len(reads)):
                        read_idx = read_idxs[i,1]
                        if read_idx == -1:
                            continue
                        if read_idx + 1 != len(reads[i].corrected_path):
                            if reads[i].corrected_path[read_idx + 1] == next_end_idx:
                                read_idxs[i,1] = read_idx + 1
                            else:
                                read_idxs[i,1] = -1
                                excluded_reads.add(reads[i].name)
                        else:
                            read_idxs[i,1] = -1
                    
                    end_idx = node_indexes[('b',next_end_idx)]
                    if len(fwd_edges[next_end_idx]) == 0:
                        end_idx = None
                    
                elif end_idx is None:
                    # break
                    ## move along start index
                    fwd_start_idx = node_indexes[('f',start_idx)]
                    rev_weights = []
                    for rev_node in rev_edges[fwd_start_idx]:
                        weight = 0
                        for read_id in weights[(rev_node,fwd_start_idx)]:
                            if read_id in excluded_reads:
                                continue
                            weight += 1
                        rev_weights.append(weight)
                    
                    if max(rev_weights) == 0:
                        break
                    next_start_idx = rev_edges[fwd_start_idx][max(range(len(rev_weights)),key=rev_weights.__getitem__)]
                    next_node_id = next_start_idx
                    for i in range(len(reads)):
                        read_idx = read_idxs[i,0]
                        if read_idx is -1:
                            continue
                        if read_idx != 0:
                            if reads[i].corrected_path[read_idx - 1] == next_start_idx:
                                read_idxs[i,0] = read_idx - 1
                            else:
                                read_idxs[i,0] = -1
                                excluded_reads.add(reads[i].name)
                        else:
                            read_idxs[i,0] = -1
                    
                    start_idx = node_indexes[('b',next_start_idx)]
                    if len(rev_edges[next_start_idx]) == 0:
                        start_idx = None
                else:
                    ## decide which to move along based on the max weighted path.
                    fwd_start_idx = node_indexes[('f',start_idx)]
                    fwd_end_idx = node_indexes[('f',end_idx)]

                    rev_weights = []
                    for rev_node in rev_edges[fwd_start_idx]:
                        weight = 0
                        for read_id in weights[(rev_node,fwd_start_idx)]:
                            if read_id in excluded_reads:
                                continue
                            weight += 1
                        rev_weights.append(weight)
                    fwd_weights = []
                    for fwd_node in  fwd_edges[fwd_end_idx]:
                        weight = 0
                        for read_id in weights[(fwd_end_idx,fwd_node)]:
                            if read_id in excluded_reads:
                                continue
                            weight += 1
                        fwd_weights.append(weight)
                    # print(len(rev_weights))
                    # print(len(fwd_weights))
                    if max(rev_weights) < max(fwd_weights):
                        #move along end index
                        if max(fwd_weights) == 0:
                            break
                        next_end_idx = fwd_edges[fwd_end_idx][max(range(len(fwd_weights)),key=fwd_weights.__getitem__)]
                        next_node_id = next_end_idx
                        for i in range(len(reads)):
                            read_idx = read_idxs[i,1]
                            if read_idx == -1:
                                continue
                            if read_idx + 1 != len(reads[i].corrected_path):
                                if reads[i].corrected_path[read_idx + 1] == next_end_idx:
                                    read_idxs[i,1] = read_idx + 1
                                else:
                                    read_idxs[i,1] = -1
                                    excluded_reads.add(reads[i].name)
                            else:
                                read_idxs[i,1] = -1
                        
                        end_idx = node_indexes[('b',next_end_idx)]
                        if len(fwd_edges[next_end_idx]) == 0:
                            end_idx = None
                    else:
                        #move along start index
                        if max(rev_weights) == 0:
                            break
                        next_start_idx = rev_edges[fwd_start_idx][max(range(len(rev_weights)),key=rev_weights.__getitem__)]
                        next_node_id = next_start_idx
                        for i in range(len(reads)):
                            read_idx = read_idxs[i,0]
                            if read_idx == -1:
                                continue
                            if read_idx != 0:
                                if reads[i].corrected_path[read_idx - 1] == next_start_idx:
                                    read_idxs[i,0] = read_idx - 1
                                else:
                                    read_idxs[i,0] = -1
                                    excluded_reads.add(reads[i].name)
                            else:
                                read_idxs[i,0] = -1
                        
                        start_idx = node_indexes[('b',next_start_idx)]
                        if len(rev_edges[next_start_idx]) == 0:
                            start_idx = None
                representative_path.append(next_node_id)
            accounted_for_reads = all_reads - excluded_reads
            # print(excluded_reads)
            # print(accounted_for_reads)
            representative_paths.append(sorted(representative_path))
            to_delete_reads = []
            for read_id in accounted_for_reads:
                to_delete_reads.append(read_id_to_idx[read_id])
            for idx in sorted(to_delete_reads,reverse=True):
                del remaining_reads[idx]
        return representative_paths
    
    def getStringtieDataStructures(self,reads):
        fwd_edges = {}
        rev_edges = {}
        weights = {}
        u_set = set()
        v_set = set()
        node_support_counts = {}
        for read in reads:
            for i in range(len(read.corrected_path) - 1):
                if (read.corrected_path[i],read.corrected_path[i+1]) in weights:
                    weights[(read.corrected_path[i],read.corrected_path[i+1])].append(read.name)
                else:
                    weights[(read.corrected_path[i],read.corrected_path[i+1])] = [read.name]
                    if read.corrected_path[i] in fwd_edges:
                        fwd_edges[read.corrected_path[i]].append(read.corrected_path[i+1])
                    else:
                        fwd_edges[read.corrected_path[i]] = [read.corrected_path[i+1]]
                    if read.corrected_path[i+1] in rev_edges:
                        rev_edges[read.corrected_path[i+1]].append(read.corrected_path[i])
                    else:
                        rev_edges[read.corrected_path[i+1]] = [read.corrected_path[i]]
                    u_set.add(read.corrected_path[i])
                    v_set.add(read.corrected_path[i+1])
                if read.corrected_path[i] in node_support_counts:
                    node_support_counts[read.corrected_path[i]] += 1
                else:
                    node_support_counts[read.corrected_path[i]] = 1
            if read.corrected_path[len(read.corrected_path) - 1] in node_support_counts:
                node_support_counts[read.corrected_path[len(read.corrected_path) - 1]] += 1
            else:
                node_support_counts[read.corrected_path[len(read.corrected_path) - 1]] = 1
        source_nodes = list( u_set - v_set )
        end_nodes = list( v_set - u_set )
        for end_node in  end_nodes:
            fwd_edges[end_node] = []
        for source_node in source_nodes:
            rev_edges[source_node] = []
        nodes = []
        node_indexes = sorted( list( u_set | v_set ) )
        node_indexes2 = {}
        node_support_list = []
        for i,j in enumerate(node_indexes):
            nodes.append(self.nodes[j])
            node_indexes2[('f',i)] = j
            node_indexes2[('b',j)] = i
            node_support_list.append(node_support_counts[j]) 
        # for node in nodes:
        #     node.indegree = 0
        # for j in edges:
        #     for k in edges[j]:
        #         nodes[node_indexes2['b',k]].indegree += 1
        return nodes,node_indexes2,fwd_edges,rev_edges,weights,node_support_list
    
    def getRepresentativePaths(self,outfile1=None,outfile2=None):
        if outfile1 is not None:
            read_paths = []
        for i in range(len(self.reads)):
            self.reads[i].corrected_path = self.reads[i].path.copy()
            if outfile1 is not None:
                read_paths.append(self.reads[i].path)
        # print(len(self.reads))
        if outfile1 is not None:
            nodes,node_indexes,source_nodes,edges,weights = self.constructGraphFromPaths(read_paths)
            del read_paths
            self.htmlOutput(outfile1,nodes,node_indexes,source_nodes,edges,weights)
            
        representative_paths = []
        differed_reads = list(range(len(self.reads)))
        index = 0
        flag = False
        while differed_reads:
            new_differed_reads = []
            # nodes,source_nodes,edges,weights = self.constructNewGraph()
            reads = []
            for i in differed_reads:
                reads.append(copy.deepcopy(self.reads[i]))
            nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructNewGraph(reads) ## This is changing every time because we're using the corrected paths.
            # if index != 0:
            #     assert (last_nodes,last_node_indexes,last_source_nodes,last_edges,last_weights) == (nodes,node_indexes,source_nodes,edges,weights)
            # last_nodes,last_node_indexes,last_source_nodes,last_edges,last_weights = nodes,node_indexes,source_nodes,edges,weights 
            ### - Something wrong with the logic above, NON-DETERMINISTIC!!!
            new_rep_path = self.getLongestPath2(nodes,node_indexes,source_nodes,edges,weights)
            # new_rep_path = self.getMostLikelyPath(nodes,node_indexes,source_nodes,end_nodes,edges,log_likelihoods)
            # representative_paths.append(new_rep_path)
            if flag:
            # if index != 0 and new_rep_path == representative_paths[index - 1]:
                #Do what? - For now select a random read's path and use that as new possibility (need to fix this in future with better idea)
                representative_paths.append(self.reads[differed_reads[0]].corrected_path)
                # print("here!")
            else:
                representative_paths.append(new_rep_path)
            # print(representative_paths[index])
            # print(differed_reads)
            for i in differed_reads:
                differed,self.reads[i].corrected_path = self.correctPaths(self.reads[i].corrected_path,representative_paths[index],psi= 10)
                print(differed)
                if differed:
                    new_differed_reads.append(i)
            if differed_reads == new_differed_reads:
                flag = True
            else:
                differed_reads = new_differed_reads
                flag = False
            index += 1
            print(differed_reads)
            print(index)
        if outfile2 is not None:
            nodes,node_indexes,source_nodes,edges,weights = self.constructGraphFromPaths(representative_paths)
            self.htmlOutput(outfile2,nodes,node_indexes,source_nodes,edges,weights)
        return representative_paths
    
    def constructNewGraph(self,reads, n = 1):
        edges = {}
        weights = {}
        u_set = set()
        v_set = set()
        # print(len(reads))
        for read in reads:
            for i in range(len(read.corrected_path) - 1):
                # if len(reads) == 1:
                #     print(i)
                if (read.corrected_path[i],read.corrected_path[i+1]) in weights:
                    weights[(read.corrected_path[i],read.corrected_path[i+1])] += 1
                else:
                    weights[(read.corrected_path[i],read.corrected_path[i+1])] = 1
                    if read.corrected_path[i] in edges:
                        # if len(reads) == 1:
                        #     print("WTF")
                        edges[read.corrected_path[i]].append(read.corrected_path[i+1])
                    else:
                        edges[read.corrected_path[i]] = [read.corrected_path[i+1]]
                    u_set.add(read.corrected_path[i])
                    v_set.add(read.corrected_path[i+1])
            # for i in range(len(read.corrected_path)):
            #     print(read.corrected_path[i])
        source_nodes = list( u_set - v_set )
        end_nodes = list(v_set - u_set)
        for end_node in  end_nodes:
            edges[end_node] = []
        nodes = []
        node_indexes = sorted(list(u_set | v_set))
        node_indexes2 = {}
        # print(node_indexes)
        log_likelihoods = {}
        for i,j in enumerate(node_indexes):
            nodes.append(self.nodes[j])
            node_indexes2[('f',i)] = j
            node_indexes2[('b',j)] = i
            if j in edges:
                total_read_count = 0
                for k in edges[j]:
                    # total_read_count += weights[(j,k)]
                    total_read_count += weights[(j,k)] ** n 
                for k in edges[j]:
                    # log_likelihoods[(j,k)] = - np.log(weights[(j,k)] / total_read_count)
                    log_likelihoods[(j,k)] = - np.log(weights[(j,k)] ** n / total_read_count)
        for node in nodes:
            node.indegree = 0
        for j in edges:
            for k in edges[j]:
                nodes[node_indexes2['b',k]].indegree += 1
        # for j in sorted(edges.keys()):
        #     for k in sorted(edges[j]):
        #         print(j,k,weights[(j,k)])
        # print(source_nodes)
        # print(end_nodes)
        # if len(reads) == 1:
        #     print(edges)
        return nodes,node_indexes2,source_nodes,end_nodes,edges,weights,log_likelihoods
    
    def constructNewGraph2(self,reads,deletedNodes, n = 1):
        edges = {}
        weights = {}
        u_set = set()
        v_set = set()
        # print(len(reads))
        for read in reads:
            for i in range(len(read.corrected_path) - 1):
                # if len(reads) == 1:
                #     print(i)
                if (read.corrected_path[i],read.corrected_path[i+1]) in weights:
                    weights[(read.corrected_path[i],read.corrected_path[i+1])] += 1
                else:
                    weights[(read.corrected_path[i],read.corrected_path[i+1])] = 1
                    if read.corrected_path[i] in edges:
                        # if len(reads) == 1:
                        #     print("WTF")
                        if (read.corrected_path[i],read.corrected_path[i+1]) == (5172, 5207):
                            print("Hmm")
                        edges[read.corrected_path[i]].append(read.corrected_path[i+1])
                    else:
                        edges[read.corrected_path[i]] = [read.corrected_path[i+1]]
                    u_set.add(read.corrected_path[i])
                    v_set.add(read.corrected_path[i+1])
            # for i in range(len(read.corrected_path)):
            #     print(read.corrected_path[i])

        source_nodes = list( u_set - v_set )
        end_nodes = list(v_set - u_set)
        for end_node in  end_nodes:
            edges[end_node] = []
        nodes = []
        node_indexes = sorted(list(u_set | v_set))
        node_indexes2 = {}
        # print(node_indexes)
        log_likelihoods = {}
        for i in range(len(self.nodes)):
            if i not in u_set and i not in v_set:
                deletedNodes.add(i)
        # print(len(deletedNodes))
        for i,j in enumerate(node_indexes):
            nodes.append(self.nodes[j])
            node_indexes2[('f',i)] = j
            node_indexes2[('b',j)] = i
            if j in edges:
                total_read_count = 0
                for k in edges[j]:
                    # total_read_count += weights[(j,k)]
                    total_read_count += weights[(j,k)] ** n 
                for k in edges[j]:
                    # log_likelihoods[(j,k)] = - np.log(weights[(j,k)] / total_read_count)
                    log_likelihoods[(j,k)] = - np.log(weights[(j,k)] ** n / total_read_count)
        for node in nodes:
            node.indegree = 0
        for j in edges:
            for k in edges[j]:
                nodes[node_indexes2['b',k]].indegree += 1
        # for key in sorted(node_indexes2.keys()):
        #     print("idx:",key,node_indexes2[key])
        # for j in sorted(edges.keys()):
        #     for k in sorted(edges[j]):
        #         print(j,k,weights[(j,k)])
        # print(sorted(list(deletedNodes)))
        # print(source_nodes)
        # print(end_nodes)
        # if len(reads) == 1:
        #     print(edges)
        return nodes,node_indexes2,source_nodes,end_nodes,edges,weights,log_likelihoods,deletedNodes
    
    def constructGraphFromPaths(self,paths,n = 1):
        edges = {}
        weights = {}
        u_set = set()
        v_set = set()
        for path in paths:
            for i in range(len(path) - 1):
                if (path[i],path[i+1]) in weights:
                    weights[(path[i],path[i+1])] += 1
                else:
                    weights[(path[i],path[i+1])] = 1
                    if path[i] in edges:
                        edges[path[i]].append(path[i+1])
                    else:
                        edges[path[i]] = [path[i+1]]
                    u_set.add(path[i])
                    v_set.add(path[i+1])

        # source_nodes = list( u_set - v_set )
        source_nodes = []
        end_nodes = list(v_set - u_set)
        for end_node in  end_nodes:
            edges[end_node] = []
        nodes = []
        node_indexes = sorted(list(u_set | v_set))
        node_indexes2 = {}
        # print(node_indexes)
        log_likelihoods = {}
        for i,j in enumerate(node_indexes):
            nodes.append(self.nodes[j])
            if j in self.source_node_indices:
                source_nodes.append(j)
            node_indexes2[('f',i)] = j
            node_indexes2[('b',j)] = i
            if j in edges:
                total_read_count = 0
                for k in edges[j]:
                    # total_read_count += weights[(j,k)]
                    total_read_count += weights[(j,k)] ** n
                for k in edges[j]:
                    # log_likelihoods[(j,k)] = - np.log(weights[(j,k)] / total_read_count)
                    log_likelihoods[(j,k)] = - np.log(weights[(j,k)] ** n / total_read_count)
        # print(source_nodes)
        end_nodes = []
        for i,node in enumerate(nodes):
            node.indegree = 0
            if node.end_node_flag:
                end_nodes.append(node_indexes2[('f',i)])
        for j in edges:
            for k in edges[j]:
                nodes[node_indexes2['b',k]].indegree += 1
        return nodes,node_indexes2,source_nodes,end_nodes,edges,weights,log_likelihoods
    
    def s(self,x,y,a,b):
        if x == y:
            return a
        else:
            return -b
    
    def alignPaths2(self,x,y):
        #x - first seq (query for the purposes of defining insertion / deletion)
        #y - second seq (reference for the purposes of defining insertion / deletion)
        ### Don't need to do full pairwise alignment, as paths are indices of topo sorted graph. Can use that fact to just compare indices instead
        alignment = []
        i = 0
        j = 0
        while True:
            if x[i] == y[j]:
                alignment.append(0)
                i += 1
                j += 1
            elif x[i] < y[j]:
                alignment.append(1)
                i += 1
            elif x[i] > y[j]:
                alignment.append(2)
                j += 1
            if i == len(x) or j == len(y):
                break
        for _ in range(i,len(x)):
            alignment.append(1)
        for _ in range(j,len(y)):
            alignment.append(2)
        return alignment
    
    def scorePaths(self, x, y):
        #x - first seq (query for the purposes of defining insertion / deletion)
        #y - second seq (reference for the purposes of defining insertion / deletion)
        ### Don't need to do full pairwise alignment, as paths are indices of topo sorted graph. Can use that fact to just compare indices instead
        score = 0
        i = 0
        j = 0
        while True:
            if x[i] == y[j]:
                score += 1
                i += 1
                j += 1
            elif x[i] < y[j]:
                i += 1
                score -= 1
            elif x[i] > y[j]:
                j += 1
                score -= 1
            if i == len(x) or j == len(y):
                break
        return score
    
    def collapseLinearStretches(self,nodes,node_indexes,edges,weights,log_likelihoods):
        to_delete = set()
        new_end_nodes = []
        for idx in sorted(node_indexes.keys()):
            print("idx,",idx,node_indexes[idx])
        for i,node in enumerate(nodes):
            node.visited = False
            if node.end_node_flag:
                print("flag",i)
        for i,node in enumerate(nodes):
            if node.visited:
                continue
            node_idx = node_indexes['f',i]
            # print(node_idx)
            # if node_idx == 4949:
            #     print("HERE1!")
            if len(edges[node_idx]) == 1 and not node.end_node_flag:
                # if node_idx == 4949:
                #     print("HERE2!")
                next_node_idx = edges[node_idx][0]
                back_next_node_idx = node_indexes['b',next_node_idx]
                while len(edges[node_idx]) == 1 and nodes[back_next_node_idx].indegree == 1 and not nodes[back_next_node_idx].start_node_flag and not nodes[back_next_node_idx].end_node_flag:
                    node.char += nodes[back_next_node_idx].char
                    nodes[back_next_node_idx].visited = True
                    to_delete.add(back_next_node_idx)
                    edges[node_idx] = copy.deepcopy(edges[next_node_idx])
                    del edges[next_node_idx]
                    del log_likelihoods[(node_idx,next_node_idx)]
                    del weights[(node_idx,next_node_idx)]
                    if len(edges[node_idx]) == 1:
                        k = edges[node_idx][0]
                        log_likelihoods[(node_idx,k)] = log_likelihoods[(next_node_idx,k)]
                        del log_likelihoods[(next_node_idx,k)]
                        weights[(node_idx,k)] = weights[(next_node_idx,k)]
                        del weights[(next_node_idx,k)]
                        next_node_idx = edges[node_idx][0]
                        back_next_node_idx = node_indexes['b',next_node_idx]
                    elif len(edges[node_idx]) > 1:
                        for k in edges[node_idx]:
                            log_likelihoods[(node_idx,k)] = log_likelihoods[(next_node_idx,k)]
                            del log_likelihoods[(next_node_idx,k)]
                            weights[(node_idx,k)] = weights[(next_node_idx,k)]
                            del weights[(next_node_idx,k)]
                        break
                    # else:
                    #     new_end_nodes.append(node_idx)
                    #     break
                ### New code
                # if node_idx == 5429:
                #     print("HERE4!")
                #     print(back_next_node_idx,next_node_idx,nodes[back_next_node_idx].end_node_flag)
                #     print(len(edges[node_idx]))
                #     print(nodes[back_next_node_idx].indegree)
                if nodes[back_next_node_idx].end_node_flag:
                    if len(edges[node_idx]) == 1 and nodes[back_next_node_idx].indegree == 1 and not nodes[back_next_node_idx].start_node_flag:
                        node.char += nodes[back_next_node_idx].char
                        nodes[back_next_node_idx].visited = True
                        to_delete.add(back_next_node_idx)
                        edges[node_idx] = copy.deepcopy(edges[next_node_idx])
                        del edges[next_node_idx]
                        del log_likelihoods[(node_idx,next_node_idx)]
                        del weights[(node_idx,next_node_idx)]
                        if len(edges[node_idx]) == 1:
                            k = edges[node_idx][0]
                            log_likelihoods[(node_idx,k)] = log_likelihoods[(next_node_idx,k)]
                            del log_likelihoods[(next_node_idx,k)]
                            weights[(node_idx,k)] = weights[(next_node_idx,k)]
                            del weights[(next_node_idx,k)]
                            next_node_idx = edges[node_idx][0]
                            back_next_node_idx = node_indexes['b',next_node_idx]
                        elif len(edges[node_idx]) > 1:
                            for k in edges[node_idx]:
                                log_likelihoods[(node_idx,k)] = log_likelihoods[(next_node_idx,k)]
                                del log_likelihoods[(next_node_idx,k)]
                                weights[(node_idx,k)] = weights[(next_node_idx,k)]
                                del weights[(next_node_idx,k)]
                        new_end_nodes.append(node_idx)
                        print("here",node_idx)
                    elif len(edges[node_idx]) != 0:
                        new_end_nodes.append(node_idx)
                        # print("here3",node_idx)
                    else:
                        new_end_nodes.append(next_node_idx)
                        # print("here2",next_node_idx)
            else:
                if nodes[i].end_node_flag:
                    new_end_nodes.append(node_idx)
                ### End new code
                
        new_nodes = []
        new_node_indexes = {}
        counter = 0
        for i,node in enumerate(nodes):
            # if node_indexes[('f',i)] == 457:
            #     print('deleting', i in to_delete)
            if i in to_delete:
                continue
            new_nodes.append(node)
            new_node_indexes[('f',counter)] = node_indexes[('f',i)]
            new_node_indexes[('b',node_indexes[('f',i)])] = counter
            counter += 1
        return new_nodes,new_node_indexes,edges,weights,log_likelihoods,new_end_nodes
    
    def popBubbles(self,nodes,node_indexes,edges,weights,log_likelihoods):
        #Placeholder for now
        # print("here")
        to_collapse = []
        for node in nodes:
            node.visited = False
        for i, node in enumerate(nodes):
            if node.visited:
                continue
            # print("here")
            if node.align_ring_partner is not None:
                # print("here")
                node.visited = True
                node = node_indexes[('f',i)]
                next_node = i
                align_ring = [node]
                while nodes[next_node].align_ring_partner != node:
                    align_ring.append(nodes[next_node].align_ring_partner)
                    next_node = node_indexes[('b',nodes[next_node].align_ring_partner)]
                    nodes[next_node].visited = True
                # print(align_ring)
                to_collapse.append(align_ring)
        self.collapse(to_collapse,nodes,node_indexes,edges,weights,log_likelihoods)
        return None
    
    def collapse(self,align_rings,nodes,node_indexes,edges,weights,log_likelihoods):
        for align_ring in align_rings:
            pass
        return None
        # return new_nodes,node_indexes,edges,weights,log_likelihoods
    
    def compareAlignments(self,a1,a2):
        if len(a1) != len(a2):
            print("NOT EQUAL")
        match_counts1 = 0
        match_counts2 = 0
        continuous_ins1 = 0
        continuous_ins2 = 0
        continuous_del1 = 0
        continuous_del2 = 0
        continuous_ins1_list = []
        continuous_ins2_list = []
        continuous_del1_list = []
        continuous_del2_list = []
        for i in range(len(a1)):
            if a1[i] == 0:
                match_counts1 += 1
                continuous_ins1_list.append(continuous_ins1)
                continuous_del1_list.append(continuous_del1)
                continuous_ins1 = 0
                continuous_del1 = 0
            if a2[i] == 0:
                match_counts2 += 1
                continuous_ins2_list.append(continuous_ins2)
                continuous_del2_list.append(continuous_del2)
                continuous_ins2 = 0
                continuous_del2 = 0
        assert continuous_ins1_list == continuous_ins2_list
        assert continuous_del1_list == continuous_del2_list
        print(match_counts1)
        print(match_counts2)
    
    def alignPaths(self,x,y,a=1,b=100,o=0,e=1,anchor="both"):
        #x - first seq (query for the purposes of defining insertion / deletion)
        #y - second seq (reference for the purposes of defining insertion / deletion)
        #a - match score
        #b - mismatch penalty
        #o - gap open penalty
        #e - gap extend penalty
        inf = np.inf
        # x = x.upper()
        # y = y.upper()
        v = np.zeros([len(x) + 1, len(y) + 1])
        m = np.zeros([len(x) + 1, len(y) + 1])
        x_gap = np.zeros([len(x) + 1, len(y) + 1])
        y_gap = np.zeros([len(x) + 1, len(y) + 1])
        walkback_matrix = np.zeros([len(x) + 1, len(y) + 1])
    
        ### Initialize matrix values
        if anchor is None: #Semi-global alignment, allow gaps at the begining and end of y
            ## Both end gaps are free. (for y only)
            for i in range(1,len(x)+1): 
                v[i,0] = -o - i*e
                #v[i,0] = 0. #Redundant, but do it for now
                #x_gap[i,0] = -o - i*e
                y_gap[i,0] = -o - i*e
                #y_gap[i,0] = -inf
                #x_gap[i,0] = -inf
            for j in range(1,len(y)+1): 
                #v[0,j] = -o - j*e
                x_gap[0,j] = -o - j*e
                #y_gap[0,j] = -o - j*e
                v[0,j] = 0. #Redundant, but do it for now
                #x_gap[0,j] = -inf
                #y_gap[0,j] = -inf
        elif anchor == "both": ## Working as intended
            ## Both end gaps are penalized
            for i in range(1,len(x)+1): 
                v[i,0] = -o - i*e
                x_gap[i,0] = -o - i*e
                y_gap[i,0] = -inf
            for j in range(1,len(y)+1):
                v[0,j] = -o - j*e
                y_gap[0,j] = -o - j*e
                x_gap[0,j] = -inf
        elif anchor == "left" or anchor == "l" or anchor == "L":
            ## Gaps on the left side are penalized
            for i in range(1,len(x)+1): 
                v[i,0] = -o - i*e
                #v[i,0] = 0. #Redundant, but do it for now
                x_gap[i,0] = -o - i*e
                #y_gap[i,0] = -o - i*e
                y_gap[i,0] = -inf
                #x_gap[i,0] = -inf
            for j in range(1,len(y)+1): 
                v[0,j] = -o - j*e
                #x_gap[0,j] = -o - j*e
                y_gap[0,j] = -o - j*e
                #v[0,j] = 0. #Redundant, but do it for now
                x_gap[0,j] = -inf
                #y_gap[0,j] = -inf
        elif anchor == "right" or anchor == "r" or anchor == "R":
            ## Gaps on the right side are penalized
            for i in range(1,len(x)+1): 
                v[i,0] = -o - i*e
                #v[i,0] = 0. #Redundant, but do it for now
                x_gap[i,0] = -o - i*e
                #y_gap[i,0] = -o - i*e
                y_gap[i,0] = -inf
                #x_gap[i,0] = -inf
            for j in range(1,len(y)+1): 
                #v[0,j] = -o - j*e
                #x_gap[0,j] = -o - j*e
                y_gap[0,j] = -o - j*e
                v[0,j] = 0. #Redundant, but do it for now
                x_gap[0,j] = -inf
                #y_gap[0,j] = -inf
        ## Fill in the score matrices:
        for i in range(1,len(x)+1):
            for j in range(1,len(y)+1):
                m[i,j] = v[i-1,j-1] + self.s(x[i-1],y[j-1],a,b)
                x_gap[i,j] = max( x_gap[i-1,j], v[i-1,j] - o) - e
                y_gap[i,j] = max( y_gap[i,j-1], v[i,j-1] - o) - e
                v[i,j] = max([m[i,j],x_gap[i,j],y_gap[i,j]])
                possible_origins = 0
                if v[i,j] == m[i,j]:
                    possible_origins += 1
                if v[i,j] == x_gap[i,j]:
                    possible_origins += 2
                if v[i,j] == y_gap[i,j]:
                    possible_origins += 4
                walkback_matrix[i,j] = possible_origins
        ## Walkback
        if anchor is None or anchor == "left":
            i = len(x)
            #print i
            optimal_score = max(v[i,:])
            #print optimal_score
            for k in range(len(y)+1):
                if v[i,k] == optimal_score:
                    j = k
                    break
            #print "j"
            #print j
            walkback = [ 2 for x in range(len(y)- k)] 
        elif anchor == "both" or anchor == "right":
            i = len(x)
            j = len(y)
            walkback = []
    
        while i != 0 and j != 0:
            direction = walkback_matrix[i,j]
            if direction == 1:
                i = i-1
                j = j-1
                walkback.append(0)
            elif direction == 2:
                i= i-1
                walkback.append(1)
            elif direction == 4:
                j=j-1
                walkback.append(2)
            elif direction == 3: #Ambiguous alignment
                i = i-1 #Arbitrarily decide to use match
                j = j-1
                walkback.append(0)
            elif direction == 5: #Ambiguous alignment
                i = i-1 #Arbitrarily decide to use match
                j = j-1
                walkback.append(0)
            elif direction == 6: #Ambiguous alignment
                j=j-1 #Arbitrarily decide to use deletion
                walkback.append(2)
            elif direction == 7: #Ambiguous alignment
                i = i-1 #Arbitrarily decide to use match
                j = j-1
                walkback.append(0)
        if i > 0:
            walkback =  walkback + [1 for x in range(i)] # 1 is insertion
        elif j > 0:
            walkback =  walkback + [2 for x in range(j)] # 2 is deletion
        walkback.reverse()
        # print(walkback)
        return walkback
        
    def detectDifferentPathsAndEnds(self,q_path,r_path,psi=10):
        # alignment1 = self.alignPaths(q_path,r_path,a,b,o,e,anchor)
        # alignment1 = self.alignPaths(q_path,r_path)
        # og_q_path = copy.deepcopy(q_path)
        alignment = self.alignPaths2(q_path,r_path)
        # print("comparison",self.compareAlignments(alignment1,alignment))
        q_idx = 0
        r_idx = 0
        continuous_ins = 0
        continuous_del = 0
        diff_flag = False
        q_ops = []
        for op in alignment:
            if op == 0:
                if continuous_ins >= psi:
                    diff_flag = True
                    break
                # elif continuous_ins > 0:
                #     #Correct short insertion\
                #     for i in range(q_idx - continuous_ins,q_idx):
                #         q_ops.append(("ins",i,0))
                if continuous_del >= psi:
                    diff_flag = True
                    break
                # elif continuous_del > 0:
                #     #Correct short deletion
                #     for i in range(r_idx - continuous_del,r_idx):
                #         q_ops.append(("del",q_idx-continuous_ins,i))
                q_idx += 1
                r_idx += 1
                continuous_ins = 0
                continuous_del = 0
            elif op == 1:
                continuous_ins += 1
                q_idx += 1
                
            elif op == 2:
                continuous_del += 1
                r_idx += 1
                
        if continuous_ins >= psi:
            diff_flag = True
        # elif continuous_ins > 0:
        #     #Correct short insertion
        #     for i in range(q_idx - continuous_ins,q_idx):
        #         q_ops.append(("ins",i,0))
        if continuous_del >= psi:
            diff_flag = True
        # elif continuous_del > 0:
        #     #Correct short deletion
        #     for i in range(r_idx - continuous_del,r_idx):
        #         q_ops.append(("del",q_idx-continuous_ins,i))
        assert q_idx == len(q_path)
        assert r_idx == len(r_path)
        # print(list(reversed(q_ops)))
        # last_idx = np.inf
        # last_op = None
        # for op,q_idx,r_idx in reversed(sorted(q_ops,key=lambda t: (t[1],t[0]))):
        #     # print("here")
        # # for op,q_idx,r_idx in reversed(q_opss):
        #     # if last_idx <= q_idx:
        #     #     print("PANIC")
        #     #     print(last_idx,q_idx)
        #     #     print(last_op,op)
        #     if op == "ins":
        #         assert q_path[q_idx-1] < q_path[q_idx+1]
        #         del q_path[q_idx]
        #     elif op == "del":
        #         # try
        #         assert q_path[q_idx-1] < r_path[r_idx] < q_path[q_idx]
        #         # except:
        #         #     print("panic!")
        #         #     print(q_path[q_idx-5])
        #         #     print(q_path[q_idx-4])
        #         #     print(q_path[q_idx-3])
        #         #     print(q_path[q_idx-2])
        #         #     print(q_path[q_idx-1])
        #         #     print(q_path[q_idx])
        #         #     # print(q_path[q_idx+1])
        #         #     print(r_path[r_idx])
        #         q_path.insert(q_idx,r_path[r_idx])
            # last_idx = q_idx
            # last_op = op
                
        # alignment = self.alignPaths(q_path,r_path,a,b,o,e,anchor)
        # for i in range(len(q_path)-1):
        #     try:
        #         assert(q_path[i] < q_path[i+1])
        #     except:
        #         print(q_path[i])
        #         print(q_path[i+1])
        # if og_q_path == q_path:
        #     print("NOT CORRECTING 2")
        # for i in range(len(q_path)-1):
        #     assert q_path[i] < q_path[i+1]
        return diff_flag
        
    def correctPaths(self,q_path,r_path,psi=10):
        # alignment1 = self.alignPaths(q_path,r_path,a,b,o,e,anchor)
        # alignment1 = self.alignPaths(q_path,r_path)
        # og_q_path = copy.deepcopy(q_path)
        alignment = self.alignPaths2(q_path,r_path)
        # print("comparison",self.compareAlignments(alignment1,alignment))
        q_idx = 0
        r_idx = 0
        continuous_ins = 0
        continuous_del = 0
        diff_flag = False
        q_ops = []
        start_flag = True
        ignore_deletion_flag = False
        for op in alignment:
            if op == 0:
                # if not start_flag:
                #     last_q_match_idx = q_idx - continuous_ins - 1
                #     last_r_match_idx = r_idx - continuous_del - 1
                #     assert q_path[last_q_match_idx] == r_path[last_r_match_idx]
                # print(r_path[r_idx],q_path[q_idx],sep='\t')
                if continuous_ins >= psi:
                    diff_flag = True
                    ignore_deletion_flag = True
                elif continuous_ins > 0:
                    #Correct short insertion\
                    for i in range(q_idx - continuous_ins,q_idx):
                        q_ops.append(("ins",i,0))
                if continuous_del >= psi:
                    diff_flag = True
                elif continuous_del > 0 and not ignore_deletion_flag:
                    #Correct short deletion
                    for i in range(r_idx - continuous_del,r_idx):
                        q_ops.append(("del",q_idx-continuous_ins,i))
                ignore_deletion_flag = False
                start_flag = False
                q_idx += 1
                r_idx += 1
                continuous_ins = 0
                continuous_del = 0
            elif op == 1:
                if not start_flag:
                    continuous_ins += 1
                    # print("----",q_path[q_idx],sep = '\t')
                q_idx += 1
                
            elif op == 2:
                if not start_flag:
                    continuous_del += 1
                    # print(r_path[r_idx],"----",sep = '\t')
                r_idx += 1
                
        
        assert q_idx == len(q_path)
        assert r_idx == len(r_path)
        # print(list(reversed(q_ops)))
        last_q_idx = np.inf
        last_r_idx = np.inf
        last_op = None
        for op,q_idx,r_idx in reversed(sorted(q_ops,key=lambda t: (t[1],t[0],t[2]))):
            # print("here")
            # if last_q_idx <= q_idx:
            #     print("PANIC")
            #     print(last_q_idx,q_idx)
            #     print(last_r_idx,r_idx)
            #     print(last_op,op)
            if op == "ins":
                assert q_path[q_idx-1] < q_path[q_idx+1]
                del q_path[q_idx]
            elif op == "del":
                assert q_path[q_idx-1] < r_path[r_idx] < q_path[q_idx]
                q_path.insert(q_idx,r_path[r_idx])
            last_q_idx = q_idx
            last_r_idx = r_idx
            last_op = op
                
        # alignment = self.alignPaths(q_path,r_path,a,b,o,e,anchor)
        # for i in range(len(q_path)-1):
        #     try:
        #         assert(q_path[i] < q_path[i+1])
        #     except:
        #         print(q_path[i])
        #         print(q_path[i+1])
        # if og_q_path == q_path:
        #     print("NOT CORRECTING 2")
        for i in range(len(q_path)-1):
            assert q_path[i] < q_path[i+1]
        return (diff_flag,q_path)
    
    def correctPathsAndEnds3(self,q_path,r_path,nodes,node_indexes,psi=10,debug=False):
        # alignment1 = self.alignPaths(q_path,r_path,a,b,o,e,anchor)
        # alignment1 = self.alignPaths(q_path,r_path)
        # og_q_path = copy.deepcopy(q_path)
        for i in range(len(q_path)-1):
            assert q_path[i] < q_path[i+1]
        for i in range(len(r_path)-1):
            try:
                assert r_path[i] < r_path[i+1]
            except:
                print(r_path)
                raise
        # if r_path[-1] == 2021 and q_path[-1] == 5181:
        #     debug=True
        alignment = self.alignPaths2(q_path,r_path)
        # print("comparison",self.compareAlignments(alignment1,alignment))
        q_idx = 0
        r_idx = 0
        continuous_ins = 0
        continuous_del = 0
        diff_flag = False
        q_ops = []
        start_flag = True
        ignore_deletion_flag = False
        for op in alignment:
            if op == 0:
                # if not start_flag: 
                #     last_q_match_idx = q_idx - continuous_ins - 1
                #     last_r_match_idx = r_idx - continuous_del - 1
                #     assert q_path[last_q_match_idx] == r_path[last_r_match_idx]
                if debug:
                    print(r_path[r_idx],q_path[q_idx],sep='\t')
                if continuous_ins >= psi:
                    diff_flag = True
                    ignore_deletion_flag = True
                elif continuous_ins > 0:
                    #Correct short insertion\
                    # if debug:
                    #     print(continuous_ins)
                    for i in range(q_idx - continuous_ins,q_idx):
                        q_ops.append(("ins",i,0))
                    # if start_flag:
                    #     nodes[node_indexes[('b',q_path[q_idx])]].start_node_flag = True
                if continuous_del >= psi:
                    diff_flag = True
                elif continuous_del > 0 and not ignore_deletion_flag:
                    #Correct short deletion
                    # if debug:
                    #     print(continuous_del)
                    for i in range(r_idx - continuous_del,r_idx):
                        q_ops.append(("del",q_idx-continuous_ins,i))
                ignore_deletion_flag = False
                start_flag = False
                q_idx += 1
                r_idx += 1
                continuous_ins = 0
                continuous_del = 0
            elif op == 1:
                # if not start_flag:
                #     continuous_ins += 1
                continuous_ins += 1
                if debug:
                    print("----",q_path[q_idx],sep = '\t')
                q_idx += 1
                
            elif op == 2:
                if not start_flag:
                    continuous_del += 1
                if debug:
                    print(r_path[r_idx],"----",sep = '\t')
                r_idx += 1
                

        if continuous_ins >= psi:
            diff_flag = True
            ignore_deletion_flag = True
        elif continuous_ins > 0:
            #Correct short insertion\
            if debug:
                print("Hitting here for some reason")
            for i in range(q_idx - continuous_ins,q_idx):
                q_ops.append(("ins",i,0))
            # nodes[node_indexes[('b',q_path[q_idx])]].start_node_flag = True
        assert q_idx == len(q_path)
        assert r_idx == len(r_path)
        # print(list(reversed(q_ops)))
        last_q_idx = np.inf
        last_r_idx = np.inf
        last_op = None
        # if debug:
        #     print("q_ops",q_ops)
        for op,q_idx,r_idx in reversed(sorted(q_ops,key=lambda t: (t[1],t[0],t[2]))):
            # print("here")
            # if last_q_idx <= q_idx:
            #     print("PANIC")
            #     print(last_q_idx,q_idx)
            #     print(last_r_idx,r_idx)
            #     print(last_op,op)
            if op == "ins":
                if q_idx != len(q_path) - 1 and q_idx != 0:
                    # try:
                    assert q_path[q_idx-1] < q_path[q_idx+1]
                    # except:
                    #     print(q_path[q_idx-1])
                    #     print(q_path[q_idx])
                    #     print(q_path[q_idx+1])
                    #     raise
                del q_path[q_idx]
            elif op == "del":
                if q_idx != 0:
                    try:
                        assert q_path[q_idx-1] < r_path[r_idx] < q_path[q_idx]
                    except:
                        # if debug:
                        print(q_path[-1])
                        print(r_path[-1])
                        print("q_ops",q_ops)
                        print(q_path[q_idx-1], r_path[r_idx], q_path[q_idx])
                        raise
                q_path.insert(q_idx,r_path[r_idx])
            last_q_idx = q_idx
            last_r_idx = r_idx
            last_op = op
                
        # alignment = self.alignPaths(q_path,r_path,a,b,o,e,anchor)
        for i in range(len(q_path)-1):
            try:
                assert(q_path[i] < q_path[i+1])
            except:
                print(q_path[i])
                print(q_path[i+1])
                raise
        # if og_q_path == q_path:
        #     print("NOT CORRECTING 2")
        for i in range(len(q_path)-1):
            assert q_path[i] < q_path[i+1]
        # nodes[node_indexes[('b',q_path[0])]].start_node_flag = True
        self.source_node_indices.append(q_path[0])
        # print("end_nodes", node_indexes[('b',q_path[-1])])
        nodes[node_indexes[('b',q_path[-1])]].end_node_flag = True
        self.end_node_indices.append(q_path[-1])
        return (diff_flag,q_path)
    
    def correctPathsAndEnds2(self,q_path,r_path,psi=10,debug=False):
        # alignment1 = self.alignPaths(q_path,r_path,a,b,o,e,anchor)
        # alignment1 = self.alignPaths(q_path,r_path)
        # og_q_path = copy.deepcopy(q_path)
        for i in range(len(q_path)-1):
            assert q_path[i] < q_path[i+1]
        for i in range(len(r_path)-1):
            try:
                assert r_path[i] < r_path[i+1]
            except:
                print(r_path)
                raise
        # if r_path[-1] == 2021 and q_path[-1] == 5181:
        #     debug=True
        alignment = self.alignPaths2(q_path,r_path)
        # print("comparison",self.compareAlignments(alignment1,alignment))
        q_idx = 0
        r_idx = 0
        continuous_ins = 0
        continuous_del = 0
        diff_flag = False
        q_ops = []
        start_flag = True
        ignore_deletion_flag = False
        for op in alignment:
            if op == 0:
                # if not start_flag:
                #     last_q_match_idx = q_idx - continuous_ins - 1
                #     last_r_match_idx = r_idx - continuous_del - 1
                #     assert q_path[last_q_match_idx] == r_path[last_r_match_idx]
                if debug:
                    print(r_path[r_idx],q_path[q_idx],sep='\t')
                if continuous_ins >= psi:
                    diff_flag = True
                    ignore_deletion_flag = True
                elif continuous_ins > 0:
                    #Correct short insertion\
                    # if debug:
                    #     print(continuous_ins)
                    for i in range(q_idx - continuous_ins,q_idx):
                        q_ops.append(("ins",i,0))
                if continuous_del >= psi:
                    diff_flag = True
                elif continuous_del > 0 and not ignore_deletion_flag:
                    #Correct short deletion
                    # if debug:
                    #     print(continuous_del)
                    for i in range(r_idx - continuous_del,r_idx):
                        q_ops.append(("del",q_idx-continuous_ins,i))
                ignore_deletion_flag = False
                start_flag = False
                q_idx += 1
                r_idx += 1
                continuous_ins = 0
                continuous_del = 0
            elif op == 1:
                # if not start_flag:
                #     continuous_ins += 1
                continuous_ins += 1
                if debug:
                    print("----",q_path[q_idx],sep = '\t')
                q_idx += 1
                
            elif op == 2:
                if not start_flag:
                    continuous_del += 1
                if debug:
                    print(r_path[r_idx],"----",sep = '\t')
                r_idx += 1
                

        if continuous_ins >= psi:
            diff_flag = True
            ignore_deletion_flag = True
        elif continuous_ins > 0:
            #Correct short insertion\
            if debug:
                print("Hitting here for some reason")
            for i in range(q_idx - continuous_ins,q_idx):
                q_ops.append(("ins",i,0))
        assert q_idx == len(q_path)
        assert r_idx == len(r_path)
        # print(list(reversed(q_ops)))
        last_q_idx = np.inf
        last_r_idx = np.inf
        last_op = None
        # if debug:
        #     print("q_ops",q_ops)
        for op,q_idx,r_idx in reversed(sorted(q_ops,key=lambda t: (t[1],t[0],t[2]))):
            # print("here")
            # if last_q_idx <= q_idx:
            #     print("PANIC")
            #     print(last_q_idx,q_idx)
            #     print(last_r_idx,r_idx)
            #     print(last_op,op)
            if op == "ins":
                if q_idx != len(q_path) - 1 and q_idx != 0:
                    # try:
                    assert q_path[q_idx-1] < q_path[q_idx+1]
                    # except:
                    #     print(q_path[q_idx-1])
                    #     print(q_path[q_idx])
                    #     print(q_path[q_idx+1])
                    #     raise
                del q_path[q_idx]
            elif op == "del":
                if q_idx != 0:
                    try:
                        assert q_path[q_idx-1] < r_path[r_idx] < q_path[q_idx]
                    except:
                        # if debug:
                        print(q_path[-1])
                        print(r_path[-1])
                        print("q_ops",q_ops)
                        print(q_path[q_idx-1], r_path[r_idx], q_path[q_idx])
                        raise
                q_path.insert(q_idx,r_path[r_idx])
            last_q_idx = q_idx
            last_r_idx = r_idx
            last_op = op
                
        # alignment = self.alignPaths(q_path,r_path,a,b,o,e,anchor)
        for i in range(len(q_path)-1):
            try:
                assert(q_path[i] < q_path[i+1])
            except:
                print(q_path[i])
                print(q_path[i+1])
                raise
        # if og_q_path == q_path:
        #     print("NOT CORRECTING 2")
        for i in range(len(q_path)-1):
            assert q_path[i] < q_path[i+1]
        return (diff_flag,q_path)
    
    def correctPathsAndEnds(self,q_path,r_path,psi=10):
        # alignment1 = self.alignPaths(q_path,r_path,a,b,o,e,anchor)
        # alignment1 = self.alignPaths(q_path,r_path)
        # og_q_path = copy.deepcopy(q_path)
        alignment = self.alignPaths2(q_path,r_path)
        # print("comparison",self.compareAlignments(alignment1,alignment))
        q_idx = 0
        r_idx = 0
        continuous_ins = 0
        continuous_del = 0
        diff_flag = False
        q_ops = []
        start_flag = True
        ignore_deletion_flag = False
        for op in alignment:
            if op == 0:
                # if not start_flag:
                #     last_q_match_idx = q_idx - continuous_ins - 1
                #     last_r_match_idx = r_idx - continuous_del - 1
                #     assert q_path[last_q_match_idx] == r_path[last_r_match_idx]
                # print(r_path[r_idx],q_path[q_idx],sep='\t')
                if continuous_ins >= psi:
                    diff_flag = True
                    ignore_deletion_flag = True
                elif continuous_ins > 0:
                    #Correct short insertion\
                    for i in range(q_idx - continuous_ins,q_idx):
                        q_ops.append(("ins",i,0))
                if continuous_del >= psi:
                    diff_flag = True
                elif continuous_del > 0 and not ignore_deletion_flag:
                    #Correct short deletion
                    # print(continuous_del)
                    # print(continuous_ins)
                    # print("q_path",q_path[q_idx - 2])
                    # print("q_path",q_path[q_idx - 1])
                    # print("q_path",q_path[q_idx - 0])
                    # print("q_path",q_path[q_idx + 1])
                    # print("q_path",q_path[q_idx + 2])
                    # print("r_path",r_path[r_idx-2])
                    # print("r_path",r_path[r_idx-1])
                    # print("r_path",r_path[r_idx])
                    # print("r_path",r_path[r_idx+1])
                    # print("r_path",r_path[r_idx+2])
                    for i in range(r_idx - continuous_del,r_idx):
                        q_ops.append(("del",q_idx-continuous_ins,i))
                ignore_deletion_flag = False
                start_flag = False
                q_idx += 1
                r_idx += 1
                continuous_ins = 0
                continuous_del = 0
            elif op == 1:
                continuous_ins += 1
                # print("----",q_path[q_idx],sep = '\t')
                q_idx += 1
                
            elif op == 2:
                continuous_del += 1
                # print(r_path[r_idx],"----",sep = '\t')
                r_idx += 1
                
        
        assert q_idx == len(q_path)
        assert r_idx == len(r_path)
        # print("END")
        if continuous_ins >= psi:
            diff_flag = True
            ignore_deletion_flag = True
        elif continuous_ins > 0:
            # print("SMALL_INSERTION")
            #Correct short insertion
            for i in range(q_idx - continuous_ins,q_idx):
                q_ops.append(("ins",i,0))
        if continuous_del >= psi:
            diff_flag = True
        elif continuous_del > 0 and not ignore_deletion_flag:
            # print("SMALL_DELETION")
            #Correct short deletion
            for i in range(r_idx - continuous_del,r_idx ):
                q_ops.append(("del",q_idx-continuous_ins,i))
        # print(list(reversed(q_ops)))
        last_q_idx = np.inf
        last_r_idx = np.inf
        last_op = None
        for op,q_idx,r_idx in reversed(sorted(q_ops,key=lambda t: (t[1],t[0],t[2]))):
            # print("here")
            # if last_q_idx <= q_idx:
            #     print("PANIC")
            #     print(last_q_idx,q_idx)
            #     print(last_r_idx,r_idx)
            #     print(last_op,op)
            # if op == "ins":
            #     print(op,q_idx,r_idx,q_path[q_idx],r_path[r_idx])
            # elif op == "del":
            #     if q_idx == len(q_path):
            #         print(op,q_idx,r_idx,"NA",r_path[r_idx])
            #     else:
            #         print(op,q_idx,r_idx,q_path[q_idx],r_path[r_idx])
            if op == "ins":
                # try:
                if q_idx != 0 and q_idx != len(q_path) - 1:
                    assert q_path[q_idx-1] < q_path[q_idx+1]
                # except:
                # print(q_idx)
                # print(q_path[q_idx])
                # print(q_path[q_idx])
                # raise
                del q_path[q_idx]
            elif op == "del":
                try:
                    if q_idx != 0 and q_idx != len(q_path):
                        assert q_path[q_idx-1] < r_path[r_idx] < q_path[q_idx]
                except:
                    
                    print("panic!")
                    print(q_idx)
                    print(r_idx)
                    print(len(q_path))
                    print(len(r_path))
                    print(q_path[q_idx-5])
                    print(q_path[q_idx-4])
                    print(q_path[q_idx-3])
                    print(q_path[q_idx-2])
                    print(q_path[q_idx-1])
                    print(q_path[q_idx])
                    # print(q_path[q_idx+1])
                    print(r_path[r_idx])
                    raise
                q_path.insert(q_idx,r_path[r_idx])
            last_q_idx = q_idx
            last_r_idx = r_idx
            last_op = op
                
        # alignment = self.alignPaths(q_path,r_path,a,b,o,e,anchor)
        # for i in range(len(q_path)-1):
        #     try:
        #         assert(q_path[i] < q_path[i+1])
        #     except:
        #         print(q_path[i])
        #         print(q_path[i+1])
        # if og_q_path == q_path:
        #     print("NOT CORRECTING 2")
        for i in range(len(q_path)-1):
            try:
                assert q_path[i] < q_path[i+1]
            except:
                print(i,i+1,q_path[i],q_path[i+1])
                raise
        return (diff_flag,q_path)
    
    def trimPaths(self,longest_path):
        longest_path_set = set(longest_path)
        for i,index in enumerate(longest_path):
            if len(self.edges[index]) > 1:
                if i + 1 < len(longest_path):
                    for new_index in set(self.edges[index]) - set([longest_path[i + 1]]):
                        if new_index in longest_path_set:
                            if longest_path.index(new_index) - i < 10:
                                self.removeEdges([index] + [new_index])
                        # new_paths = self.followPath(new_index,longest_path_set)
                        # print(new_index)
                        new_paths = self.dfs_paths_one_return(new_index,longest_path_set)
                        # print(new_path)
                        for new_path in new_paths: 
                            if len(new_path) < 10:
                                # print([index] + new_path)
                                self.removeEdges([index] + new_path)
    
    def removeEdges(self,path):
        for i in range(len(path) - 1):
            if path[i + 1] in self.edges[path[i]]:
                self.edges[path[i]].remove(path[i + 1])
    
    def jsOutput(self,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path=[],collapse=False):
        """returns a list of strings containing a a description of the graph for viz.js, http://visjs.org"""
        
        # get the consensus sequence, which we'll use as the "spine" of the
        # graph
        # path, __, __ = self.consensus()
        # path = self.getLongestPath2(nodes,node_indexes,source_nodes,edges,weights)
        # path = self.getMostLikelyPath(nodes,node_indexes,source_nodes,end_nodes,edges,log_likelihoods)
        path = sorted(self.greedyPathWalk(nodes,node_indexes,source_nodes,edges,weights),key = lambda t: len(t),reverse=True)[0]
        # print(path)
        highlighted_nodes1 = set()
        highlighted_edges1 = set()
        highlighted_nodes2 = set()
        highlighted_edges2 = set()
        pathdict = {}
        for i, nodeID in enumerate(path):
            if collapse:
                pathdict[nodeID] = i*500
            else:
                pathdict[nodeID] = i*150
            highlighted_nodes1.add(nodeID)
            if i != 0:
                highlighted_edges1.add((path[i-1],nodeID))
        
        if len(highlighted_path) != 0:
            for i,nodeID in enumerate(highlighted_path):
                highlighted_nodes2.add(nodeID)
                if i != 0:
                    highlighted_edges2.add((highlighted_path[i-1],nodeID))
        # print(highlighted_nodes)
        lines = ['var nodes = [']

        # ni = self.nodeiterator()
        count = 0
        last_j_in_path_dict = sorted(pathdict.keys())[0]
        for i,node in enumerate(nodes):
            j = node_indexes[('f',i)]
            # line = '    {id:'+str(j)+', label: "'+node.char+'"'
            # line = '    {id:'+str(j)+', label: "'+"%0.1f" %(node.log_likelihood_scaled)+'"'
            line = '    {id:'+str(j)+', label: "'+"%d,%0.1f,%s" %(j,node.log_likelihood_scaled,node.char)+'"'
            # line = '    {id:'+str(j)+', label: "'+"%0.1f" %(node.log_likelihood)+'"'
            # line = '    {id:'+str(j)+', label: "'+str(j)+'"'
            # if j == 1618:
            #     print(j in highlighted_nodes)
            if j in highlighted_nodes1:
                # print(j)
                if j in highlighted_nodes2:
                    highlight = ', color: \'purple\''
                else:
                    highlight = ', color: \'red\''
            elif j in highlighted_nodes2:
                highlight = ', color: \'blue\''
            else:
                highlight = ', color: \'#D3D3D3\''
            if j in pathdict and count % 5 == 0:
                if collapse:
                    line += ', allowedToMoveX: true, x: ' + str(pathdict[j]) + ', y: 0 ' + highlight + ', allowedToMoveY: true},'
                else:
                    line += ', allowedToMoveX: false, x: ' + str(pathdict[j]) + ', y: 0 ' + highlight + ', allowedToMoveY: true},'
                last_j_in_path_dict = j
            else:
                # print(last_j_in_path_dict)
                line += ', allowedToMoveX: true, x: ' + str(pathdict[last_j_in_path_dict]) + highlight + ', y: 0 , allowedToMoveY: true},'
            lines.append(line)

        lines[-1] = lines[-1][:-1]
        lines.append('];')

        lines.append(' ')

        lines.append('var edges = [')
        # ni = self.nodeiterator()
        # for node in ni():
        written_set = set()
        for i,node in enumerate(nodes):
            k = node_indexes[('f',i)]
            # nodeID = str(node.ID)
            nodeID = str(k)
            for edge in edges[k]:
            # for edge in node.outEdges:
                target = str(edge)
                # weight = str(len(node.outEdges[edge].labels)+1)
                weight = str(weights[(k,edge)])
                if (k,edge) not in written_set:
                    # if (k,edge) == (1618,1623):
                    #     print((k,edge) in highlighted_edges)
                    if (k,edge) in highlighted_edges1:
                        if (k,edge) in highlighted_edges2:
                            highlight = ', color: \'purple\''
                        else:
                            highlight = ', color: \'red\''
                        # print('highlighted',k,edge)
                    elif (k,edge) in highlighted_edges2:
                        highlight = ', color: \'blue\''
                    else:
                        highlight = ', color: \'black\''
                        # print('not highlighted',k,edge)
                    # if (k,edge) == (1618,1623):
                    #     print(highlight)
                    lines.append('    {from: '+nodeID+', to: '+target + highlight + ', value: '+weight+', arrows:\'to, from, middle\' },')
                    written_set.add((k,edge))
            # for alignededge in node.alignedTo:
            #     # These edges indicate alignment to different bases, and are
            #     # undirected; thus make sure we only plot them once:
            #     if node.ID > alignededge:
            #         continue
            #     target = str(alignededge)
            #     lines.append('    {from: '+nodeID+', to: '+target+', value: 1, style: "dash-line"},')
        lines[-1] = lines[-1][:-1]
        lines.append('];')
        return lines
    
    def htmlOutput(self, outfile,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path=[],collapse=False):
        header = """
                  <!doctype html>
                  <html>
                  <head>
                    <title>POA Graph Alignment</title>
                    <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/vis/3.11.0/vis.min.js"></script>
                  </head>
                  <body>
                  <div id="mynetwork"></div>
                  <script type="text/javascript">
                    // create a network
                  """
        outfile.write(textwrap.dedent(header[1:]))
        if collapse:
            nodes,node_indexes,edges,weights,log_likelihoods,end_nodes = self.collapseLinearStretches(nodes,node_indexes,edges,weights,log_likelihoods)
        lines = self.jsOutput(nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods,highlighted_path=highlighted_path,collapse=collapse)
        for line in lines:
            outfile.write(line+'\n')
        footer = """
                  var container = document.getElementById('mynetwork');
                  var data= {
                    nodes: nodes,
                    edges: edges,
                  };
                  var options = {
                    width: '100%',
                    height: '800px'
                  };
                  var network = new vis.Network(container, data, options);
                </script>
                </body>
                </html>
                """
        outfile.write(textwrap.dedent(footer))
        outfile.close()
    
    def assignReadsToPaths(self,reads,paths):
        closestPaths = []
        for read in reads:
            path_scores = []
            for path in paths:
                path_scores.append(self.scorePaths(read.path,path))
            closestPaths.append(max(range(len(path_scores)),key=path_scores.__getitem__))
        return closestPaths
    
    def getRepresentativeGraph(self,outfile1=None,outfile2=None,outfile3=None,output_every_iteration=True,outfile_prefix = "testing",psi=10):
        # paths = self.getRepresentativePaths3(outfile1,outfile2,output_every_iteration=output_every_iteration,outfile_prefix = outfile_prefix,psi=psi)
        paths = self.getRepresentativePaths4(outfile1,outfile2,psi=psi)
        closest = self.assignReadsToPaths(self.reads,paths)
        closest_set = set()
        for c in closest:
            closest_set.add(c)
        print(closest)
        selected_paths = []
        for c in closest_set:
            selected_paths.append(paths[c])
        nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructGraphFromPaths(selected_paths)
        self.htmlOutput(outfile3,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods)
    
    def getCorrectedReads(self,outfile1=None,outfile2=None,outfile3=None,outfile4=None,output_every_iteration=True,outfile_prefix = "testing",psi=10):
        paths = self.getRepresentativePaths3(outfile1,outfile2,output_every_iteration=output_every_iteration,outfile_prefix = outfile_prefix,psi=psi)
        closest = self.assignReadsToPaths(self.reads,paths)
        # closest_set = set()
        # for c in closest:
        #     closest_set.add(c)
        # print(closest)
        selected_paths = []
        for c in closest:
            selected_paths.append(paths[c])
        nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructGraphFromPaths(selected_paths)
        self.htmlOutput(outfile3,nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods)
        self.alignReadsToPaths(self.reads,selected_paths,outfile=outfile4)
    
    def alignReadsToPaths(self,reads,paths,outfile = None):
        assert len(reads) == len(paths)
        read_seqs = []
        for read in reads:
            read_seq = []
            for node_idx in read.path:
                read_seq.append(self.nodes[node_idx].char)
            read_seqs.append(read_seq)
        path_seqs = []
        for path in paths:
            path_seq = []
            for node_idx in path:
                path_seq.append(self.nodes[node_idx].char)
            path_seqs.append(path_seq)
        for i in range(len(path_seqs)):
            corrected_read = self.correctSeq(read_seqs[i],path_seqs[i])
            if outfile is not None:
                outfile.write(">" + reads[i].name + "\n")
                outfile.write(''.join(corrected_read) + "\n")
        outfile.close()
    
    def correctSeq(self,read_seq,ref_seq,psi=10):
        alignment = self.alignSeqs(read_seq,ref_seq)
        
        q_idx = 0
        r_idx = 0
        continuous_ins = 0
        continuous_del = 0
        diff_flag = False
        q_ops = []
        start_flag = True
        ignore_deletion_flag = False
        for op in alignment:
            if op == 0:
                # if not start_flag:
                #     last_q_match_idx = q_idx - continuous_ins - 1
                #     last_r_match_idx = r_idx - continuous_del - 1
                #     assert read_seq[last_q_match_idx] == ref_seq[last_r_match_idx]
                # print(ref_seq[r_idx],read_seq[q_idx],sep='\t')
                read_seq[q_idx] = ref_seq[r_idx]
                if continuous_ins >= psi:
                    diff_flag = True
                    ignore_deletion_flag = True
                elif continuous_ins > 0:
                    #Correct short insertion\
                    for i in range(q_idx - continuous_ins,q_idx):
                        q_ops.append(("ins",i,0))
                if continuous_del >= psi:
                    diff_flag = True
                elif continuous_del > 0 and not ignore_deletion_flag:
                    #Correct short deletion
                    for i in range(r_idx - continuous_del,r_idx):
                        q_ops.append(("del",q_idx-continuous_ins,i))
                ignore_deletion_flag = False
                start_flag = False
                q_idx += 1
                r_idx += 1
                continuous_ins = 0
                continuous_del = 0
            elif op == 1:
                if not start_flag:
                    continuous_ins += 1
                    # print("----",read_seq[q_idx],sep = '\t')
                q_idx += 1
                
            elif op == 2:
                if not start_flag:
                    continuous_del += 1
                    # print(ref_seq[r_idx],"----",sep = '\t')
                r_idx += 1
                
        
        assert q_idx == len(read_seq)
        assert r_idx == len(ref_seq)
        # print(list(reversed(q_ops)))
        last_q_idx = np.inf
        last_r_idx = np.inf
        last_op = None
        for op,q_idx,r_idx in reversed(sorted(q_ops,key=lambda t: (t[1],t[0],t[2]))):
            # print("here")
            # if last_q_idx <= q_idx:
            #     print("PANIC")
            #     print(last_q_idx,q_idx)
            #     print(last_r_idx,r_idx)
            #     print(last_op,op)
            if op == "ins":
                # assert read_seq[q_idx-1] < read_seq[q_idx+1]
                del read_seq[q_idx]
            elif op == "del":
                # assert read_seq[q_idx-1] < ref_seq[r_idx] < read_seq[q_idx]
                read_seq.insert(q_idx,ref_seq[r_idx])
            last_q_idx = q_idx
            last_r_idx = r_idx
            last_op = op
                
        # alignment = self.alignPaths(read_seq,ref_seq,a,b,o,e,anchor)
        # for i in range(len(read_seq)-1):
        #     try:
        #         assert(read_seq[i] < read_seq[i+1])
        #     except:
        #         print(read_seq[i])
        #         print(read_seq[i+1])
        # if og_read_seq == read_seq:
        #     print("NOT CORRECTING 2")
        return read_seq
    
    def getpHMM(self,outfile1=None,outfile2=None,outfile3=None,outfile4=None,output_every_iteration=True,outfile_prefix = "testing",psi=10):
        # for i in range(len(self.reads)):
        #     self.reads[i].corrected_path = self.reads[i].path.copy()
        # nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructNewGraph(self.reads)
        # nodes,node_indexes,edges,weights,log_likelihoods,end_nodes = self.collapseLinearStretches(nodes,node_indexes,edges,weights,log_likelihoods)
        
        # paths = self.getRepresentativePaths3(outfile1,outfile2,output_every_iteration=output_every_iteration,outfile_prefix = outfile_prefix,psi=psi)
        paths = self.getRepresentativePaths4(outfile1,outfile2,psi=psi)
        closest = self.assignReadsToPaths(self.reads,paths)
        # closest_set = set()
        # for c in closest:
        #     closest_set.add(c)
        # print(closest)
        selected_paths = []
        for c in closest:
            selected_paths.append(paths[c])
            # for idx in paths[c]:
            #     print(idx)
        # print(len(selected_paths)
        nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods = self.constructGraphFromPaths(selected_paths)
        for node in nodes:
            node.start_node_flag = False
            node.end_node_flag = False
        for read in self.reads:
            nodes[node_indexes['b',read.corrected_path[0]]].start_node_flag = True
            nodes[node_indexes['b',read.corrected_path[-1]]].end_node_flag = True
        # print(end_nodes)
        # for read in self.reads:
        #     print(read.path[-1])
        #     assert read.path[-1] in end_nodes
        # for path in selected_paths:
        #     print(path[-1])
        #     assert path[-1] in end_nodes
        # for idx in sorted(node_indexes.keys()):
        #     print("idx,", idx, node_indexes[idx])
        nodes,node_indexes,edges,weights,log_likelihoods,end_nodes = self.collapseLinearStretches(nodes,node_indexes,edges,weights,log_likelihoods)
        # print(edges)
        # print(node_indexes)
        # print()
        # self.htmlOutput(open("test_collapse.html",'w'),nodes,node_indexes,source_nodes,end_nodes,edges,weights,log_likelihoods)
        print("FUCK1")
        # for idx in sorted(node_indexes.keys()):
        #     print("idx,", idx, node_indexes[idx])
        read_seqs = []
        read_ids = []
        for read in self.reads:
            read_seqs.append("".join(read.seq))
            read_ids.append(read.name)
        # print(end_nodes)
        # for node in end_nodes:
        #     print(node_indexes[('b',node)])
        return splitProfileHMM(nodes,node_indexes,edges,selected_paths,read_seqs,read_ids,source_nodes,end_nodes,outfile4)
    
test = poGraph()
# print(test)
## PO output is already topologically sorted
# test.makeGraph(open("/Users/nproach/Documents/src/poaV2/test.po"))
# test.makeGraph(open("/Users/nproach/Documents/src/poaV2/mvk-1_reads.3.po"))
time1 = time.time()
test.makeGraph(open("/Users/nproach/Documents/src/poaV2/mvk-1_reads.4.po"))
time2 = time.time()
print(time2 - time1)
# test.getRepresentativePaths4(outfile1 = open("testing_p_queue_after.html",'w'),outfile2 = open("testing_p_queue_before.html",'w'),psi=10)
# test.makeGraph(open("/Users/nproach/Documents/src/poaV2/test_phmm.po"))
start_time = time.time()
phmm = test.getpHMM(open("before.html",'w'),open("after.html",'w'),open("representative_graph.html",'w'),open("testing_phmm.fa",'w'),output_every_iteration=True,outfile_prefix = "testing",psi=10)
print("---%s--- seconds" %(time.time() - start_time))
start_time = time.time()
phmm.viterbi_score(test.reads[0].seq,max_threshold = None,min_threshold= None)
print("---%s--- seconds" %(time.time() - start_time))

start_time = time.time()
phmm.viterbi_score2(test.reads[0].seq,max_threshold = None,min_threshold= None,window=20)
print("---%s--- seconds" %(time.time() - start_time))

# print(test.edges)
# print(test.source_node_indices)
# longest_path = test.getLongestPath()
# print(longest_path)

# paths = test.getRepresentativePaths3(open("before.html",'w'),open("after.html",'w'),output_every_iteration=True,outfile_prefix = "testing",psi=10)
# test.getRepresentativeGraph(open("before.html",'w'),open("after.html",'w'),open("representative_graph.html",'w'),output_every_iteration=True,outfile_prefix = "testing",psi=15)
# test.getCorrectedReads(open("before.html",'w'),open("after.html",'w'),open("representative_graph.html",'w'),open("testing_consensus.3.fa",'w'),output_every_iteration=True,outfile_prefix = "testing",psi=10)

# test_seq = []

# for _ in range(2555):
#     rn = np.random.random()
#     if rn < 0.25:
#         test_seq.append("A")
#     elif rn < 0.5:
#         test_seq.append("T")
#     elif rn < 0.75:
#         test_seq.append("C")
#     else:
#         test_seq.append("G")

# start_time = time.time()
# phmm.viterbi_score(''.join(test_seq))
# print("---%s--- seconds" %(time.time() - start_time))


#

# splitProfileHMM.find_windows = find_windows
# phmm.windows = {}
# for end_node in phmm.end_nodes:
#     if end_node < phmm.num_match_states:
#         print(end_node)
#         phmm.find_windows(end_node,len(test.reads[0].seq))
# closest = test.assignReadsToPaths(test.reads,paths)
# closest_set = set()
# for c in closest:
#     closest_set.add(c)
# selected_paths = []
# for c in closest_set:
#     selected_paths.append(paths[c])



# paths,read1,read2 = test.getRepresentativePaths2(open("before.html",'w'),open("after.html",'w'),debug=True)#,output_every_iteration=True,outfile_prefix = "testing")

# test.trimPaths(longest_path)
# print(test.getLongestPath())
# print(sorted(test.getPossiblePaths(),key=len,reverse=True))

# paths = sorted(test.getPossiblePaths(),key=len,reverse=True)
# print(paths)
# paths = test.trimPaths(paths[0])
# print(sorted(test.getPossiblePaths(),key=len,reverse=True))

# print(test.findLongestPath())