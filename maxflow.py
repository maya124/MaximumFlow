"""
File: maxflow_repetitive.py
Author: Maya Varma
Usage: python maxflow_repetitive.py $set_order
Inputs: set_order is the order of variant sets expressed as 1_2_3_4_5
Description: Runs the maxflow algorithm.
"""
import itertools
import pickle
import numpy as np
from collections import defaultdict 
from operator import xor
import sys
   
#This class represents a directed graph using adjacency matrix representation 
class Graph: 
   
    def __init__(self,graph): 
        self.graph = graph # residual graph 
        self. ROW = len(graph) 
           
   
    '''Returns true if there is a path from source 's' to sink 't' in 
    residual graph. Also fills parent[] to store the path '''
    def BFS(self,s, t, parent): 
  
        # Mark all the vertices as not visited 
        visited =[False]*(self.ROW) 
          
        # Create a queue for BFS 
        queue=[] 
          
        # Mark the source node as visited and enqueue it 
        queue.append(s) 
        visited[s] = True
           
        # Standard BFS Loop 
        while queue: 
  
            #Dequeue a vertex from queue and print it 
            u = queue.pop(0) 
          
            # Get all adjacent vertices of the dequeued vertex u 
            # If a adjacent has not been visited, then mark it 
            # visited and enqueue it 
            for ind, val in enumerate(self.graph[u]):
                if visited[ind] == False and val > 0 : 
                    queue.append(ind) 
                    visited[ind] = True
                    parent[ind] = u 
  
        # If we reached sink in BFS starting from source, then return 
        # true, else false 
        return True if visited[t] else False
              
      
    # Returns the maximum flow from s to t in the given graph 
    def FordFulkerson(self, source, sink, labels, order): 
  
        # This array is filled by BFS and to store path 
        parent = [-1]*(self.ROW) 
  
        max_flow = 0 # There is no flow initially 
    
        # Augment the flow while there is path from source to sink 
        while self.BFS(source, sink, parent) : 
  
            # Find minimum residual capacity of the edges along the 
            # path filled by BFS. Or we can say find the maximum flow 
            # through the path found. 
            path_flow = float("Inf") 
            s = sink 
            pathString = ''
            pathString += 'node %s, ' % labels[s]
            while(s !=  source): 
                path_flow = min (path_flow, self.graph[parent[s]][s]) 
                s = parent[s] 
                pathString += 'node %s, ' % labels[s]
            pathString += '\n'
            with open("path.txt", "a+") as f:
                f.write(pathString)

            # Add path flow to overall flow 
            max_flow +=  path_flow 
  
            # update residual capacities of the edges and reverse edges 
            # along the path 
            v = sink 
            while(v !=  source): 
                u = parent[v] 
                self.graph[u][v] -= path_flow 
                self.graph[v][u] += path_flow 
                v = parent[v] 

        return max_flow 
    
    
def loadPairToLD():
    '''
    Load in all linkage disequilibrium results. Each pairToLD variable consists of a map in this
    format: pairToLD[(chr, var1, var2)] = [ldVal1, ldVal2,...]
    '''
    pairToLD = {}
    def compilePairToLD(set1, set2, ldList):
        for elem in ldList:
            pairToLD[(str(set1)+":"+str(set2), elem[0], elem[1], elem[2])] = ldList[elem]
            pairToLD[(str(set2)+":"+str(set1), elem[0], elem[2], elem[1])] = ldList[elem]

    file = open("pairToLD_R2_1_2",'rb')
    pairToLD12 = pickle.load(file)
    compilePairToLD(1, 2, pairToLD12)
    file.close()
    file = open("pairToLD_R2_2_3",'rb')
    pairToLD23 = pickle.load(file)
    compilePairToLD(2, 3, pairToLD23)
    file.close()
    file = open("pairToLD_R2_3_4",'rb')
    pairToLD34 = pickle.load(file)
    compilePairToLD(3, 4, pairToLD34)
    file.close()
    file = open("pairToLD_R2_4_5",'rb')
    pairToLD45 = pickle.load(file)
    compilePairToLD(4, 5, pairToLD45)
    file.close()
    file = open("pairToLD_R2_1_3",'rb')
    pairToLD13 = pickle.load(file)
    compilePairToLD(1, 3, pairToLD13)
    file.close()
    file = open("pairToLD_R2_1_4",'rb')
    pairToLD14 = pickle.load(file)
    compilePairToLD(1, 4, pairToLD14)
    file.close()
    file = open("pairToLD_R2_1_5",'rb')
    pairToLD15 = pickle.load(file)
    compilePairToLD(1, 5, pairToLD15)
    file.close()
    file = open("pairToLD_R2_2_4",'rb')
    pairToLD24 = pickle.load(file)
    compilePairToLD(2, 4, pairToLD24)
    file.close()
    file = open("pairToLD_R2_2_5",'rb')
    pairToLD25 = pickle.load(file)
    compilePairToLD(2, 5, pairToLD25)
    file.close()
    file = open("pairToLD_R2_3_5",'rb')
    pairToLD35 = pickle.load(file)
    compilePairToLD(3, 5, pairToLD35)
    file.close()

    return pairToLD

def constructAdjacency(pairToLD, i1, i2, i3, i4, i5):
    '''
    Construct the adjacency matrix for the max flow network. Returns matrix.
    '''
    set1 = set()
    set2 = set()
    set3 = set()
    set4 = set()
    set5 = set()
    for elem in pairToLD:
        if(elem[0]==str(i1)+":"+str(i2)):
            set1.add(elem[1]+"_"+elem[2])
            set2.add(elem[1]+"_"+elem[3])
        if(elem[0]==str(i2)+":"+str(i3)):
            set2.add(elem[1]+"_"+elem[2])
            set3.add(elem[1]+"_"+elem[3])
        if(elem[0]==str(i3)+":"+str(i4)):
            set3.add(elem[1]+"_"+elem[2])
            set4.add(elem[1]+"_"+elem[3])
        if(elem[0]==str(i4)+":"+str(i5)):
            set4.add(elem[1]+"_"+elem[2])
            set5.add(elem[1]+"_"+elem[3])

    set1_in = list(set1)
    set1_out = list(set1)
    set2_in = list(set2)
    set2_out = list(set2)
    set3_in = list(set3)
    set3_out = list(set3)
    set4_in = list(set4)
    set4_out = list(set4)
    set5_in = list(set5)
    set5_out = list(set5)

    #Create adjacency matrix from two sets
    adj = np.zeros((len(set1_in)+len(set1_out)+len(set2_in)+len(set2_out)+len(set3_in)+len(set3_out)+ \
                len(set4_in)+len(set4_out)+len(set5_in)+len(set5_out)+2, \
                len(set1_in)+len(set1_out)+len(set2_in)+len(set2_out)+len(set3_in)+len(set3_out)+ \
                len(set4_in)+len(set4_out)+len(set5_in)+len(set5_out)+2))
    labels = ['1_%s_in'% x for x in set1_in]+['1_%s_out'% x for x in set1_out]+\
    ['2_%s_in'%x for x in set2_in]+['2_%s_out'%x for x in set2_out] + \
    ['3_%s_in'%x for x in set3_in]+['3_%s_out'%x for x in set3_out]+ \
    ['4_%s_in'%x for x in set4_in]+['4_%s_out'%x for x in set4_out] + \
    ['5_%s_in'%x for x in set5_in]+['5_%s_out'%x for x in set5_out]

    labels.insert(0,'source')
    labels.append('sink')

    #Construct an edge between two nodes if one of the following conditions is met:
    #     - LD score is greater than 0.8, suggesting that the two variants are in LD
    #     - the same variant exists in both sets

    def addEdgesAcrossSets(set1_out, set2_in, pairToLD, setNum1, setNum2, id1, id2):
        '''
        Input: set1_out, set2_in - lists of variants (we are connecting out_nodes to in_nodes)
               pairToLD - mapping from variant pairs to LD scores
               setNum1, setNum2 - strings indicating set numbers
               id1, id2 - strings representing original cross val set orders for pairToLD lookup
        '''
        for pair in list(itertools.product(set1_out, set2_in)):
            chrom1 = pair[0].split("_")[0]
            chrom2 = pair[1].split("_")[0]
            if(chrom1 != chrom2): continue
            index1 = labels.index(setNum1+'_'+pair[0]+'_out')
            index2 = labels.index(setNum2+'_'+pair[1]+'_in')
            var1 = pair[0].split("_")[1]
            var2 = pair[1].split("_")[1]
            if((id1+':'+id2, chrom1, var1, var2) in pairToLD and min(pairToLD[(id1+':'+id2,chrom1, var1, var2)])>0.8):
                adj[index1, index2] = 1
            if(var1==var2):
                adj[index1, index2] = 1
            
    def addEdgesWithinSets(set_in, set_out, id1):
        '''
        Input: set_in, set_out - lists of variants (we are connecting in_nodes to out_nodes in the same set)
        id1 - string indicating set number
        '''
        for elem in set_in:
            index1 = labels.index(id1+'_'+elem+'_in')
            index2 = labels.index(id1+'_'+elem+'_out')
            adj[index1, index2] = 1
    
    addEdgesWithinSets(set1_in, set1_out, '1')
    addEdgesAcrossSets(set1_out, set2_in, pairToLD, '1', '2', str(i1), str(i2))
    addEdgesWithinSets(set2_in, set2_out, '2')
    addEdgesAcrossSets(set2_out, set3_in, pairToLD, '2', '3', str(i2), str(i3))
    addEdgesWithinSets(set3_in, set3_out, '3')
    addEdgesAcrossSets(set3_out, set4_in, pairToLD, '3', '4', str(i3), str(i4))
    addEdgesWithinSets(set4_in, set4_out, '4')
    addEdgesAcrossSets(set4_out, set5_in, pairToLD, '4', '5', str(i4), str(i5))
    addEdgesWithinSets(set5_in, set5_out, '5')

    #Add edges from the source to all variants in set1_in
    for elem in set1_in:
        adj[labels.index('source'), labels.index('1_'+elem+'_in')] = 1
    
    #Add edges from the variants in set5_out to the sink
    for elem in set5:
        adj[labels.index('5_'+elem+'_out'), labels.index('sink')] = 1

    print(adj.shape)
    print(np.count_nonzero(adj))
    return adj, labels


def main():
    pairToLD = loadPairToLD()
    order = sys.argv[1]
    order_list = order.split("_")
    adj, labels = constructAdjacency(pairToLD, int(order_list[0]), int(order_list[1]), int(order_list[2]), int(order_list[3]), int(order_list[4]))
    g = Graph(adj)

    source = 0; sink = labels.index('sink')
    result = "The maximum possible flow is %d " % g.FordFulkerson(source, sink, labels, order) 
    with open("results.txt", "w") as f:
        f.write(result)
    
if __name__=="__main__":
    main()
