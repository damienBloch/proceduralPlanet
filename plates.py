import numpy as np
import networkx as nx
import sys

sys.setrecursionlimit(10**6) 

class Plates:
    def __init__(self,centers,seed=0):
        self.centers=centers.copy()
        self.rng=np.random.RandomState(seed)
    def _seed(self,numberPlates):
        seeds=self.rng.choice(np.arange(self.centers.number_of_nodes()),size=numberPlates,replace=False)
        for i in range(self.centers.number_of_nodes()):
            self.centers.nodes[i]["plate"]=-1
        for i in range(numberPlates):
            self.centers.nodes[seeds[i]]["plate"]=i
        return seeds
    def _colorPlates(self,queue):
        if(len(queue)==0):
            return
        queue=list(queue)
        if(len(queue)>0):
            index=self.rng.choice(queue)
            neighbors=list(self.centers.adj[index])
            uncoloredNeighbors=[]
            for n in neighbors:
                if(self.centers.nodes[n]["plate"]<0):
                    uncoloredNeighbors.append(n)
            if(len(uncoloredNeighbors)>0):
                newIndex=self.rng.choice(uncoloredNeighbors)
                self.centers.nodes[newIndex]["plate"]=self.centers.nodes[index]["plate"]
                queue.append(newIndex)
            else:
                queue.remove(index)
        self._colorPlates(queue)
    def generatePlates(self,numberPlates):
        seeds=self._seed(numberPlates)
        self._colorPlates(seeds)
        
        self.plates=[]
        platesNumber=np.array([p for _,p in self.centers.nodes.data("plate")],dtype=int)
        for i in range(numberPlates):
            sub=[i for i,p in enumerate(platesNumber) if p==i]
            self.plates.append(nx.subgraph(self.centers,list(sub)))        
        return self.centers,self.plates
        
        