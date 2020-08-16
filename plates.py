import numpy as np
import networkx as nx
import sys
from tqdm import tqdm

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
            sub=[k for k,p in enumerate(platesNumber) if p==i]
            self.plates.append(nx.subgraph(self.centers,list(sub)))        
        return self.centers,self.plates
    def randomSpeed(self,plates,centers):
        numberPlates=len(plates)
        platesRotationAxes=self.rng.normal(size=(numberPlates,3))
        norm=np.linalg.norm(platesRotationAxes,axis=-1)
        for i in range(3):
            platesRotationAxes[:,i]/=norm
        for i,a in enumerate(platesRotationAxes):
            plates[i].rotationAxis=platesRotationAxes[i]
            plates[i].rotationSpeed=self.rng.rand()*100
        for i in centers.nodes:
            n=centers.nodes[i]
            p=plates[n["plate"]]
            n["speed"]=np.cross(n["center"],p.rotationAxis)*p.rotationSpeed
        return centers,plates
    def platesElevation(self,plates,centers,flood=0.7):
        for p in plates:
            if(self.rng.rand()<flood):
                p.elevation=self.rng.normal(loc=-7.5,scale=0.5)
            else:
                p.elevation=self.rng.normal(loc=1.5,scale=0.2)
        for i in centers.nodes:
            n=centers.nodes[i]
            n["elevation"]=plates[n["plate"]].elevation
        return centers,plates
    def computePressure(self,plates,centers,corners,parameters):
        border=[(i,j,k,l) for (i,j,(k,l)) in corners.edges.data("separates") if centers.nodes[k]["plate"]!=centers.nodes[l]["plate"]]
        for (i,j,k,l) in tqdm(border,desc="Computing plates collision"):
            normal=centers.nodes[k]["center"]-centers.nodes[l]["center"]
            normal/=np.linalg.norm(normal)
            vrel=centers.nodes[k]["speed"]-centers.nodes[l]["speed"]
            vrelNorm=np.dot(normal,vrel)
            if("pressure" in centers.nodes[k]):
                centers.nodes[k]["pressure"]-=vrelNorm/200*corners.edges[i,j]["length"]*parameters.radius
            else:
                centers.nodes[k]["pressure"]=-vrelNorm/200*corners.edges[i,j]["length"]*parameters.radius

            if("pressure" in centers.nodes[l]):
                centers.nodes[l]["pressure"]-=vrelNorm/200*corners.edges[i,j]["length"]*parameters.radius
            else:
                centers.nodes[l]["pressure"]=-vrelNorm/200*corners.edges[i,j]["length"]*parameters.radius
        return centers,plates
    def computeElevation(self,plates,centers,corners,parameters):
        border=[(i,j,k,l) for (i,j,(k,l)) in corners.edges.data("separates") if centers.nodes[k]["plate"]!=centers.nodes[l]["plate"]]
        for (i,j,k,l) in tqdm(border,desc="Compute elevation"):
            node_k=centers.nodes[k]
            node_l=centers.nodes[l]
            plate_k=plates[node_k["plate"]]
            plate_l=plates[node_l["plate"]]
            #continental collision
            if(plate_k.elevation>0 and plate_l.elevation>0):
                addedElevation=0
                if node_k["pressure"]>0:
                    addedElevation+=node_k["pressure"]/(node_k["area"]*parameters.radius**2)**.5*3/2            
                if node_l["pressure"]>0:
                    addedElevation+=node_l["pressure"]/(node_l["area"]*parameters.radius**2)**.5*3/2
                for onSamePlate in node_k["distances"]:
                        centers.nodes[onSamePlate]["elevation"]+=np.exp(-(node_k["distances"][onSamePlate]*parameters.radius)**2/2/400**2)*addedElevation
                for onSamePlate in node_l["distances"]:
                    centers.nodes[onSamePlate]["elevation"]+=np.exp(-(node_l["distances"][onSamePlate]*parameters.radius)**2/2/400**2)*addedElevation
            #subduction
            def subduction(node_k,node_l):
                if node_k["pressure"]>0:
                    addedElevation=node_k["pressure"]/(node_k["area"]*parameters.radius**2)**.5*3
                    for onSamePlate in node_k["distances"]:
                        centers.nodes[onSamePlate]["elevation"]+=np.exp(-(node_k["distances"][onSamePlate]*parameters.radius)**2/2/300**2)*addedElevation
                if node_l["pressure"]>0:
                    addedElevation=-node_l["pressure"]/(node_l["area"]*parameters.radius**2)**.5
                    for onSamePlate in node_l["distances"]:
                        centers.nodes[onSamePlate]["elevation"]+=1/(1+node_l["distances"][onSamePlate]*parameters.radius/600)*addedElevation
                    addedElevation=node_k["elevation"]-node_l["elevation"]   
                    for onSamePlate in node_l["distances"]:
                        centers.nodes[onSamePlate]["elevation"]+=1/(1+node_l["distances"][onSamePlate]*parameters.radius/200)*addedElevation
            if(plate_k.elevation>0 and plate_l.elevation<=0):
                subduction(node_k,node_l)
            if(plate_l.elevation>0 and plate_k.elevation<=0):
                subduction(node_l,node_k)
            #dorsale
            addedElevation=1
            if(plate_k.elevation<0 and plate_l.elevation<0):
                if(node_k["pressure"]<0 and node_l["pressure"]<0):
                    elevation_difference=(node_k["elevation"]-node_l["elevation"])*.5
                    addedElevation-=node_k["pressure"]/(node_k["area"]*parameters.radius**2)**.5/2*.3
                    addedElevation-=node_l["pressure"]/(node_l["area"]*parameters.radius**2)**.5/2*.3
                    for onSamePlate in node_l["distances"]:
                        centers.nodes[onSamePlate]["elevation"]+=1/(1+node_l["distances"][onSamePlate]*parameters.radius/200)**1.5*(elevation_difference+addedElevation)
                    for onSamePlate in node_k["distances"]:
                        centers.nodes[onSamePlate]["elevation"]+=1/(1+node_k["distances"][onSamePlate]*parameters.radius/200)**1.5*(-elevation_difference+addedElevation)
            
        return centers,plates

        