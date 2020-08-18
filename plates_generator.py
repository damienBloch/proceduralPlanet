import numpy as np
import networkx as nx
import sys
import types
from utils import Projection
from matplotlib import collections as mc
import cartopy.crs as ccrs
from ortools.algorithms import pywrapknapsack_solver

sys.setrecursionlimit(10**6) 

def seedsPlates(planet,planetGenerator):
    planet.plates=nx.Graph()
    planet.plates.add_nodes_from(planet.meshCenters.nodes)
    seeds=planetGenerator.random.choice(np.arange(planet.plates.number_of_nodes()),size=planetGenerator.numberPlates,replace=False)  
    for i in planet.plates.nodes:
        planet.plates.nodes[i]["plate"]=-1
    for i in range(planetGenerator.numberPlates):
        planet.plates.nodes[seeds[i]]["plate"]=i
        
def _colorPlates(planet,random,queue):
    if(len(queue)==0):
        return
    else:
        index=random.choice(queue)
        neighbors=list(planet.meshCenters.adj[index])
        uncoloredNeighbors=[]
        for n in neighbors:
            if(planet.plates.nodes[n]["plate"]<0):
                uncoloredNeighbors.append(n)
        if(len(uncoloredNeighbors)>0):
            newIndex=random.choice(uncoloredNeighbors)
            planet.plates.nodes[newIndex]["plate"]=planet.plates.nodes[index]["plate"]
            queue.append(newIndex)
        else:
            queue.remove(index)
    _colorPlates(planet,random,queue)
def colorPlates(planet,planetGenerator):
    _colorPlates(planet,planetGenerator.random,[i for i,p in planet.plates.nodes.data("plate") if p>=0])
    
    planet.platesSubgraphs=[]
    platesNumber=np.array([p for _,p in planet.plates.nodes.data("plate")],dtype=int)
    for i in range(planetGenerator.numberPlates):
        sub=[k for k,p in enumerate(platesNumber) if p==i]
        planet.platesSubgraphs.append(nx.subgraph(planet.meshCenters,list(sub)))   
    
    def plotPlatesBorders(self,ax,*args,**kwargs):
        segments=np.array([[self.meshCorners.nodes[i]["position"],self.meshCorners.nodes[j]["position"]] for i,j,(pa,pb) in self.meshCorners.edges.data("separates") if self.plates.nodes[pa]["plate"]!=self.plates.nodes[pb]["plate"]])
        (x,y,z)=np.shape(segments)
        segments=np.reshape(segments,(x*y,z))
        lat,long=Projection().toLatLong(segments)
        lat=np.reshape(lat,(x,y))
        long=np.reshape(long,(x,y))
        lines=[[(slong[0],slat[0]),(slong[1],slat[1])] for slat,slong in zip(lat,long)]
        ax.add_collection(mc.LineCollection(lines,*args,transform=ccrs.Geodetic(),**kwargs))
    planet.plotPlatesBorders=types.MethodType(plotPlatesBorders,planet)
def distanceMatrixPlates(planet,planetGenerator):
    pass

def _solveKnapsack(values,weights,capacity):
    solver = pywrapknapsack_solver.KnapsackSolver(
        pywrapknapsack_solver.KnapsackSolver.
        KNAPSACK_MULTIDIMENSION_BRANCH_AND_BOUND_SOLVER, '')

    capacities = [capacity]

    solver.Init(values, [weights], capacities)
    computed_value = solver.Solve()

    packed_items = []
    packed_weights = []
    total_weight = 0
    for i in range(len(values)):
        if solver.BestSolutionContains(i):
            packed_items.append(i)
    return packed_items

def platesElevation(planet,planetGenerator):
    """select the plates to be ocean or land such that the immerged fraction is the maximal value smaller than the flood parameter.
    This is done by solving a knapsack problem."""
    size = [plate.number_of_nodes() for plate in planet.platesSubgraphs]
    weight=size
    capacity=int(planetGenerator.flood*sum(size))
    immergedPlates=_solveKnapsack(size,weight,capacity)    
    for i,p in enumerate(planet.platesSubgraphs):
        if(i in immergedPlates):
            p.elevation=planetGenerator.random.normal(loc=-7.5,scale=0.5)
        else:
            p.elevation=planetGenerator.random.normal(loc=1.5,scale=0.2)
    for i in planet.plates.nodes:
        n=planet.plates.nodes[i]
        n["elevation"]=planet.platesSubgraphs[n["plate"]].elevation
    def plotPlatesTypes(self,ax,*args,**kwargs):
        polygons=[]
        colors=[]
        for node,corners in self.meshCenters.nodes.data("corners"):
            points=np.array([self.meshCorners.nodes[k]["position"] for k in corners])    
            lat,long=Projection().toLatLong(points,"deg")            
            polygon=list(zip(long,lat))
            polygons.append(polygon)
            if self.plates.nodes[node]["elevation"]<0:
                colors.append((0,0,1))
            else:
                colors.append((0,1,0))
        ax.add_collection(mc.PolyCollection(polygons,*args,transform=ccrs.Geodetic(),facecolors=colors,**kwargs))
    planet.plotPlatesTypes=types.MethodType(plotPlatesTypes,planet)
    
def platesDistanceMatrix(planet,planetGenerator):
    for plate in planet.platesSubgraphs:
        for (source,targets) in nx.algorithms.shortest_paths.generic.shortest_path_length(plate,weight="length"):
            plate.nodes[source]["distances"]=targets