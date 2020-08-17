import numpy as np
from utils import Projection,pairwise,pairwiseRoll,lengthSpherical
import cartopy.crs as ccrs
import types
from scipy.spatial import SphericalVoronoi
import networkx as nx

def meshSeeds(planet,planetGenerator):
    
    phi=1.324717957244746025960908854
    offset=planetGenerator.random.randint(2**32)
    
    N=planetGenerator.numberTiles
    u=np.modf((offset+np.arange(1,N+1))/phi)[0]
    v=np.modf((offset+np.arange(1,N+1))/phi**2)[0]
    l=np.arccos(2*u-1)+np.pi/2
    p=2*np.pi*v
    points=np.array([np.cos(l)*np.cos(p),np.cos(l)*np.sin(p),np.sin(l)]).transpose()
    planet.meshPoints=points
    
    def plotMeshPoints(self,ax,*args,**kwargs):
        lat,long=Projection().toLatLong(self.meshPoints,"deg")
        ax.scatter(long,lat,*args,transform=ccrs.PlateCarree(),**kwargs)
    planet.plotMeshPoints=types.MethodType(plotMeshPoints,planet)
    
def _relaxPoints(points,stagnation=1e-2):
    """Relaxes a set of point to have a more organic and uniform distribution.
    This is done by doing a couple of iteration of Lloyd relaxation by replacing the Voronoi cell points by the baricenter of the corners of the Voronoi cells.
    """
    while True:
        sv=SphericalVoronoi(points)
        flatRegions = [item for sublist in sv.regions for item in sublist]
        groupSize=[len(sublist) for sublist in sv.regions]
        startIndex=np.zeros(len(groupSize)+1,dtype=int)
        startIndex[1:]=np.cumsum(groupSize)
        barycenters=[np.mean(sv.vertices[flatRegions[i:j]],axis=0) for (i,j) in pairwise(startIndex)]
        newPoints=barycenters/np.linalg.norm(barycenters,axis=-1)[:,None]
        if(np.max(np.linalg.norm(newPoints-points,axis=-1))<stagnation):
            break
        points=newPoints
    return newPoints
    
def relaxMesh(planet,planetGenerator):
    planet.meshPoints=_relaxPoints(planet.meshPoints)
    
def createGraph(planet,planetGenerator):
    sv=SphericalVoronoi(planet.meshPoints)

    centers=nx.Graph()
    centers.add_nodes_from(list(zip(np.arange(len(sv.regions)),[{"corners":p,"center":center,"area":A} for p,center,A in zip(sv.regions,sv.points,sv.calculate_areas())])))

    for k in centers.nodes:
        centers.nodes[k]["label"]=k

    corners=nx.Graph()
    corners.add_nodes_from(zip([i for i,_ in enumerate(sv.vertices)],[{"position":p,"touches":[]} for p in sv.vertices]))

    #ensure the the polygon is clockwise oriented and correct it if not
    for i,(region,point) in enumerate(zip(sv.regions,sv.points)):
        if(np.dot(np.cross(point,sv.vertices[region[0]]),sv.vertices[region[1]])>0):
            sv.regions[i].reverse()

    for region in sv.regions:
        n=len(region)
        cornerLinks=[(region[i],region[(i+1)%n]) for i in range(len(region))]
        corners.add_edges_from(cornerLinks)

    for i in range(centers.number_of_nodes()):
        for corner in centers.nodes[i]["corners"]:
            corners.nodes[corner]["touches"].append(i)
    for region in sv.regions:
        n=len(region)
        cornerLinks=[(i,j) for i,j in pairwiseRoll(region)]
        corners.add_edges_from(cornerLinks)

    #compute length of each border 
    for (i,j) in corners.edges:
        corners.edges[i,j]["length"]=lengthSpherical(corners.nodes[i]["position"],corners.nodes[j]["position"])
        corners.edges[i,j]["separates"]=set(corners.nodes[i]["touches"]).intersection(corners.nodes[j]["touches"])

    for i in range(corners.number_of_nodes()):
        l=corners.nodes[i]["touches"]
        n=len(l)
        centerLinks=[(l[j],l[(j+1)%n]) for j in range(n)]
        centers.add_edges_from(centerLinks) 
    for (i,j) in centers.edges:
        p1=centers.nodes[i]["center"]
        p2=centers.nodes[j]["center"]
        centers.edges[i,j]["length"]=lengthSpherical(p1,p2)
    planet.meshCenters=centers
    planet.meshCorners=corners
    
    def plotTilesBorders(self,ax,*args,**kwargs):
        [[planet.meshCorners.nodes[i]["position"],planet.meshCorners.nodes[j]["position"]] for (i,j) in planet.meshCorners.edges]
        #lat,long=Projection().toLatLong(self.meshPoints,"deg")
        #ax.scatter(long,lat,*args,transform=ccrs.PlateCarree(),**kwargs)
    planet.plotTilesBorders=types.MethodType(plotTilesBorders,planet)
        
    