import numpy as np
from utils import Projection,pairwise,pairwiseRoll,lengthSpherical
import cartopy.crs as ccrs
import types
from scipy.spatial import SphericalVoronoi
from graph_tool.all import *
from matplotlib import collections as mc

def meshSeeds(planet,planetGenerator):
    
    phi=1.324717957244746025960908854
    offset=planetGenerator.random.randint(2**32)
    
    N=planetGenerator.numberTiles
    u=np.modf((offset+np.arange(1,N+1))/phi)[0]
    v=np.modf((offset+np.arange(1,N+1))/phi**2)[0]
    l=np.arccos(2*u-1)+np.pi/2
    p=2*np.pi*v
    points=np.array([np.cos(l)*np.cos(p),np.cos(l)*np.sin(p),np.sin(l)]).transpose()
    #add some random displacement otherwise tiles near the equator are aligned
    points+=planetGenerator.random.normal(loc=0,scale=1e-2,size=np.shape(points))
    points/=np.linalg.norm(points,axis=-1)[:,None]
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
    #the sphere is spanned into regions of closer points to the seeds with a spherical Voronoi mapping.
    sv=SphericalVoronoi(planet.meshPoints)
    sv.sort_vertices_of_regions()
    
    #ensure that the polygon is clockwise oriented and correct it if not. 
    #For plotting, all the polygons must have the same orientation. 
    for i,region in enumerate(sv.regions):
        if(np.dot(np.cross(sv.points[i],sv.vertices[region[0]]),sv.vertices[region[1]])>0):
            sv.regions[i].reverse()
    
    #the graph of regions centers is created. This is useful for checking travel between regions. 
    #At this point it is not yet possible to set the connections between regions.
    centers=Graph(directed=False)
    centers.add_vertex(len(sv.regions))
    
    #for each cell, register the corners points index delimiting it.
    centers.vertex_properties["corners"] = centers.new_vertex_property("vector<int>")
    for v,corners in zip(centers.vertices(),sv.regions):
        centers.vp.corners[v]=corners
        
    #for each cell, set its center point.
    centers.vertex_properties["center"] = centers.new_vertex_property("vector<float>")
    for v,center in zip(centers.vertices(),sv.points):
        centers.vp.center[v]=center
      
    #for each cell, set its area in km^2. This is often used as a weight for the cell.
    centers.vertex_properties["area"] = centers.new_vertex_property("float")
    for v,area in zip(centers.vertices(),sv.calculate_areas()):
        centers.vp.area[v]=area*planet.parameters.radius**2
        
    #create the graph of corners points. It is the dual of the graph centers. 
    #this graph is useful to plot the borders and check touching regions.
    corners=Graph(directed=False)
    corners.add_vertex(len(sv.vertices))
    
    #for each corner, register its position.
    corners.vertex_properties["position"] = corners.new_vertex_property("vector<float>")
    for v,position in zip(corners.vertices(),sv.vertices):
        corners.vp.position[v]=position
    

            
    #set the edges between corners
    cornerLinks=[]
    for region in sv.regions:
        cornerLinks+=[(min(i,j),max(i,j)) for (i,j) in pairwiseRoll(region)]
    #the edges are converted to dictionnary and back to a list to remove the doublons
    corners.add_edge_list(list(dict.fromkeys(cornerLinks)))
    
    #set the polygons that each corner touches
    corners.vertex_properties["touches"] = corners.new_vertex_property("vector<int>")
    for i,cornerList in enumerate(centers.vp.corners):
        for corner in cornerList:
            corners.vp.touches[corner].append(i)
            
    #compute length of each border in km
    corners.edge_properties["length"] = corners.new_edge_property("float")
    edges=corners.get_edges()
    positions=np.array(list(corners.vp.position))[edges]
    lengths=np.arccos(np.sum(positions[:,0]*positions[:,1],axis=-1))*planet.parameters.radius    
    for i,e in enumerate(corners.edges()):
        corners.ep.length[e]=lengths[i]
    corners.edge_properties["separates"] = corners.new_edge_property("vector<int>")
    for e in corners.edges(): 
        corners.ep.separates[e]=list(set(corners.vp.touches[e.source()]).intersection(corners.vp.touches[e.target()]))
        
    #add links between regions
    centerLinks=[]
    for c in corners.vertices():
        l=corners.vp.touches[c]
        centerLinks+=[(min(i,j),max(i,j)) for (i,j) in pairwiseRoll(l)]
    centers.add_edge_list(list(dict.fromkeys(centerLinks))) 
        
    #compute distance between regions in km
    centers.edge_properties["length"] = centers.new_edge_property("float")
    edges=centers.get_edges()
    positions=np.array(list(centers.vp.center))[edges]
    lengths=np.arccos(np.sum(positions[:,0]*positions[:,1],axis=-1))*planet.parameters.radius    
    for i,e in enumerate(centers.edges()):
        centers.ep.length[e]=lengths[i]       

    planet.meshCenters=centers    
    planet.meshCorners=corners    
    
    def plotTilesBorders(self,ax,*args,**kwargs):
        edges=self.meshCorners.get_edges()
        segments=np.array(list(self.meshCorners.vp.position))[edges]
        (x,y,z)=np.shape(segments)
        segments=np.reshape(segments,(x*y,z))
        lat,long=Projection().toLatLong(segments)
        lat=np.reshape(lat,(x,y))
        long=np.reshape(long,(x,y))
        lines=[[(slong[0],slat[0]),(slong[1],slat[1])] for slat,slong in zip(lat,long)]
        ax.add_collection(mc.LineCollection(lines,*args,transform=ccrs.Geodetic(),**kwargs))
    planet.plotTilesBorders=types.MethodType(plotTilesBorders,planet)
    
    def plotTilesJunctions(self,ax,*args,**kwargs):
        edges=self.meshCenters.get_edges()
        segments=np.array(list(self.meshCenters.vp.center))[edges]
        (x,y,z)=np.shape(segments)
        segments=np.reshape(segments,(x*y,z))
        lat,long=Projection().toLatLong(segments)
        lat=np.reshape(lat,(x,y))
        long=np.reshape(long,(x,y))
        lines=[[(slong[0],slat[0]),(slong[1],slat[1])] for slat,slong in zip(lat,long)]
        ax.add_collection(mc.LineCollection(lines,*args,transform=ccrs.Geodetic(),**kwargs))
    planet.plotTilesJunctions=types.MethodType(plotTilesJunctions,planet)   