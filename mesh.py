import numpy as np
from scipy.spatial import SphericalVoronoi
import networkx as nx

def lengthSpherical(a,b):
    return np.arccos(np.sum(a*b,axis=-1))

class Tiling:
    """Generates tiling on the sphere.
    
    The tiling is composed of a set of unit vectors randomly placed on the sphere and there closest area definig a Voronoi map.
    """
    
    def __init__(self,seed=0):
        self._rng=np.random.RandomState(seed)
    def _samplePoints(self,N):
        """Generates N points on the sphere according to a quasirandom sequence. These generated points are already well distributed on the sphere.
        
        See http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/ for more info.
        """
        phi=1.324717957244746025960908854
        offset=self._rng.randint(2**32)
        u=np.modf((offset+np.arange(1,N+1))/phi)[0]
        v=np.modf((offset+np.arange(1,N+1))/phi**2)[0]
        l=np.arccos(2*u-1)+np.pi/2
        p=2*np.pi*v
        points=np.array([np.cos(l)*np.cos(p),np.cos(l)*np.sin(p),np.sin(l)]).transpose()
        return points
    def _relaxPoints(self,points,iterate=5):
        """Relaxes a set of point to have a more organic and uniform distribution.
        
        This is done by doing a couple tens of iteration of Lloyd relaxation by replacing the Voronoi cell points by the baricenter of the corners of the Voronoi cells.
        """
        for i in range(iterate):
            sv=SphericalVoronoi(points)
            for i in range(len(sv.points)):
                points[i]=np.mean(sv.vertices[sv.regions[i]],axis=0)
                points[i]/=np.linalg.norm(points[i])
        return points
    def generate(self,N,iterate=5):
        points=self._samplePoints(N)
        relaxedPoints=self._relaxPoints(points,iterate)
        sv=SphericalVoronoi(relaxedPoints)
        sv.sort_vertices_of_regions()
        return sv
    def createGraph(self,N,iterate=5):
        sv=self.generate(N,iterate)
        
        centers=nx.Graph()
        centers.add_nodes_from(list(zip(np.arange(len(sv.regions)),[{"corners":p,"center":center,"area":A} for p,center,A in zip(sv.regions,sv.points,sv.calculate_areas())])))
        
        corners=nx.Graph()
        corners.add_nodes_from(list(zip(np.arange(len(sv.vertices)),[{"position":p,"touches":[]} for p in sv.vertices])))
        
        for region in sv.regions:
            n=len(region)
            cornerLinks=[(region[i],region[(i+1)%n]) for i in range(len(region))]
            corners.add_edges_from(cornerLinks)
            
        for i in range(centers.number_of_nodes()):
            for corner in centers.nodes[i]["corners"]:
                corners.nodes[corner]["touches"].append(i)
        for region in sv.regions:
            n=len(region)
            cornerLinks=[(region[i],region[(i+1)%n]) for i in range(len(region))]
            corners.add_edges_from(cornerLinks)
            
        for (i,j) in corners.edges:
            corners.edges[i,j]["length"]=lengthSpherical(corners.nodes[i]["position"],corners.nodes[j]["position"])
            
        for i in range(corners.number_of_nodes()):
            l=corners.nodes[i]["touches"]
            n=len(l)
            centerLinks=[(l[j],l[(j+1)%n]) for j in range(n)]
            centers.add_edges_from(centerLinks)            
        return centers,corners
        
        
        
        