import numpy as np
import sys
import types
from utils import Projection
from matplotlib import collections as mc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from ortools.algorithms import pywrapknapsack_solver
from graph_tool.all import *
from scipy.linalg import pinv

def seedsPlates(planet,planetGenerator):
    planet.plates=Graph()
    planet.plates.add_vertex(planet.meshCenters.num_vertices())
    seeds=planetGenerator.random.choice(np.arange(planet.plates.num_vertices()),size=planetGenerator.numberPlates,replace=False)
    
    planet.plates.vertex_properties["color"] = planet.plates.new_vertex_property("int")
    for i in planet.plates.vertices():
        planet.plates.vp.color[i]=-1
    for i in range(planetGenerator.numberPlates):
        planet.plates.vp.color[seeds[i]]=i
        
def colorPlates(planet,planetGenerator):
    random=planetGenerator.random
    queue=[i for i,p in enumerate(planet.plates.vp.color.a) if p>=0]
    while(len(queue)>0):
        index=random.choice(queue)
        uncoloredNeighbors=[n for n in planet.meshCenters.vertex(index).out_neighbors() if planet.plates.vp.color[n]<0]
        if(len(uncoloredNeighbors)>0):
            newIndex=random.choice(uncoloredNeighbors)
            planet.plates.vp.color[newIndex]=planet.plates.vp.color[index]
            queue.append(newIndex)
        else:
            queue.remove(index)
    
    def plotPlatesBorders(self,ax,*args,**kwargs):
        edges=[[int(e.source()),int(e.target())] for e in self.meshCorners.edges() 
       if self.plates.vp.color[self.meshCorners.ep.separates[e][0]]
       != self.plates.vp.color[self.meshCorners.ep.separates[e][1]]]
        segments=np.array(list(self.meshCorners.vp.position))[edges]
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
    size=[]
    for k in range(planetGenerator.numberPlates):
        plate = GraphView(planet.meshCenters,planet.plates.vp.color.a==k)
        size.append(plate.num_vertices())
    weight=size
    capacity=int(planetGenerator.flood*sum(size))
    immergedPlates=_solveKnapsack(size,weight,capacity)
    planet.meshCenters.vertex_properties["elevation"] = planet.meshCenters.new_vertex_property("float")
    
    def plateElevation(oceanic):
        if oceanic:
            return planetGenerator.random.normal(loc=-7.5,scale=0.5)
        else:
            return planetGenerator.random.normal(loc=1.5,scale=0.2)
    elevation=np.array([plateElevation(p in immergedPlates) for p in range(planetGenerator.numberPlates)])
    planet.meshCenters.vp.elevation.a=elevation[planet.plates.vp.color.a]

    def plotPlatesTypes(self,ax,*args,**kwargs):
        positions=np.array(list(self.meshCorners.vp.position))
        polygons=[]
        for k,corners in enumerate(self.meshCenters.vp.corners):
            lat,long=Projection().toLatLong(positions[corners])
            polygon=list(zip(long,lat))
            polygons.append(polygon)
        elevations=np.array(list(self.meshCenters.vp.elevation))
        r=np.piecewise(elevations,[elevations<0],[lambda x:0,lambda x:0])
        g=np.piecewise(elevations,[elevations<0],[lambda x:0,lambda x:1])
        b=np.piecewise(elevations,[elevations<0],[lambda x:0.8,lambda x:0])
        colors=np.array([r,g,b]).T
        ax.add_collection(mc.PolyCollection(polygons,*args,transform=ccrs.Geodetic(),facecolors=colors,**kwargs))
    planet.plotPlatesTypes=types.MethodType(plotPlatesTypes,planet)
    
def platesResistanceMatrix(planet,planetGenerator):
    planet.plates.vertex_properties["resistance"] = planet.plates.new_vertex_property("object")
    for k in range(planetGenerator.numberPlates):
        u=GraphView(planet.meshCenters,vfilt=planet.plates.vp.color.a==k)
        l=u.ep.length.a
        conductance=1/l
        cprop=u.new_edge_property("float")
        cprop.a=conductance
        
        L=pinv(laplacian(u,weight=cprop).todense())
        P=np.full_like(L,planet.parameters.radius/u.num_vertices())
        G=L+P
        diag=G[np.arange(len(G)),np.arange(len(G))]
        d1,d2=np.meshgrid(diag,diag)
        R=d1+d2-2*G
        
        
        
        d=dict(zip(u.get_vertices(),np.arange(u.num_vertices())))
        for n in u.vertices():
            planet.plates.vp.resistance[int(n)]={k:R[d[int(n)],d[k]] for k in d}
def relief(planet,planetGenerator):
    def randomSpeed():
            platesRotationAxes=planetGenerator.random.normal(size=(planetGenerator.numberPlates,3))
            norm=np.linalg.norm(platesRotationAxes,axis=-1)
            platesRotationAxes/=norm[:,None]
            platesRotationSpeed=planetGenerator.random.rand(planetGenerator.numberPlates)*200-100
            position=np.array(list(planet.meshCenters.vp.center))
            plate=np.array(list(planet.plates.vp.color))
            speeds=np.cross(position,platesRotationAxes[plate]*platesRotationSpeed[plate,None])
            planet.plates.vertex_properties["speed"] = planet.plates.new_vertex_property("vector<float>")
            for i,n in enumerate(planet.plates.vertices()):
                planet.plates.vp.speed[n]=speeds[i]
    randomSpeed()
    def plotPlatesSpeed(self,ax,*args,**kwargs):
        dt=2/100./np.sqrt(planet.meshCenters.num_vertices())
        p0=np.array(list(planet.meshCenters.vp.center))
        p1=p0+dt*np.array(list(planet.plates.vp.speed))
        lat0,long0=Projection().toLatLong(p0,"deg")
        lat1,long1=Projection().toLatLong(p1,"deg")
        lines=[[(long0[i],lat0[i]),(long1[i],lat1[i])] for i,_ in enumerate(lat0)]
        ax.add_collection(mc.LineCollection(lines,*args,transform=ccrs.Geodetic(),**kwargs))
    planet.plotPlatesSpeed=types.MethodType(plotPlatesSpeed,planet)
    def computePressure():
        edges=[[edge.source(),edge.target(),planet.meshCorners.ep.separates[edge][0],planet.meshCorners.ep.separates[edge][1],planet.meshCorners.ep.length[edge]] for edge in planet.meshCorners.edges() 
       if planet.plates.vp.color[planet.meshCorners.ep.separates[edge][0]]
       != planet.plates.vp.color[planet.meshCorners.ep.separates[edge][1]]]
        planet.plates.vertex_properties["pressure"] = planet.plates.new_vertex_property("float")
        for corner_i,corner_j,tile_a,tile_b,length in edges:
            normal=np.array(planet.meshCenters.vp.center[tile_a])-np.array(planet.meshCenters.vp.center[tile_b])
            normal/=np.linalg.norm(normal)
            vrel=np.array(planet.plates.vp.speed[tile_a])-np.array(planet.plates.vp.speed[tile_b])
            vrelNorm=np.dot(normal,vrel)
            pressure=-vrelNorm/200*length
            planet.plates.vp.pressure[tile_a]+=pressure
            planet.plates.vp.pressure[tile_b]+=pressure
    computePressure()
    def plotPressure(self,ax,*args,**kwargs):
        lat,long=Projection().toLatLong(self.meshPoints,"deg")
        color=list(map(lambda p: [1,0,0] if p>0 else [0.1,0.1,1],self.plates.vp.pressure.a))
        ax.scatter(long,lat,*args,s=np.abs(self.plates.vp.pressure.a),transform=ccrs.PlateCarree(),c=color,**kwargs)      
    planet.plotPressure=types.MethodType(plotPressure,planet)
    def computeRelief():
        edges=[[edge.source(),edge.target(),planet.meshCorners.ep.separates[edge][0],planet.meshCorners.ep.separates[edge][1],planet.meshCorners.ep.length[edge]] for edge in planet.meshCorners.edges() 
       if planet.plates.vp.color[planet.meshCorners.ep.separates[edge][0]]
       != planet.plates.vp.color[planet.meshCorners.ep.separates[edge][1]]]
        for corner_i,corner_j,tile_a,tile_b,length in edges:
            #continental interaction
            if(planet.meshCenters.vp.elevation[tile_a]>0 and planet.meshCenters.vp.elevation[tile_b]>0):                
                #continental collision, creates mountain range
                addedElevation=0
                if(planet.plates.vp.pressure[tile_a]>0):
                    addedElevation+=planet.plates.vp.pressure[tile_a]/planet.meshCenters.vp.area[tile_a]**.5*1.5
                if(planet.plates.vp.pressure[tile_b]>0):
                    addedElevation+=planet.plates.vp.pressure[tile_b]/planet.meshCenters.vp.area[tile_b]**.5*1.5
                if(planet.plates.vp.pressure[tile_a]>0):
                    for onSamePlate,distance in planet.plates.vp.resistance[tile_a].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=np.exp(-(distance)/125)*addedElevation
                if(planet.plates.vp.pressure[tile_b]>0):
                    for onSamePlate,distance in planet.plates.vp.resistance[tile_b].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=np.exp(-(distance)/125)*addedElevation
                        
                #continental rift
                addedElevation=0
                if(planet.plates.vp.pressure[tile_a]<0):
                    addedElevation+=planet.plates.vp.pressure[tile_a]/planet.meshCenters.vp.area[tile_a]**.5*1.5
                if(planet.plates.vp.pressure[tile_b]<0):
                    addedElevation+=planet.plates.vp.pressure[tile_b]/planet.meshCenters.vp.area[tile_b]**.5*1.5   
                if(planet.plates.vp.pressure[tile_a]<0):
                    for onSamePlate,distance in planet.plates.vp.resistance[tile_a].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=1/(1+(distance/200)**5)*addedElevation
                if(planet.plates.vp.pressure[tile_b]<0):
                    for onSamePlate,distance in planet.plates.vp.resistance[tile_b].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=1/(1+(distance/200)**5)*addedElevation
            #oceanic/continental interaction
            #subduction
            def subduction(continent,ocean):
                pressureContinent=planet.plates.vp.pressure[continent]
                areaContinent=planet.meshCenters.vp.area[continent]
                if(pressureContinent>0):
                    addedElevation=1.5*pressureContinent/areaContinent**.5
                    for onSamePlate,distance in planet.plates.vp.resistance[continent].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=np.exp(-(distance/100)**2/2)*addedElevation
                pressureOcean=planet.plates.vp.pressure[ocean]
                areaOcean=planet.meshCenters.vp.area[ocean]
                if(pressureOcean>0):
                    addedElevation=-pressureOcean/areaOcean**.5
                    for onSamePlate,distance in planet.plates.vp.resistance[ocean].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=1/(1+distance/600)*addedElevation
                    addedElevation=(planet.meshCenters.vp.elevation[continent]-planet.meshCenters.vp.elevation[ocean])*.8
                    for onSamePlate,distance in planet.plates.vp.resistance[ocean].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=1/(1+(distance/100)**3)*addedElevation
            if(planet.meshCenters.vp.elevation[tile_a]>0 and planet.meshCenters.vp.elevation[tile_b]<0):
                subduction(tile_a,tile_b)
            if(planet.meshCenters.vp.elevation[tile_b]>0 and planet.meshCenters.vp.elevation[tile_a]<0):
                subduction(tile_b,tile_a)
            #emergence
            def emergence(continent,ocean):
                pressureContinent=planet.plates.vp.pressure[continent]
                areaContinent=planet.meshCenters.vp.area[continent]
                if(pressureContinent<0):
                    addedElevation=-.1*pressureContinent/areaContinent**.5
                    for onSamePlate,distance in planet.plates.vp.resistance[continent].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=np.exp(-(distance/100)**2/2)*addedElevation
                pressureOcean=planet.plates.vp.pressure[ocean]
                areaOcean=planet.meshCenters.vp.area[ocean]
                if(pressureOcean<0):
                    addedElevation=pressureOcean/areaOcean**.5*.2
                    for onSamePlate,distance in planet.plates.vp.resistance[ocean].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=1/(1+distance/600)*addedElevation
                    addedElevation=(planet.meshCenters.vp.elevation[continent]-planet.meshCenters.vp.elevation[ocean])*.2
                    for onSamePlate,distance in planet.plates.vp.resistance[ocean].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=1/(1+(distance/100)**4)*addedElevation*2
            if(planet.meshCenters.vp.elevation[tile_a]>0 and planet.meshCenters.vp.elevation[tile_b]<0):
                emergence(tile_a,tile_b)
            if(planet.meshCenters.vp.elevation[tile_b]>0 and planet.meshCenters.vp.elevation[tile_a]<0):
                emergence(tile_b,tile_a)
            #oceanic interaction
            
            elevationA=planet.meshCenters.vp.elevation[tile_a]
            elevationB=planet.meshCenters.vp.elevation[tile_b]
            if(elevationA<0 and elevationB<0):
                addedElevation=1
                pressureA=planet.plates.vp.pressure[tile_a]
                pressureB=planet.plates.vp.pressure[tile_b]
                #dorsale
                if(pressureA<0 and pressureB<0):
                    elevationDifference=(elevationA-elevationB)*.5
                    areaA=planet.meshCenters.vp.area[tile_a]
                    areaB=planet.meshCenters.vp.area[tile_b]
                    addedElevation-=pressureA/areaA**.5*.15
                    addedElevation-=pressureB/areaB**.5*.15
                    for onSamePlate,distance in planet.plates.vp.resistance[tile_a].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=1/(1+(distance/100)**2)*(addedElevation-elevationDifference)
                    for onSamePlate,distance in planet.plates.vp.resistance[tile_b].items():
                        planet.meshCenters.vp.elevation[onSamePlate]+=1/(1+(distance/100)**2)*(addedElevation+elevationDifference)
            
    computeRelief()
    def plotRelief(self,ax,*args,**kwargs):
        positions=np.array(list(self.meshCorners.vp.position))
        polygons=[]
        for k,corners in enumerate(self.meshCenters.vp.corners):
            lat,long=Projection().toLatLong(positions[corners])
            polygon=list(zip(long,lat))
            polygons.append(polygon)
        elevations=np.array(list(self.meshCenters.vp.elevation))
        
        cmap_terrain = plt.cm.get_cmap('gist_earth')
        cmap_ocean = plt.cm.get_cmap('ocean')
        def color(e):
            if e>0:
                y=np.interp(e,[0,10],[0.5,1])
                return cmap_terrain(y)
            else:
                y=np.interp(e,[-10,0],[0.25,0.8])
                return cmap_ocean(y)
        colors=list(map(color,elevations))
        ax.add_collection(mc.PolyCollection(polygons,*args,transform=ccrs.Geodetic(),facecolors=colors,**kwargs))
    planet.plotRelief=types.MethodType(plotRelief,planet)