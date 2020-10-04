from planet import Planet,PlanetParameters
import numpy as np
from mesh_generator import meshSeeds,relaxMesh,createGraph
from plates_generator import seedsPlates,colorPlates,platesElevation,platesResistanceMatrix,relief,randomRelief
from copy import deepcopy
from timeit import timeit
from si_prefix import si_format
from warnings import warn

class PlanetGenerator:  
    
    callbacks=[meshSeeds,relaxMesh,createGraph,seedsPlates,colorPlates,platesElevation,platesResistanceMatrix,relief,randomRelief]
    callbacksDescription=["Generate mesh seeds","Relax mesh points","Create graph","Pick plates seeds","Color plates","Pick plates type","Compute resistance matrix within plates","Generate relief","Add random relief"]
    numberSteps=len(callbacks)    
    
    def __init__(self,parameters,numberTiles,seed=0,numberPlates=30,flood=0.6,saveIntermediates=False):
        self.planet=None
        self.construction=[{} for i in range(PlanetGenerator.numberSteps)]   
        self.parameters=parameters
        self.random=np.random.RandomState(seed)
        self.seed=seed
        self.numberTiles=numberTiles
        self.numberPlates=numberPlates
        self.flood=flood
        self._saveIntermediates=saveIntermediates
        
    def build(self,start=0,stop=None):
        
        if(stop is None):
                stop=PlanetGenerator.numberSteps
        
        for k in range(stop,PlanetGenerator.numberSteps):
            self.construction[k]={}            
        
        if start==0:
            self.planet=Planet(self.parameters)
            self.random=np.random.RandomState(self.seed)
        else:
            if(not self._saveIntermediates):
                warn("Intermediate states were not saved. Building the structure from scratch.")
                self.random=np.random.RandomState(self.seed)
                self.build(start=0,stop=stop)     
                return
            else:
                self.planet=self.construction[start-1]["planet"]
                self.random=self.construction[start-1]["random"]
            
        for i,(callback,desc) in enumerate(zip(PlanetGenerator.callbacks[start:stop],PlanetGenerator.callbacksDescription[start:stop])):
            print(desc+"...",end="")
            duration=timeit(lambda:callback(self.planet,self),number=1)
            if(self._saveIntermediates):
                self.construction[i]["planet"]=deepcopy(self.planet)
                self.construction[i]["random"]=deepcopy(self.random)
            self.construction[i]["duration"]=duration
            print(" Done in {}s".format(si_format(duration)))
        #print("Total time: {}s".format(si_format(sum([d["duration"] for d in self.construction]))))