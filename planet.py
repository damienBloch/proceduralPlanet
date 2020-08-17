import numpy as np
from dataclasses import dataclass

@dataclass
class PlanetParameters:
    radius: float
    starTemperature: float
    starRadius: float
    starDistance: float
        
class Planet:
    def __init__(self,parameters):
        self.parameters=parameters
        