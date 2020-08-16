import numpy as np
class Illumination():
    def __init__(self,day=1,year=365,inclination=0):
        self.dayFreq=2*np.pi/day
        self.yearFreq=2*np.pi/year
        self.inclination=inclination
    def light(self,lon,lat,t):
        lat2d,t2d=np.meshgrid(lat,t,indexing="ij")
        lon2d,t2d=np.meshgrid(lon,t,indexing="ij")
        normal=np.array([np.cos(lat2d)*np.cos(lon2d),np.cos(lat2d)*np.sin(lon2d),np.sin(lat2d)])
        u=np.array([-np.cos(self.yearFreq*t2d)*np.cos(self.dayFreq*t2d)-np.sin(self.yearFreq*t2d)*np.sin(self.dayFreq*t2d)*np.cos(self.inclination),
            +np.cos(self.yearFreq*t2d)*np.sin(self.dayFreq*t2d)-np.sin(self.yearFreq*t2d)*np.cos(self.dayFreq*t2d)*np.cos(self.inclination),
            -np.sin(self.yearFreq*t2d)*np.sin(self.inclination)])
        lightMap=-(u[0]*normal[0]+u[1]*normal[1]+u[2]*normal[2])
        lightMap=np.maximum(0,lightMap)
        return lightMap