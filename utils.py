import numpy as np
from numba import njit
import itertools

@njit
def toLatLong(r,ex,ey,ez):
    lat=np.zeros(r.shape[0])
    long=np.zeros(r.shape[0])
    for i in range(r.shape[0]):
        lat[i]=np.arccos(r[i][2]/np.sqrt(r[i][0]**2+r[i][1]**2+r[i][2]**2))
        long[i]=np.arctan2(r[i][1],r[i][0])
    return lat,long

class Projection:
        def __init__(self):
            self.ex=np.array([1,0,0])
            self.ey=np.array([0,1,0])
            self.ez=np.array([0,0,1])
        def toLatLong(self,r,units="deg"):
            lat,long=toLatLong(r,self.ex,self.ey,self.ez)
            if(units=="rad"):
                return lat,long
            elif(units=="deg"):
                return (lat-np.pi/2)/np.pi*180,long/np.pi*180            

def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b) 

def lengthSpherical(a,b):
    return np.arccos(np.sum(a*b,axis=-1))

def pairwiseRoll(border):
    n=len(border)
    for i in range(n):
        yield((border[i],border[(i+1)%n]))
    