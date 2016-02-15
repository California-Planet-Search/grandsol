from __future__ import division
import numpy as np
from scipy.constants import c

def z2v(z):
    """
    Convert redshift to velocity using relativistic equations

    Parameters:
    z : float or array-like : redshift (V/c)

    Returns:
    v : float or array-like : velocity in m/s
    """

    v = c * ((1 + z)**2 - 1) / ((1 + z)**2 +1)

    return v

def redshift_to_vel(z, uz=0):
    """
    Convert redshift and uncertainty to velocity and uncertainty using relativistic equations

    Parameters:
    z : float or array-like : redshift (V/c)
    uz : float or array-like : uncertainty on redshift (V/c)

    Returns:
    (v, uv) : tuple : velocity and uncertainty in m/s
    """
    v = z2v(z)
    
    uv = np.sqrt(uz**2 * (4 * (1+z)) / ((1+z)**2 + 1)**2)

    return (v, uv)

class RV():
    """
    Radial Velocity object designed to handle relativistic arithmatic.
    Can be constructed by giving velocities or Z. Velocities should be in m/s

    Example: V = RV(vel=mnvel)  or   V = RV(z=zn)

    Attributes
    ----------
    vel : float or array-like : radial velocities

    Methods
    -------
    __add__, __sub__, __mul__, __div__ : handle the +, -, *, and / operators respectively and return new RV objects
    with the relativistic sum or difference of the original RV.vel objects

    sum : sums an array of velocities over the specified axis using relativistic addition

    """
    
    def __init__(self, **kwargs):#, vel=None, z=None):
        kwkeys = kwargs.keys()
        if 'z' in kwkeys: self.vel = z2v(kwargs['z'])
        if 'vel' in kwkeys: self.vel = kwargs['vel']

    def __repr__(self):
        return "<RV Object>\n%s" % repr(self.vel).replace("<RV Object>\n","")

    def __neg__(self):
        return RV( vel=-self.vel)
    
    def __add__(self, u):
        v = self.vel
        u = u.vel
        return RV( vel = (u + v) / (1 + (u*v/(c**2))) )

    def __sub__(self, u):
        return self.__add__(-u)

    def __mul__(self, N):
        "From Talukder & Ahmad (2012), published in AJASE"
        
        if not isinstance(N, float) and not isinstance(N, int):
            raise TypeError, "RV objects can only be multiplied or divided by scalars."
        else:
            g = self.vel/c
            m = c * ((1 + g)**N - (1 - g)**N) / ((1 + g)**N + (1 - g)**N)
            return RV( vel = m )
        
    def __div__(self, N):
        return self.__mul__(1/N)

    def sum(self, axis=0):
        l = self.vel.shape[axis]
        out = np.zeros_like(np.sum(self.vel, axis=axis))
        for i in range(l):
            out = (RV( vel=out ) + RV( vel=self.vel[i,:] )).vel

        return RV( vel=out )

    def mean(self, axis=0):
        N = self.vel.shape[axis]
        s = self.sum(axis=axis)
        g = s.vel/c
        
        mnvel = s.__div__(N)

        return RV( vel=mnvel )
    
    def values(self):
        return np.array(self.vel)
