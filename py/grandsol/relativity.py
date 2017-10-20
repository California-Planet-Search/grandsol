from __future__ import division
import numpy as np
from scipy.constants import c

def z2v(z):
    """
    Convert redshift to velocity using equation :eq:`z2v`


    Args:
        z (float or array-like): redshift (V/c)

    Returns:
        float or array-like: velocity in m/s
    """

    v = c * ((1 + z)**2 - 1) / ((1 + z)**2 +1)

    return v

def redshift_to_vel(z, uz=0):
    """
    Convert redshift and uncertainty to velocity and uncertainty using relativistic equations.

    .. math::
        :label: z2v
        
        v = c\\frac{(1 + z)^2 - 1)}{(1 + z)^2 +1}
        
        \\delta v = \\sqrt{\\delta z^2\\frac{4(1 + z)}{((1 + z)^2 + 1)^2}}

    Args:
        z (float or array-like): redshift (V/c)
        uz (float or array-like): uncertainty on redshift (V/c)

    Returns:
        tuple : (velocity,  uncertainty in m/s)
    """
    v = z2v(z)
    
    uv = np.sqrt(uz**2 * (4 * (1+z)) / ((1+z)**2 + 1)**2)

    return (v, uv)

class RV():
    """
    Radial Velocity object designed to handle relativistic arithmatic.
    Can be constructed by giving velocities or Z. Velocities should be in m/s

    Args:
        vel (float or array): velocity or array of velocities
        z (float or array): Z or array of Z where :math:`Z = v/c`

    Note:
        * ``__add__``, ``__sub__``, ``__mul__``, and ``__div__`` handle the +, -, *,
        and / operators respectively and return new RV objects with the relativistic
        arithmatic performed on the original RV.vel objects
        * Cannot define `both` vel and z when constructing the RV object
    
    Attributes:
        vel (float or array-like): attribute of an RV object containing the radial velocities
        

    Examples:
        >>> V1 = RV(vel=mnvel)   # construct RV object by giving the velocities
        >>> V2 = RV(z=zn)        # construct RV object by giving Z
        >>> V3 = V1 + V2         # perform relativistic addition

    """
    
    def __init__(self, vel=None, z=None):
        if vel is not None: self.vel = vel
        if z is not None: self.vel = z2v(z)  
        if (vel is None) and (z is None): raise ValueError, \
            "vel or z must be specified when constructing RV object"
        if (vel is not None) and (z is not None): raise ValueError, \
            "cannot specifiy both vel and z when constructing RV object"

    def __repr__(self):
        return "<RV Object>\n%s" % repr(self.vel).replace("<RV Object>\n","")

    def __neg__(self):
        return RV( vel=-self.vel)
    
    def __add__(self, u):
        """
        Relativistic addition of two RV objects.

        .. math:: u' = \\frac{v + u}{1 + \\frac{vu}{c^2}}
        is the relativistic equivalent of :math:`u' = u + v`
        """
        
        v = self.vel / c
        u = u.vel / c
        new = (u+v) / (1 + u*v)
        return RV( vel = (new * c) )

    def __sub__(self, u):
        """
        Relativistic subtraction of two RV objects.
        """

        return self.__add__(-u)

    def __mul__(self, N):
        """
        Relativistic multiplication of an RV with a scalar.
        From Talukder & Ahmad (2012), published in AJASE.

        .. math:: u' = c\\frac{(1 + z)^N - (1 - z)^N}{(1 + z)^N + (1 - g)^N}
        where :math:`z=v/c` is the relativistic equivalent of :math:`u' = Nv` where :math:`N` is a scalar.

        Note:
            Multiplication of two RV objects is not currently implemented.

        """
        
        if not isinstance(N, float) and not isinstance(N, int):
            raise TypeError, "RV objects can only be multiplied or divided by scalars."
        else:
            g = self.vel/c
            m = c * ((1 + g)**N - (1 - g)**N) / ((1 + g)**N + (1 - g)**N)
            return RV( vel = m )
        
    def __div__(self, N):
        """
        Relativistic division of an RV with a scalar.

        Note:
            Division of two RV objects is not currently implemented.

        """
        
        return self.__mul__(1/N)

    def sum(self, axis=0):
        """
        Calculate the relativistic sum of the velocities contained
        in an RV object along the specified axis.

        Args:
            axis (int): perform sum along this axis

        Returns:
            grandsol.relativity.RV: new RV object colapsed by summing over the specified axis
        """
        
        l = self.vel.shape[axis]
        out = np.zeros_like(np.sum(self.vel, axis=axis))
        for i in range(l):
            out = (RV( vel=out ) + RV( vel=self.vel[i,:] )).vel

        return RV( vel=out )

    def mean(self, axis=0):
        """
        Calculate the relativistic mean of velocities contained in
        an RV object along the specified axis.

        Args:
            axis (int): perform sum along this axis

        Returns:
            grandsol.relativity.RV: new RV object colapsed by averaging over the specified axis
        """

        
        N = self.vel.shape[axis]
        s = self.sum(axis=axis)
        g = s.vel/c
        
        mnvel = s.__div__(N)

        return RV( vel=mnvel.values() )
    
    def values(self):
        """
        Return a Numpy array of the velocities contained in an RV object.

        Returns:
            numpy.array: radial velocity array
        """
        
        return np.array(self.vel)
