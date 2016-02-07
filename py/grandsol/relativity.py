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

