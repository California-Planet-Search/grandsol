import numpy.ma as ma
import numpy as np

def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation (MAD) along given axis of an array:

    .. math:: \\frac{median(|a - median(a)|)}{c}

    Args:
        a (array): input array
        c (float): The constant to convert from MAD to std, 0.6745 is used by default
        axis (int): axis along which to calculate the MAD

    Returns:
        float: MAD
        
    """
    

    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m

def clipped_mean(s, sigma=5, inweights=None, iterations=10, cenfunc=np.median, varfunc=MAD, verbose=False):
    """
    Calculate a weighted sigma-clipped mean along the 0th axis of a 2D array.

    Args:
        s (float array): input array
        sigma (float): Number of standard deviations to use as the clipping limit
        inweights (float array): Weights to use in the mean calculation. Defaults to equal weights
        iterations (int): Number of sigma-clipping iterations to run
        cenfunc (callable): The technique to compute the center for the clipping.
                            Must be a callable that takes in a 1D data array and outputs the central value.
                            Defaults to numpy.median
        varfunc (callable): The technique to compute the variance about the center. Must be a callable that
                            takes in a 1D data array and outputs the width estimator that will be interpreted
                            as a variance. Defaults to the MAD.
        verbose (bool): Print message indicating number of masked outliers
    Returns:
        numpy.masked.MaskedArray: A masked array with a shape matching the
        input that is masked where the algorithm has rejected those values.
        
    """
    # Calcualte and subtract reference
    ref = cenfunc(s, axis=0)

    if inweights == None:
        inweights = np.ones(s.shape[0], dtype=float)
        
    for n in range(iterations):
        weights = np.ones_like(s, dtype=float) * inweights[:,np.newaxis]
        for i in range(s.shape[0]):
                                    
            diff = s[i,:] - ref

            bad = (np.abs(diff) >= sigma*varfunc(diff))
            
            w = np.ones_like(diff, dtype=float)
            w[bad] = 0
            weights[i,:] *= w
            if verbose: print "Masking %d bad measurements." % np.where(w==0)[0].shape[0]
        
            m = np.average(s, axis=0, weights=weights)
        
        ref = m.copy()

    return m
