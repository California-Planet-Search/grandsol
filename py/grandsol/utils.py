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

def orbfit(vdf, tc, per):
    """
    Wrapper function to perform a maximum-liklihood fit of an orbital model to RVs contained in a DataFrame output by
    ``grandsol.io.combine_orders``. Depends on the ``radvel`` radial velocity fitting package be installed in in the PYTHONPATH.

    Args:
        vdf (DataFrame): as output by ``grandsol.io.combine_orders``
        tc (float): time of inferior conjunction (i.e. time of transit if transiting)
        per (float): orbital period in days

    Returns:
        radvel.Posterior: Posterior object containing the best-fitting parameters.

    """
    
    import radvel
    import copy
    from scipy import optimize
    
    time_base = 2450000
    params = radvel.RVParameters(1,basis='per tc secosw sesinw logk')
    params['per1'] = per
    params['tc1'] = tc
    params['secosw1'] = 0.00 
    params['sesinw1'] = 0.00
    params['logk1'] = np.log(25.)
    params['dvdt'] = 0
    params['curv'] = 0
    mod = radvel.RVModel(params, time_base=time_base)

    like = radvel.likelihood.RVLikelihood(mod, vdf['jd'], vdf['mnvel'], vdf['errvel'])
    like.params['gamma'] = 0.0
    like.params['logjit'] = np.log(1)

    like.vary['dvdt'] = False
    like.vary['curv'] = False

    post = radvel.posterior.Posterior(like)
    post0 = copy.deepcopy(post)

    post.priors += [radvel.prior.EccentricityPrior( 1 )] # Keeps eccentricity < 1

    res  = optimize.minimize(post.neglogprob_array, post.get_vary_params(), method='Powell',
                         options=dict(maxiter=100000,maxfev=100000,xtol=1e-8) )

    print "Initial loglikelihood = %f" % post0.logprob()
    print "Final loglikelihood = %f" % post.logprob()
    #post.params = post.params.basis.to_cps(post.params)
    print post

    #print post.params.basis.to_cps(post.params)

    return post
    
def foldData(bjd,e,p,cat=False):
    """
    Phase fold an RV time series.

    Args:
        bjd (array): JD timestamps
        e (float): reference epoch for phase folding. Phase is 0 at this time. We commonly use the time of transit.
        p (float): orbital period in days
        cat (bool): (optional) contatenate phase and phase+1 for situations whre you might want to plot two orbits

    Returns:
        array: orbital phase corresponding to the JD timestamps contained in the bjd argument
    """

    
    tt = bjd - e
    ncycles = tt/p
    roundcycles = np.floor(ncycles)
    phase = (tt-(roundcycles*p))/p

    if cat:
        phase = np.concatenate((phase,phase+1))

    return phase
