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
        tuple: (A numpy array collapsed along the 0th axis, and weight matrix)
        
    """
    # Calcualte and subtract reference
    ref = cenfunc(s, axis=0)

    if inweights is None:
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

    return m,weights

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
    params = radvel.RVParameters(1,basis='per tc secosw sesinw k')
    params['per1'] = per
    params['tc1'] = tc
    params['secosw1'] = 0.00 
    params['sesinw1'] = 0.00
    params['k1'] = 55.0
    params['dvdt'] = 0
    params['curv'] = 0
    mod = radvel.RVModel(params, time_base=time_base)

    like = radvel.likelihood.RVLikelihood(mod, vdf['jd'], vdf['mnvel'], vdf['errvel'])
    like.params['gamma'] = 0.0
    like.params['jit'] = 1.0

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

def timebin(time, meas, meas_err, binsize):
#  This routine bins a set of times, measurements, and measurement errors 
#  into time bins.  All inputs and outputs should be floats or double. 
#  binsize should have the same units as the time array.
#  - from Andrew Howard, ported to Python by BJ Fulton

    ind_order = np.argsort(time)
    time = time[ind_order]
    meas = meas[ind_order]
    meas_err = meas_err[ind_order]
    ct=0
    while ct < len(time):
        ind = np.where((time >= time[ct]) & (time < time[ct]+binsize))[0]
        num = len(ind)
        wt = (1./meas_err[ind])**2.     #weights based in errors
        wt = wt/np.sum(wt)              #normalized weights
        if ct == 0:
            time_out = [np.sum(wt*time[ind])]
	    meas_out = [np.sum(wt*meas[ind])]
	    meas_err_out = [1./np.sqrt(np.sum(1./(meas_err[ind])**2))]
        else:
            time_out.append(np.sum(wt*time[ind]))
	    meas_out.append(np.sum(wt*meas[ind]))
	    meas_err_out.append(1./np.sqrt(np.sum(1./(meas_err[ind])**2)))
        ct += num

    return time_out, meas_out, meas_err_out


def wp(data, wt, percentiles):
    """Compute weighted percentiles. 

    If the weights are equal, this is the same as normal percentiles. 
    Elements of the C{data} and C{wt} arrays correspond to 
    each other and must have equal length (unless C{wt} is C{None}). 
   
    Args:
        data (array-like): The data. 
        wt (None or array-like): How important is a given piece of data. 
            All the weights must be non-negative and the sum must be 
            greater than zero.
        percentiles: what percentiles to use.  (Not really percentiles, 
            as the range is 0-1 rather than 0-100.) 

    Returns:
        float: the weighted percentiles of the data. 
    """
    import numpy
    assert numpy.greater_equal(percentiles, 0.0).all(), "Percentiles less than zero" 
    assert numpy.less_equal(percentiles, 1.0).all(), "Percentiles greater than one" 
    data = numpy.asarray(data) 
    assert len(data.shape) == 1 
    if wt is None: 
          wt = numpy.ones(data.shape, numpy.float) 
    else: 
          wt = numpy.asarray(wt, numpy.float) 
          assert wt.shape == data.shape 
          assert numpy.greater_equal(wt, 0.0).all(), "Not all weights are non-negative." 
    assert len(wt.shape) == 1 
    n = data.shape[0] 
    assert n > 0 
    i = numpy.argsort(data) 
    sd = numpy.take(data, i, axis=0) 
    sw = numpy.take(wt, i, axis=0) 
    aw = numpy.add.accumulate(sw) 
    if not aw[-1] > 0: 
          raise ValueError, "Nonpositive weight sum" 
    w = (aw-0.5*sw)/aw[-1] 
    spots = numpy.searchsorted(w, percentiles) 
    o = [] 
    for (s, p) in zip(spots, percentiles): 
          if s == 0: 
                  o.append(sd[0]) 
          elif s == n: 
                  o.append(sd[n-1]) 
          else: 
                  f1 = (w[s] - p)/(w[s] - w[s-1]) 
                  f2 = (p - w[s-1])/(w[s] - w[s-1]) 
                  assert f1>=0 and f2>=0 and f1<=1 and f2<=1 
                  assert abs(f1+f2-1.0) < 1e-6 
                  o.append(sd[s-1]*f1 + sd[s]*f2) 
    return o


def running_median(x,y, width, percentile=0.5):
    """
    Return a running-median of a time-series.

    Args:
        x (array): independant variable
        y (array): dependant variable
        width (float): width of median window in units of x

    Returns:
        array: median filtered version of the data
    """

    filtered = np.zeros_like(x)
    for i in range(len(x)):
        spot = x[i]
        diff = x - spot
        w = np.exp(-((x-spot)**2)/(2*width**2))
        win = np.where(np.abs(diff) <= 5*width)[0]
        #med = np.median(y[win])
        med = wp(y,w, [percentile])[0]
        #print y[i], med
        filtered[i] = med

    return filtered
