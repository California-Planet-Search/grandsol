import pylab as pl
import numpy as np
import os
import grandsol

default_size = (12,8)

def fit51peg(vdf):
    import radvel
    import copy
    from scipy import optimize
    
    time_base = 2450000
    params = radvel.RVParameters(1,basis='per tc secosw sesinw logk')
    params['per1'] = 4.230785
    params['tc1'] = 2450003.9156
    params['secosw1'] = 0.00 
    params['sesinw1'] = 0.00
    params['logk1'] = np.log(55.)
    params['dvdt'] = 0
    params['curv'] = 0
    mod = radvel.RVModel(params, time_base=time_base)

    like = radvel.likelihood.RVLikelihood(mod, vdf['jd'], vdf['mnvel'], vdf['errvel'])
    like.params['gamma'] = 0.0
    like.params['logjit'] = np.log(1)

    #like.vary['logjit'] = False
    #like.vary['secosw1'] = False
    #like.vary['sesinw1'] = False
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
    tt = bjd - e
    ncycles = tt/p
    roundcycles = np.floor(ncycles)
    phase = (tt-(roundcycles*p))/p

    if cat:
        phase = np.concatenate((phase,phase+1))

    return phase

def velplot_mean(vdf, fmt='s'):
    pl.errorbar(vdf['jd'], vdf['mnvel'], yerr=vdf['errvel'], fmt=fmt, markersize=10, markeredgewidth=1)
    pl.ylabel('RV [m$^{-1}$]')
    pl.xlabel('HJD$_{\\rm UTC}$ - 2440000')

def velplot_by_order(runname, obdf, orders, outfile=None):
    fig = pl.figure(figsize=default_size)
    vdf, relvel = grandsol.io.combine_orders(runname, obdf, orders, varr_byorder=True)

    sigmas = []
    for i,o in enumerate(orders):
        pl.plot(vdf['jd'], relvel[i,:], 'o')
        sigmas.append(np.std(relvel[i,:]))

    velplot_mean(vdf)

    legendlabels = ["order %d $\sigma=%.2f$ m s$^{-1}$" % (i, s) for i,s in zip(orders,sigmas)] + ['Mean $\sigma=%.2f$' % np.std(vdf['mnvel'])]
    
    pl.legend(legendlabels, loc='best')
    pl.title(runname + " orders")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)

def velplot_by_iter(runname, obdf, orders, iters=[1,2,3,4,5,6,7,8,9,10], outfile=None):
    fig = pl.figure(figsize=default_size)
    workdir = os.getcwd()
    prev = 0.0
    sigmas = []
    for i in iters:
        idir = "iter%02d" % i
        if os.path.exists(idir): os.chdir(idir)
        else: continue
            
        try:
            vdf = grandsol.io.combine_orders(runname, obdf, orders)
            diff = np.sum(((vdf['mnvel'] - prev) / vdf['errvel'])**2)
            prev = vdf['mnvel']
            print i, diff
        except IOError:
            continue
        velplot_mean(vdf, fmt='s')
        sigmas.append(np.std(vdf['mnvel']))
        
        os.chdir(workdir)

    legendlabels = ["iteration %d\n$\sigma=%.2f$ m s$^{-1}$" % (i, s) for i,s in zip(iters,sigmas)]
    
    pl.legend(legendlabels, loc='best')
    pl.title(runname + " iterations")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)

def phaseplot_by_iter(runname, obdf, orders, tc, per, iters=[1,2,3,4,5,6,7,8,9,10], outfile=None):
    fig = pl.figure(figsize=default_size)
    workdir = os.getcwd()
    prev = 0.0
    sigmas = []
    for i in iters:
        idir = "iter%02d" % i
        if os.path.exists(idir): os.chdir(idir)
        else: continue
            
        try:
            vdf = grandsol.io.combine_orders(runname, obdf, orders)
            diff = np.sum(((vdf['mnvel'] - prev) / vdf['errvel'])**2)
            prev = vdf['mnvel']
            print i, diff
        except IOError:
            continue

        vdf['jd'] += 2440000
        phase = foldData(vdf['jd'], tc, per, cat=True) - 1
        vcat = np.append(vdf['mnvel'], vdf['mnvel'])
        ecat = np.append(vdf['errvel'], vdf['errvel'])

        post = fit51peg(vdf)
        modt = np.linspace(-0.5, 0.5, 10000) * per + tc
        mod = post.likelihood.model(modt)
        modp = foldData(modt, tc, per, cat=True) - 1
        omod = np.argsort(modp)
        mod = np.append(mod,mod)[omod]

        ebar = pl.errorbar(phase, vcat, yerr=ecat, fmt='s')
        pl.plot(modp[omod], mod, color=ebar[0].get_color(), lw=2, label='_nolegend_')
        pl.xlim(-0.5, 0.5)
        pl.xlabel('Phase')
        pl.ylabel('RV m s$^{-1}$')

        sigmas.append(np.exp(post.params['logk1']))
        
        os.chdir(workdir)

    legendlabels = ["iteration %d\n$K=%.2f$ m s$^{-1}$" % (i, s) for i,s in zip(iters,sigmas)]
    
    pl.legend(legendlabels, loc='best')
    pl.title(runname + " iterations")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)
