import pylab as pl
from matplotlib import cm
import matplotlib
import numpy as np
import pandas as pd
import os
import grandsol

default_size = (14,10)
cmap = cm.jet

def fit51peg(vdf):
    import radvel
    import copy
    from scipy import optimize
    
    time_base = 2450000
    params = radvel.RVParameters(1,basis='per tc secosw sesinw logk')
    params['per1'] = 4.230785
    params['tc1'] = 2450001.9156
    params['secosw1'] = 0.00 
    params['sesinw1'] = 0.00
    params['logk1'] = np.log(55.)
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
    tt = bjd - e
    ncycles = tt/p
    roundcycles = np.floor(ncycles)
    phase = (tt-(roundcycles*p))/p

    if cat:
        phase = np.concatenate((phase,phase+1))

    return phase

def velplot_mean(vdf, fmt='s', color='k'):
    pl.errorbar(vdf['jd'], vdf['mnvel'], yerr=vdf['errvel'], fmt=fmt, color=color, markersize=10, markeredgewidth=1)
    pl.ylabel('RV [m$^{-1}$]')
    pl.xlabel('HJD$_{\\rm UTC}$ - 2440000')

def velplot_by_order(runname, obdf, orders, outfile=None, vsbc=False):
    fig = pl.figure(figsize=default_size)
    colors = [ cmap(x) for x in np.linspace(0.05, 0.95, len(orders))]
    
    vdf, relvel = grandsol.io.combine_orders(runname, obdf, orders, varr_byorder=True)

    sigmas = []
    for i,o in enumerate(orders):
        if vsbc:
            pl.plot(obdf['bc'], relvel[i,:], 'o', color=colors[i])
        else:
            pl.plot(vdf['jd'], relvel[i,:], 'o', color=colors[i])
        #sigmas.append(np.std(relvel[i,:]))
        sigmas.append(grandsol.utils.MAD(relvel[i,:]))

    if vsbc:
        pl.errorbar(obdf['bc'], vdf['mnvel'], yerr=vdf['errvel'], fmt='s', color=colors[i], markersize=10, markeredgewidth=1)
        pl.ylabel('RV [m$^{-1}$]')
        pl.xlabel('BC [m$^{-1}$]')
    else:
        velplot_mean(vdf, color=colors[i])

    legendlabels = ["order %d $\sigma_m=%.2f$ m s$^{-1}$" % (i, s) for i,s in zip(orders,sigmas)] + ['Mean $\sigma=%.2f$' % grandsol.utils.MAD(vdf['mnvel'])]
    
    pl.legend(legendlabels, loc='best')
    pl.title(runname + " orders")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)

def velplot_by_iter(runname, orders, iters=[1,2,3,4,5,6,7,8,9,10], outfile=None):
    fig = pl.figure(figsize=default_size)
    colors = [ cmap(x) for x in np.linspace(0.05, 0.95, len(iters))]
    
    workdir = os.getcwd()
    prev = 0.0
    sigmas = []
    for i in iters:
        idir = "iter%02d" % i
        if os.path.exists(idir):
            os.chdir(idir)
            obdf = grandsol.io.read_obslist('obslist_%02d' % i)
        else: continue
            
        try:
            vdf = grandsol.io.combine_orders(runname, obdf, orders)
            diff = np.sum(((vdf['mnvel'] - prev) / vdf['errvel'])**2)
            prev = vdf['mnvel']
            print i, diff
        except IOError:
            continue

        velplot_mean(vdf, fmt='s', color=colors[i-1])
        #sigmas.append(np.std(vdf['mnvel']))
        sigmas.append(grandsol.utils.MAD(vdf['mnvel']))
        
        os.chdir(workdir)

    legendlabels = ["iteration %d\n$\sigma_m=%.2f$ m s$^{-1}$" % (i, s) for i,s in zip(iters,sigmas)]
    
    pl.legend(legendlabels, loc='best')
    pl.title(runname + " iterations")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)

def phaseplot_by_iter(runname, obdf, orders, tc, per, iters=[1,2,3,4,5,6,7,8,9,10], outfile=None):
    fig = pl.figure(figsize=default_size)
    colors = [ cmap(x) for x in np.linspace(0.05, 0.95, len(iters))]
    
    workdir = os.getcwd()
    prev = 0.0
    sigmas = []
    Klist = []
    for i in iters:
        idir = "iter%02d" % i
        if os.path.exists(idir):
            os.chdir(idir)
            obdf = grandsol.io.read_obslist('obslist_%02d' % i)
        else: continue
            
        try:
            vdf = grandsol.io.combine_orders(runname, obdf, orders)
            diff = np.sum(((vdf['mnvel'] - prev) / vdf['errvel'])**2)
            prev = vdf['mnvel']
            print i, diff
        except IOError:
            continue

        vdf['jd'] += 2440000
        post = fit51peg(vdf)
        tc = post.params['tc1']
        per = post.params['per1']
        
        phase = foldData(vdf['jd'], tc, per, cat=True) - 1
        vcat = np.append(vdf['mnvel'], vdf['mnvel'])
        ecat = np.append(vdf['errvel'], vdf['errvel'])

        modt = np.linspace(-0.5, 0.5, 10000) * per + tc
        mod = post.likelihood.model(modt)
        modp = foldData(modt, tc, per, cat=True) - 1
        omod = np.argsort(modp)
        mod = np.append(mod,mod)[omod]

        ebar = pl.errorbar(phase, vcat, yerr=ecat, fmt='s', color=colors[i-1])
        pl.plot(modp[omod], mod, color=ebar[0].get_color(), lw=2, label='_nolegend_')
        pl.xlim(-0.5, 0.5)
        pl.xlabel('Phase')
        pl.ylabel('RV m s$^{-1}$')

        Klist.append(np.exp(post.params['logk1']))
        sigmas.append(grandsol.utils.MAD(post.likelihood.residuals()))
        
        os.chdir(workdir)

    legendlabels = ["iteration %d\n$K=%.1f$ MAD=%.1f m s$^{-1}$" % (i,k,s) for i,k,s in zip(iters,Klist,sigmas)]
    
    pl.legend(legendlabels, loc='best')
    pl.title(runname + " iterations")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)


def plot_mean_residuals(modfile):

    outdf = pd.DataFrame()
        
    model = grandsol.io.read_modfile(modfile)
    
    pixmean = model.groupby('pixel', as_index=False).mean()
    pixmean['residuals'] = pixmean['residuals'] = (pixmean['spec'] - (pixmean['model']*pixmean['cont'])) / pixmean['smooth_cont']
    pixmean['residuals_x10'] = pixmean['residuals'] * 10
    pixmean['residuals_percent'] = pixmean['residuals'] * 100
    pixmean['normspec'] = pixmean['spec'] / pixmean['smooth_cont']
    pixmean['normmod'] = pixmean.residuals + pixmean.normspec
    
    #pixmean.plot('wav_star', 'normspec', color='k', lw=0.5)
    #ax = pl.gca()
    #pixmean.plot('wav_star', 'normmod', color='b', linestyle='--', ax=ax)
    pixmean.plot('wav_star', 'residuals_percent', color='r', lw=1)
    #ax = pl.gca()
    #pixmean.plot('wav_obs', 'residuals_percent', color='b', lw=1, ax=ax)

    pl.annotate('$\sigma_m$ residuals = %.3g %%' % grandsol.utils.MAD(pixmean.residuals_percent.values), xy=(0.7, 0.05), xycoords='axes fraction')
    pl.ylabel('Relative Flux [%]')
    #pl.ylim(-0.003, 0.003)
    pl.ylim(-0.3, 0.3)
    pl.title(modfile)

def plot_residuals_byobs(modfile, outfile=None):
    print "Plotting residuals contained in %s" % modfile
    
    model = grandsol.io.read_modfile(modfile)

    model['residuals'] = (model['spec'] - (model['model']*model['cont'])) / model['smooth_cont']
    model['residuals_percent'] = model['residuals'] * 100

    byobs = model.groupby('ind', as_index=False)

    fig = pl.figure(figsize=(16,12))
    pl.subplot(211)
    pl.subplots_adjust(hspace=0.25)

    for group in byobs.groups:
        singleobs = pd.DataFrame(byobs.get_group(group))

        pl.subplot(211)
        pl.plot(singleobs['wav_obs'],singleobs['residuals_percent'], 'k.', markersize=0.4, rasterized=True)
        ax_obs = pl.gca()
    
        pl.subplot(212)
        pl.plot(singleobs['wav_star'],singleobs['residuals_percent'], 'k.', markersize=0.4, rasterized=True)
        ax_star = pl.gca()

    
    rms = model.residuals_percent.std()
    mad = grandsol.utils.MAD(model.residuals_percent.values)
    
    waverange = model['wav_obs'].max() - model['wav_obs'].min()
    crop_low = model['wav_obs'].min() + 0.05*waverange
    crop_high = model['wav_obs'].max() - 0.05*waverange
    rms_crop = model[(model['wav_obs'] > crop_low) & (model['wav_obs'] < crop_high)].residuals_percent.std()

    ax_obs.axvline(crop_low, color='r', linestyle='dashed')
    ax_obs.axvline(crop_high, color='r', linestyle='dashed')
    ax_star.axvline(crop_low, color='r', linestyle='dashed')
    ax_star.axvline(crop_high, color='r', linestyle='dashed')


    ax_star.annotate("$\sigma$ = %.3f, $\sigma_{\\rm crop}$ = %.3f , MAD = %.3f %%" % (rms, rms_crop, mad) , xy=(0.2, 0.05), xycoords='axes fraction')
    
    ax_obs.set_ylim(-5*mad, 5*mad)
    ax_obs.set_xlim(singleobs['wav_obs'].min(), singleobs['wav_obs'].max())
    ax_obs.set_xlabel('Wavelength in observatory frame [$\AA$]')
    ax_obs.set_ylabel('Relative flux residuals [%]')

    ax_star.set_ylim(-5*mad, 5*mad)
    ax_star.set_xlim(singleobs['wav_star'].min(), singleobs['wav_star'].max())
    ax_star.set_xlabel('Wavelength in stellar frame [$\AA$]')
    ax_star.set_ylabel('Relative flux residuals [%]')

    if outfile == None: pl.show()
    else: pl.savefig(outfile)


def plot_resMAD_byiter(runname, obdf, orders, iters=[1,2,3,4,5,6,7,8,9,10], outfile=None):

    colors = [ cmap(x) for x in np.linspace(0.05, 0.95, len(orders))]
    MADarr = np.zeros((len(iters), len(orders)))
        
    for i,iteration in enumerate(iters):
        os.chdir('iter%02d' % iteration)
        for j,o in enumerate(orders):
            modfile = "%s.%02d.99.mod" % (runname, o)
            model = grandsol.io.read_modfile(modfile)
            model['residuals'] = (model['spec'] - (model['model']*model['cont'])) / model['smooth_cont']
            model['residuals_percent'] = model['residuals'] * 100

            mad = grandsol.utils.MAD(model['residuals_percent'])
            
            MADarr[i, j] = mad
            print iteration, o, mad
        os.chdir('..')
        
    #MADarr /= np.mean(MADarr, axis=0)
    #MADarr -= 1
    
    fig = pl.figure(figsize=default_size)

    for i,c in enumerate(colors):
        pl.plot(iters, MADarr[:,i], 'o-', markersize=10, color=c)

    pl.legend(['order %d' % o for o in orders])
    pl.xlabel('iteration')
    pl.ylabel('MAD of residuals [%]')
            
    if outfile == None: pl.show()
    else: pl.savefig(outfile)
