import pylab as pl
import numpy as np
import os
import grandsol

def velplot_mean(vdf, fmt='ks'):
    pl.errorbar(vdf['jd'], vdf['mnvel'], yerr=vdf['errvel'], fmt=fmt, markersize=10)
    pl.ylabel('RV [m$^{-1}$]')
    pl.xlabel('HJD$_{\\rm UTC}$ - 2440000')

def velplot_by_order(runname, obdf, orders, outfile=None):
    fig = pl.figure()
    vdf, relvel = grandsol.io.combine_orders(runname, obdf, orders, varr_byorder=True)

    sigmas = []
    for i,o in enumerate(orders):
        pl.plot(vdf['jd'], relvel[i,:], 'o')
        sigmas.append(np.std(relvel[i,:]))

    velplot_mean(vdf)
    pl.legend(orders+['Mean'], loc='best')
    pl.title(runname + " orders")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)

def velplot_by_iter(runname, obdf, orders, iters=[1,2,3], outfile=None):
    fig = pl.figure()
    workdir = os.getcwd()
    sigmas = []
    for i in iters:
        os.chdir("iter%02d" % i)
        
        vdf = grandsol.io.combine_orders(runname, obdf, orders)
        velplot_mean(vdf, fmt='s')
        sigmas.append(np.std(vdf['mnvel']))
        
        os.chdir(workdir)

    legendlabels = ["iteration %d\n$\sigma=%.2f$ m s$^{-1}$" % (i, s) for i,s in zip(iters,sigmas)]
    
    pl.legend(legendlabels, loc='best')
    pl.title(runname + " iterations")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)

        
