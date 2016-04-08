import pylab as pl
from matplotlib import cm
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import os
import subprocess
import grandsol

default_size = (14,10)
cmap = cm.jet

def velplot_mean(vdf, fmt='s', color='k'):
    """
    The standard RV time series plot for `iGrand`

    Args:
        vdf (DataFrame): containing 'jd', 'mnvel', and 'errvel' columns at a minumum
        fmt (string): (optional) matplotlib.plot format code to determine marker shape
        color (string): (optional) marker color

    Returns:
        None
    """
    
    pl.errorbar(vdf['jd'], vdf['mnvel'], yerr=vdf['errvel'], fmt=fmt, color=color, markersize=10, markeredgewidth=1)
    pl.ylabel('RV [m$^{-1}$]')
    pl.xlabel('HJD$_{\\rm UTC}$ - 2440000')

def velplot_by_order(runname, obdf, orders, outfile=None, vsbc=False):
    """

    Plot the RV time series for all orders from a single ``iGrand`` iteration.

    Args:
        runname (string): name of the iGrand run (e.g. iGrand_sun or iGrand_4628)
        obdf (DataFrame): observation list data frame as output by ``grandsol.io.read_obslist``
        orders (list): list of orders to combine to derive the RVs for each iteration
        outfile (string): (optional) name of the output file. If not given the plot
        will be displayed in an interactive window.
        vsbc (bool): (optional) plot RVs with barycentric correction on the x-axis instead of time

    Returns:
        None

    """

    
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

def truthplot(runname, truthvel, orders, iters=[1,2,3,4,5,6,7,8,9,10], outfile=None):
    """

    Plot the RV residuals relative to the input `truth` velocity as provided in the input obslist.

    Args:
        runname (string): name of the iGrand run (e.g. iGrand_sun or iGrand_4628)
        truthvel (array): known radial velocity for each observation
        orders (list): list of orders to combine to derive the RVs for each iteration
        iters (list): (optional) list of iteration numbers (starting with 1, not 0)
        outfile (string): (optional) name of the output file. If not given the plot
        will be displayed in an interactive window.

    Returns:
        None

    """

    
    fig = pl.figure(figsize=default_size)
    colors = [ cmap(x) for x in np.linspace(0.05, 0.95, len(iters))]
    
    workdir = os.getcwd()
    sigmas = []
    for i in iters:
        idir = "iter%02d" % i
        if os.path.exists(idir):
            os.chdir(idir)
            obdf = grandsol.io.read_obslist('obslist_%02d' % i)
        else:
            print "WARNING: %s does not exist" % idir
            continue
            
        try:
            vdf = grandsol.io.combine_orders(runname, obdf, orders)
            vdf['truthvel'] = truthvel
            vdf['diff'] = vdf['mnvel'] - vdf['truthvel']
        except IOError:
            continue

        sdf = vdf.sort_values('truthvel')

        pl.plot(sdf['truthvel'], sdf['diff'], '-', color=colors[i-1], lw=2)
        sigmas.append(grandsol.utils.MAD(vdf['diff']))
        
        os.chdir(workdir)

    legendlabels = ["iter%d: $\sigma_m=%.2f$ m s$^{-1}$" % (i, s) for i,s in zip(iters,sigmas)]
    
    pl.legend(legendlabels, loc='upper right', fontsize=14)
    pl.title(runname + "  residuals relative to input")
    pl.ylabel('measured velocity - input velocity [m/s]')
    pl.xlabel('input velocity [m/s]')

    if outfile == None: pl.show()
    else: pl.savefig(outfile)


def velplot_by_iter(runname, orders, iters=[1,2,3,4,5,6,7,8,9,10], outfile=None):
    """

    Plot the RV time-series for ``iGrand`` iterations.

    Args:
        runname (string): name of the iGrand run (e.g. iGrand_sun or iGrand_4628)
        orders (list): list of orders to combine to derive the RVs for each iteration
        iters (list): (optional) list of iteration numbers (starting with 1, not 0)
        outfile (string): (optional) name of the output file. If not given the plot
        will be displayed in an interactive window.

    Returns:
        None

    """

    
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
            #print i, diff
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
    """

    Fit an RV orbit model and plot the ``grand`` RVs as a function of orbital phase for each ``iGrand`` iteration.

    Args:
        runname (string): name of the iGrand run (e.g. iGrand_sun or iGrand_4628)
        obdf (DataFrame): observation list data frame as output by ``grandsol.io.read_obslist``
        orders (list): list of orders to combine to derive the RVs for each iteration
        tc (float): time of inferior conjunction (i.e. time of transit if transiting)
        per (float): orbital period in days
        iters (list): (optional) list of iteration numbers (starting with 1, not 0)
        outfile (string): (optional) name of the output file. If not given the plot
        will be displayed in an interactive window.

    Returns:
        None

    """

    
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
            #print i, diff
        except IOError:
            continue

        vdf['jd'] += 2440000
        post = grandsol.utils.orbfit(vdf, tc, per)
        tc = post.params['tc1']
        per = post.params['per1']
        
        phase = grandsol.utils.foldData(vdf['jd'], tc, per, cat=True) - 1
        vcat = np.append(vdf['mnvel'], vdf['mnvel'])
        ecat = np.append(vdf['errvel'], vdf['errvel'])

        modt = np.linspace(-0.5, 0.5, 10000) * per + tc
        mod = post.likelihood.model(modt)
        modp = grandsol.utils.foldData(modt, tc, per, cat=True) - 1
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
    
    pl.legend(legendlabels, loc='best', fontsize=12)
    pl.title(runname + " iterations")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)

def plot_residuals_byobs(modfile, outfile=None):
    """

    Plot spectral flux residuals for all observations on a single plot for each order.

    Args:
        modfile (string): name of the .mod file that contains the residuals for all observations and a single order (e.g. iGrand_sun.08.99.mod)
        outfile (string): (optional) name of the output file. If not given the plot
        will be displayed in an interactive window.

    Returns:
        None

    """

    
    print "Plotting residuals contained in %s" % modfile
    
    model = grandsol.io.read_modfile(modfile)
    order = int(modfile.split('.')[1])
    
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
        pl.title('order %d' % order)
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
    """

    Plot the MAD of the flux residuals as a function of ``iGrand`` iteration number.

    Args:
        runname (string): name of the iGrand run (e.g. iGrand_sun or iGrand_4628)
        orders (list): list of orders to include on the plot
        iters (list): (optional) list of iteration numbers (starting with 1, not 0)
        outfile (string): (optional) name of the output file. If not given the plot
        will be displayed in an interactive window.

    Returns:
        None

    """

    
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


def plot_template_byiter(runname, orders, iters=[1,2,3,4,5,6,7,8,9,10]):
    """

    Plot template evolution as a funciton of ``iGrand`` iteration.
    This will create a multi-page PDF showing the changes in the template
    relative to the last iteration.

    Args:
        runname (string): Name of the iGrand run (e.g. iGrand_sun or iGrand_4628)
        orders (list): list of orders to make plots. Will create a separate multi-page PDF file for each order
        iters (list): list of iteration numbers (starting with 1, not 0)

    Returns:
        None

    """

    
    pl.clf()
    
    zoomwidth = 5.0
    
    def _specplot(temp):
        pl.plot(temp['wav'], temp['solar'], '-', color='b', alpha=0.3, lw=2)
        pl.plot(temp['wav'], temp['temp'], 'k-', lw=2)
        pl.plot(temp['wav'], temp['temp_prev'], '--', color='0.7', lw=2)
    
        pl.plot(temp['wav'], temp['diff'], '-', color='red', lw=2)

        pl.ylim(-0.2, 1.05)
    
        ax = pl.gca()
    
        return ax

    for o in orders:
        with PdfPages('%s_%02d_temp_byiter.pdf' % (runname, o)) as pdf:
            for i in iters:

                tempfile = 'iter%02d/%s.%02d.99.tem' % (i, runname, o)
                prev_tempfile = 'iter%02d/%s.%02d.99.tem' % (i-1, runname, o)
                temp = grandsol.io.read_temfile(tempfile)
        
                if i == 1: temp.temp_prev = np.ones_like(temp.temp)
                else: temp.temp_prev = grandsol.io.read_temfile(prev_tempfile).temp
        
                deep_line = temp.wav[temp.temp.argmin()]

                zoomreg = np.array([-1,1]) * zoomwidth + deep_line

                temp_lastiter = temp.temp_prev
        
                temp['diff'] = (temp.temp - temp.temp_prev) * 100
                temp['diff'] -= temp['diff'].mean()

                pl.subplot(2,1,1)
                ax = _specplot(temp)
                pl.axvspan(*zoomreg, color='0.8')
                pl.xlim(temp.wav[1:-1].min(), temp.wav[1:-1].max())
                pl.title(tempfile)

                pl.subplot(2,1,2)
                ax = _specplot(temp)
                pl.xlim(zoomreg)
                pl.xlabel('Wavelength [$\AA$]')

                pl.legend(['solar spectrum', 'current template', 'previous template', '100$\\times$(current-previous)'], loc='best')

                pdf.savefig()
                pl.close()
                
                
def plot_lsf_byiter(runname, iobs, order, iters=[1,2,3,4,5,6,7,8,9,10]):
    """

    Plot LSF evolution for a single echelle order and single observation as a funciton of ``iGrand`` iteration.
    Colors move from blue to red with increasing iteration number.

    Args:
        runname (string): Name of the iGrand run (e.g. iGrand_sun or iGrand_4628)
        iobs (int): observation index from obslist (starting with 1, not 0)
        order (int): order index (starting with 1, not 0)
        iters (list): list of iteration numbers (starting with 1, not 0)

    Returns:
        matplotlib.figure.Figure: resulting matplotlib figure object

    """

    colors = [ cmap(x) for x in np.linspace(0.05, 0.95, len(iters))]
    
    grlsf_binary = os.path.join(os.environ['GRAND'],"bin","grlsf")

    fig = pl.figure(figsize=(20,10))
    pl.subplots_adjust(wspace=0, hspace=0, right=0.95)

    for i in iters:
        lsffile = os.path.join("iter%02d" % i, "%s.%02d.99.lsf" % (runname, order))

        cmd = [grlsf_binary, lsffile, str(iobs), str(order), '0']

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.wait()

        lsfdf = grandsol.io.read_grlsf(p.stdout)

        numnodes = lsfdf['node'].max()

        nodegroups = lsfdf.groupby('node')
        if i == iters[0]:
            prevgroups = nodegroups
            continue

        pltindex = 1

        axlist = []
        for n in nodegroups.groups.keys():
            nodelsf = nodegroups.get_group(n)

            pl.subplot(2, numnodes, pltindex)
            pl.plot(nodelsf['dj'], nodelsf['lsf'], '-', lw=2, color=colors[i-1])

            ax = pl.gca()
            axlist.append(ax)
            if n == 1:
                pl.ylabel('PSF')

            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticklabels([])

            pltindex += 1

        for n in nodegroups.groups.keys():
            nodelsf = nodegroups.get_group(n)
            prevlsf = prevgroups.get_group(n)

            pl.subplot(2, numnodes, pltindex, sharex=axlist[n-1])
            pl.plot(nodelsf['dj'], nodelsf['lsf']-prevlsf['lsf'], lw=2, color=colors[i-1])

            ax = pl.gca()
            if n == 1:
                pl.ylabel('PSF$_{i}$ - PSF$_{i-1}$')

            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticklabels([])

            pl.xlabel('node %d' % n)

            pltindex += 1

        prevgroups = nodegroups

    pl.annotate('Pixel Offset', xy=(0.5, 0.03), xycoords='figure fraction', horizontalalignment='center', fontsize=24)
    pl.suptitle('order: %d, observation index: %d' % (order, iobs))

    return pl.gcf()
