import pylab as pl
from matplotlib import cm
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.offsetbox import AnchoredOffsetbox
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import AuxTransformBox, VPacker, HPacker, TextArea, DrawingArea
import numpy as np
import pandas as pd
import os
import subprocess
import grandsol

default_size = (15,10)
cmap = cm.jet


class AnchoredScaleBar(AnchoredOffsetbox):
    """Class to support scalebars

    Draw a horizontal and/or vertical  bar with the size in data coordinate
    of the give axes. A label will be drawn underneath (center-aligned).

    Args:
        transform: the coordinate frame (typically axes.transData)
        sizex: width of x bar, in data units. 0 to omit
        sizey: width of y bar, in data units. 0 to omit
        labelx: label for x bar; None to omit
        labely: label for y bar; None to omit
        loc (int): position in containing axes
        pad (float):  padding, in fraction of the legend font size (or prop)
        borderpad (float): padding, in fraction of the legend font size (or prop)
        sep (float): separation between labels and bars in points.
        **kwargs: additional arguments passed to base class constructor (matplotlib.offsetbox.AnchoredOffsetBox)

    Returns:
        None
        
    """

    def __init__(self, transform, sizex=0, sizey=0, labelx=None, labely=None, loc=4,
                 pad=0.1, borderpad=0.1, sep=2, prop=None, **kwargs):
        bars = AuxTransformBox(transform)
        if sizex:
            bars.add_artist(Rectangle((0,0), sizex, 0, fc="none"))
        if sizey:
            bars.add_artist(Rectangle((0,0), 0, sizey, fc="none"))

        if sizex and labelx:
            bars = VPacker(children=[bars, TextArea(labelx, minimumdescent=False)],
                           align="center", pad=0, sep=sep)
        if sizey and labely:
            bars = HPacker(children=[TextArea(labely), bars],
                            align="center", pad=0, sep=sep)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=bars, prop=prop, frameon=False, **kwargs)

def add_scalebar(ax, matchx=True, matchy=True, hidex=True, hidey=True, **kwargs):
    """ Add scalebars to axes

    Adds a set of scale bars to *ax*, matching the size to the ticks of the plot
    and optionally hiding the x and y axes

    Args:
        ax (matplotlib.axes): the axis to attach ticks to
        matchx: if True, set size of scale bars to spacing between ticks, if False, size should be set using sizex
        matchy: if True, set size of scale bars to spacing between ticks, if False, size should be set using sizey
        hidex: if True, hide x-axis of parent
        hidey: if True, hide y-axis of parent
        **kwargs: additional arguments passed to AnchoredScaleBars

    Returns:
        AnchoredScaleBar: created scalebar object

    """
    def f(axis):
        l = axis.get_majorticklocs()
        return len(l)>1 and (l[1] - l[0])
    
    if matchx:
        kwargs['sizex'] = f(ax.xaxis)
        kwargs['labelx'] = str(kwargs['sizex'])
    if matchy:
        kwargs['sizey'] = f(ax.yaxis)
        kwargs['labely'] = str(kwargs['sizey'])
        
    sb = AnchoredScaleBar(ax.transData, **kwargs)
    ax.add_artist(sb)

    if hidex : ax.xaxis.set_visible(False)
    if hidey : ax.yaxis.set_visible(False)

    return sb



def velplot_mean(vdf, fmt='s', color='k', vsbc=False):
    """
    The standard RV time series plot for `iGrand`

    Args:
        vdf (DataFrame): containing 'jd', 'mnvel', and 'errvel' columns at a minumum
        fmt (string): (optional) matplotlib.plot format code to determine marker shape
        color (string): (optional) marker color
        vsbc (bool): (optional) BC on x-axis instead of JD

    Returns:
        None
    """
    if vsbc:
        pl.errorbar(vdf['bc_y'], vdf['mnvel_corr'], yerr=vdf['errvel'], fmt=fmt, color=color, markersize=10, markeredgewidth=1)
        pl.xlabel('BC [m s$^{-1}$]')
    else:
        pl.errorbar(vdf['jd'], vdf['mnvel_corr'], yerr=vdf['errvel'], fmt=fmt, color=color, markersize=10, markeredgewidth=1)
        pl.xlabel('HJD$_{\\rm UTC}$ - 2440000')

    pl.ylabel('RV [m s$^{-1}$]')

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
    weights = grandsol.io.combine_orders(runname, obdf, orders, get_weights=True)
    allbad = weights == 0
        
    sigmas = []
    plist = []
    for i,o in enumerate(orders):
        bad = allbad[i,:]
        if vsbc:
            p, = pl.plot(obdf['bc'], relvel[i,:], 'o', color=colors[i])
            pl.plot(obdf['bc'][bad], relvel[i,:][bad], 'x', markersize=24, color=colors[i])
        else:
            p, = pl.plot(vdf['jd'], relvel[i,:], 'o', color=colors[i])
            pl.plot(vdf['jd'][bad], relvel[i,:][bad], 'x', markersize=24, color=colors[i])
        #sigmas.append(np.std(relvel[i,:]))
        sigmas.append(grandsol.utils.MAD(relvel[i,:]))
        plist.append(p)

    if vsbc:
        #pl.errorbar(obdf['bc'], vdf['mnvel'], yerr=vdf['errvel'], fmt='s', color=colors[i], markersize=10, markeredgewidth=1)
        pl.ylabel('RV [m s$^{-1}$]')
        pl.xlabel('BC [m s$^{-1}$]')
    else:
        velplot_mean(vdf, color=colors[i])

    legendlabels = ["order %d $\sigma_m=%.2f$ m s$^{-1}$" % (i, s) for i,s in zip(orders,sigmas)] + ['Mean $\sigma=%.2f$' % grandsol.utils.MAD(vdf['mnvel'])]
    
    pl.legend(plist, legendlabels, loc='best', fontsize=12)
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


def velplot_by_iter(runname, orders, iters=[1,2,3,4,5,6,7,8,9,10], outfile=None, binsize=2.0, vsbc=False):
    """

    Plot the RV time-series for ``iGrand`` iterations.

    Args:
        runname (string): name of the iGrand run (e.g. iGrand_sun or iGrand_4628)
        orders (list): list of orders to combine to derive the RVs for each iteration
        iters (list): (optional) list of iteration numbers (starting with 1, not 0)
        outfile (string): (optional) name of the output file. If not given the plot
            will be displayed in an interactive window.

    Returns:
        DataFrame: DataFrame with velocities from all orders combined

    """

    
    fig = pl.figure(figsize=default_size)
    colors = [ cmap(x) for x in np.linspace(0.05, 0.95, len(iters))]
    
    workdir = os.getcwd()
    prev = 0.0
    sigmas = []
    for i in iters:
        idir = "iter%02d" % i
        print idir
        if os.path.exists(idir):
            os.chdir(idir)
            obdf = grandsol.io.read_obslist('obslist_%02d' % i)
        else:
            print "WARNING: Could not read obslist_%02d" % i
            continue
            
        try:
            vdf = grandsol.io.combine_orders(runname, obdf, orders)
            vdf_nobin = vdf.copy()
            diff = np.sum(((vdf['mnvel'] - prev) / vdf['errvel'])**2)
            prev = vdf['mnvel']
            #print i, diff
        except IOError:
            print "WARNING: Could not read velocities for iteration %d" % i
            continue

        if binsize > 0.0:
            bintimes, binvels, binerr = grandsol.utils.timebin(vdf['jd'].values,
                                                               vdf['mnvel_corr'].values,
                                                               vdf['errvel'].values,
                                                               binsize=binsize)
            _, binbc, _ = grandsol.utils.timebin(vdf['jd'].values,
                                                   vdf['bc_y'].values,
                                                   vdf['errvel'].values,
                                                   binsize=binsize)

            vdf = pd.DataFrame([])
            vdf['jd'] = bintimes
            vdf['mnvel_corr'] = binvels
            vdf['errvel'] = binerr
            vdf['bc_y'] = binbc

        
        velplot_mean(vdf, fmt='s', color=colors[i-1], vsbc=vsbc)
        #sigmas.append(np.std(vdf['mnvel']))
        sigmas.append(grandsol.utils.MAD(vdf['mnvel']))
        
        os.chdir(workdir)

    legendlabels = ["iteration %d\n$\sigma_m=%.2f$ m s$^{-1}$" % (i, s) for i,s in zip(iters,sigmas)]
    
    pl.legend(legendlabels, loc='best', fontsize=12)
    pl.title(runname + " iterations")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)

    return vdf_nobin

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

        modt = np.linspace(-0.5, 0.5, 10000) * per + tc
        mod = post.likelihood.model(modt)
        modp = grandsol.utils.foldData(modt, tc, per, cat=True) - 1
        omod = np.argsort(modp)
        mod = np.append(mod,mod)[omod]

        vcat = post.likelihood.residuals() + post.likelihood.model(post.likelihood.x)
        vcat = np.append(vcat, vcat)
        ecat = np.append(vdf['errvel'], vdf['errvel'])

        
        ebar = pl.errorbar(phase, vcat, yerr=ecat, fmt='s', color=colors[i-1])
        pl.plot(modp[omod], mod, color=ebar[0].get_color(), lw=2, label='_nolegend_')
        pl.xlim(-0.5, 0.5)
        pl.xlabel('Phase')
        pl.ylabel('RV m s$^{-1}$')

        Klist.append(post.params['k1'])
        sigmas.append(grandsol.utils.MAD(post.likelihood.residuals()))
        
        os.chdir(workdir)

    legendlabels = ["iteration %d\n$K=%.1f$ MAD=%.1f m s$^{-1}$" % (i,k,s) for i,k,s in zip(iters,Klist,sigmas)]
    
    pl.legend(legendlabels, loc='best', fontsize=12)
    pl.title(runname + " iterations")
    if outfile == None: pl.show()
    else: pl.savefig(outfile)

def plot_residuals_byobs(modfile, outfile=None, tellurics=False, iodine=False, maskfile=None):
    """

    Plot spectral flux residuals for all observations on a single plot for each order.

    Args:
        modfile (string): name of the .mod file that contains the residuals for all observations and a single order (e.g. iGrand_sun.08.99.mod)
        outfile (string): (optional) name of the output file. If not given the plot
        will be displayed in an interactive window.
        maskfile (string): (optional) path to file containing wavelength mask
    Returns:
        None

    """

    
    print "Plotting residuals contained in %s" % modfile
    
    model = grandsol.io.read_modfile(modfile)
    order = int(os.path.basename(modfile).split('.')[1])
    
    model['residuals'] = (model['spec'] - (model['model']*model['cont'])) / model['smooth_cont']
    model['residuals_percent'] = model['residuals'].values * 100
    
    byobs = model.groupby('ind', as_index=False)

    fig = pl.figure(figsize=(16,12))
    pl.subplot(211)
    pl.subplots_adjust(hspace=0.25)

    for group in byobs.groups:
        singleobs = pd.DataFrame(byobs.get_group(group))

        pl.subplot(211)
        pl.plot(singleobs['wav_obs'],singleobs['residuals_percent'], 'k.', markersize=0.6, rasterized=True)
        pl.title('order %d' % order)
        ax_obs = pl.gca()
    
        pl.subplot(212)
        pl.plot(singleobs['wav_star'],singleobs['residuals_percent'], 'k.', markersize=0.6, rasterized=True)
        ax_star = pl.gca()


    rms = model.residuals_percent.std()
    mad = grandsol.utils.MAD(model.residuals_percent.values)

    if tellurics:
        tel = pd.read_csv(os.path.join(os.environ['GRAND_IREFDIR'], 'telluric_4490_9000.csv'))
        pl.subplot(211)
        pl.plot(tel['wavl'], tel['spec']*100 - 100, 'g-', lw=2, alpha=0.5)

    if iodine:
        iod = pd.read_csv(os.path.join(os.environ['GRAND_REFDIR'], 'apfiodine.sam'), sep=' ', skipinitialspace=True, names=['index', 'spec', 'wavl'])
        iod['spec'] = (iod['spec'] / iod['spec'].mean()) - 1
        print iod['spec']
        pl.subplot(211)
        pl.plot(iod['wavl'], iod['spec'], 'b-', lw=2, alpha=0.5)

    if maskfile is not None:
        for line in open(maskfile, 'r').readlines():
            line = line.strip().split()
            if len(line) < 3:
                continue
            if int(line[0]) == order or int(line[0]) == 0:
                p1, p2 = line[1], line[2]
                iw = np.interp([p1,p2], model['pixel'], model['wav_obs'])
                print p1, p2, iw[0], iw[1]
                ax_obs.axvspan(iw[0], iw[1], color='0.5', alpha=0.5)
                #ax_obs.axvspan(line[1], line[2], color='0.5', alpha=0.5)
    
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
        
    MAD_norm = MADarr - np.mean(MADarr, axis=0)
    #MAD_norm -= 1

    meanMADs = np.mean(MADarr, axis=0)
    
    fig = pl.figure(figsize=default_size)
    pl.subplots_adjust(left=0.17, right=0.95)

    for i,c in enumerate(colors):
        pl.plot(iters, MAD_norm[:,i], 'o-', markersize=10, color=c)

    pl.legend(['order %d, mean(MAD) = %4.4f %%' % (o,meanMADs[i]) for i,o in enumerate(orders)], loc='best', fontsize=14)
    pl.xlabel('iteration')
    pl.ylabel('$\\Delta$ MAD of residuals [%]')
            
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
                #tempfile = 'iter10/%s.%02d.%02d.tem' % (runname, o, i)
                #prev_tempfile = 'iter10/%s.%02d.%02d.tem' % (runname, o, i-1)
                if not os.path.exists(tempfile):
                    tempfile = 'iter%02d/%s.%02d.10.tem' % (i, runname, o)
                    prev_tempfile = 'iter%02d/%s.%02d.10.tem' % (i-1, runname, o)


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

    fig = pl.figure(figsize=(24,12))
    pl.subplots_adjust(wspace=0, hspace=0, right=0.98, left=0.04, bottom=0.15)

    for i in iters:
        lsffile = os.path.join("iter%02d" % i, "%s.%02d.99.lsf" % (runname, order))

        cmd = [grlsf_binary, lsffile, str(iobs), str(order), '0']

        print ' '.join(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        lsfdf = grandsol.io.read_grlsf(p.stdout)

        numnodes = lsfdf['node'].max()

        nodegroups = lsfdf.groupby('node')
        if i == iters[0]:
            prevgroups = nodegroups
            continue

        pltindex = 1

        axlist = []
        centroids = []
        for n in nodegroups.groups.keys():
            nodelsf = nodegroups.get_group(n)

            pl.subplot(3, numnodes, pltindex)
            pl.plot(nodelsf['dj'], nodelsf['lsf'], '-', lw=2, color=colors[i-1])

            cen = np.sum(nodelsf['dj']*nodelsf['lsf']) / np.sum(nodelsf['lsf'])
            centroids.append(cen)
            
            ax = pl.gca()
            axlist.append(ax)
            if n == 1:
                pl.ylabel('PSF$_{i}$')

            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticklabels([])

            ylims = pl.ylim()
            pl.ylim(-0.001, ylims[1])
            pl.title('zone %d' % n)
            
            pltindex += 1

        ylimit = 0.0
        for n in nodegroups.groups.keys():
            nodelsf = nodegroups.get_group(n)
            prevlsf = prevgroups.get_group(n)

            pl.subplot(3, numnodes, pltindex, sharex=axlist[n-1])
            pl.plot(nodelsf['dj'], nodelsf['lsf']-prevlsf['lsf'], lw=2, color=colors[i-1])

            ax = pl.gca()
            if n == 1:
                pl.ylabel('PSF$_{i}$ - PSF$_{i-1}$')


            yr = np.array(pl.ylim())
            pl.ylim(-np.max(np.abs(yr)), np.max(np.abs(yr)))
            
                                
            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticklabels([])

            pltindex += 1

        for n in nodegroups.groups.keys():
            pl.subplot(3, numnodes, pltindex)
            
            pl.plot(centroids[n-1], i, 'o', color=colors[i-1], markersize=8)
            
            ax = pl.gca()
            axlist.append(ax)
            if n == 1:
                pl.ylabel('Iter [$i$]')

            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticklabels([])
            
            if i == iters[-1]:
                ar = pl.xlim()[1] - pl.xlim()[0]
                vsize = 0.0
                prec = 1.    # for m/s
                while vsize == 0:
                    arvel = ar * 1000.
                    vsize = np.round(arvel*prec / 5) / prec
                    psize = vsize / 1000.
                    label = "%s m s$^{-1}$" % vsize
                    prec *= 10.
    
                #print psize, label
                
                add_scalebar(ax, hidex=False, hidey=False, matchy=False, matchx=False, sizey=0, sizex=psize, labelx=label, loc=3, sep=5)
            
            pl.ylim(0,max(iters)+1)
            pl.xlabel('centroid$_{i}$')
            
            pltindex += 1

        prevgroups = nodegroups

    pl.annotate('Pixel Offset', xy=(0.5, 0.03), xycoords='figure fraction', horizontalalignment='center', fontsize=24)
    pl.suptitle('%s: order: %d, observation index: %d' % (runname, order, iobs))

    return pl.gcf()


def compare_wls_byorder(runname, obdf, datadir, orders=[2,3,4,5,6,7,8,9,10]):
    """

    Plot the input "truth" wavelength solution for all observations and orders from a single ``iGrand`` iteration.

    Args:
        runname (string): name of the iGrand run (e.g. iGrand_sun or iGrand_4628)
        obdf (DataFrame): observation list data frame as output by ``grandsol.io.read_obslist``
        datadir (string): path to directory containing data files
        orders (list): list of orders to combine to derive the RVs for each iteration

    Returns:
        A multi-page PDF file named `runname`_wls_byorder.pdf

    """
    
    obdf.set_index('ind', inplace=True)
    
    colors = [ cmap(x) for x in np.linspace(0.05, 0.95, len(obdf.index))]

    with PdfPages('%s_wls_byorder.pdf' % (runname)) as pdf:
        for o in orders:
            fig = pl.figure(figsize=default_size)

            gout = grandsol.io.read_modfile('%s.%02d.99.mod' % (runname,o))
            for ind in obdf.index:
                obname = obdf.loc[ind, 'obs']
                
                truthfile = os.path.join(datadir,"%s.%02d" % (obname, o))
                tdf = grandsol.io.read_truth(truthfile)

                

                diff = gout[gout['ind'] == ind]['wav_obs'].values - tdf['waveleng'].values

                pl.plot(tdf['waveleng'], diff, '-', lw=1, color=colors[ind-1])
                pl.xlim(min(tdf['waveleng']), max(tdf['waveleng']))

                pl.xlabel('input wavelength [$\\AA$]')
                pl.ylabel('output - input wavelength [$\\AA$]')
                pl.title('%s, order: %02d' % (runname,o))

            pdf.savefig()
            pl.close()

    
    
