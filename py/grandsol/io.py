import pandas as pd
import numpy as np
import os
from scipy.constants import c
import grandsol.relativity as relativity
import grandsol

kbc = pd.read_csv(os.environ['GRAND_KBCVEL'], sep=' ', skiprows=1, skip_blank_lines=True, skipinitialspace=True, names=['obs', 'name', 'bc', 'jd', 'ha', 'type'])

def get_observations(star, thin=999):
    """
    Find observations of a given star from $GRAND_KBCVEL and convert into a Pandas DataFrame.

    Parameters
    ----------
    star : str : star name (case insensitive)


    Returns
    -------
    Pandas DataFrame : containing the rows in $GRAND_KBCVEL corresponding to the requested star

    """
    
    star = kbc[(kbc['name'].str.upper() == star.upper()) & (kbc['type'] == 'o') & kbc.obs.str.startswith('rj')]

    if len(star.obs.values) > thin:
        np.random.seed(0)
        perm = np.random.permutation(len(star.obs.values))[:thin]
        star = star.iloc[perm].sort_values('jd').reset_index()
        
    return pd.DataFrame(star)

def read_obslist(infile):
    df = pd.read_csv(infile, skiprows=2, sep=' ', skipinitialspace=True, 
                      names=['ind', 'obs', 'unused', 'bc', 'vorb'])
    odf = pd.merge(df, kbc, on='obs')

    if odf.empty:
        odf = df
        odf['jd'] = odf.ind + 15500.   # Hack to fake JD timestamps in simulated data
        
    return odf

def write_obslist(df, sysvel, datadir, outfile='obslist', vorb=None, overwrite=False):
    """
    Write list of observations in the format that grand likes.

    Parameters
    ----------
    df : DataFrame : Pandas DataFrame with observation rows to be analyzed
    sysvel : float : Systemic radial velocity of the system in m/s
    outfile : str : (optional) name of output observation list
    vorb : scalar or array : initial velocity guess for each observation

    Returns
    -------
    DataFrame : Pandas DataFrame that looks like the output obslist
    """
    

    df['fill'] = 0
    if 'vorb' not in df.columns or isinstance(vorb, type(None)):
        df['vorb'] = 0
    else:
        df['vorb'] = vorb.values

    odf = df.sort_values('jd').reset_index(drop=True)
    odf['ind'] = odf.index.values + 1
    
    header = 'VSYST = %.0f m/s\nRJDIR = "%s/"\n' % (sysvel, datadir)
    body = odf.to_string(index=False, header=False,
                         columns=['ind', 'obs','fill', 'bc', 'vorb'],
                         formatters=['{:03d}'.format, '{:s}'.format, '{:d}'.format, '{:.5f}'.format, '{:.5f}'.format])
    if overwrite or not os.path.isfile(outfile):
        f = open(outfile, 'w')
        print >>f, header+body
        f.close()

    return odf

def read_vel(infile):
    """
    Read velocity file (.vel) from grand output

    Parameters
    ----------
    infile : str : input file name

    Returns
    ---------
    DataFrame : Pandas DataFrame with the data contained in the specified velocity file.
    Z is converted to velocity and appended as new columns where appropriate.
    """
    
    zdf = pd.read_csv(infile, sep=' ', names=['ind', 'zn', 'z0', 'zbarn', 'zsign']).set_index('ind')

    zdf['vbarn'],zdf['uvbarn'] = relativity.redshift_to_vel(zdf['zbarn'], zdf['zsign'])
    zdf['veln'] = relativity.z2v(zdf['zn'])
    zdf['bc'] = relativity.z2v(zdf['z0'])
    
    return zdf

def combine_orders(runname, obdf, orders, varr_byorder=False, usevln=False):
    """
    Combine velocities from multiple orders by mean and merge with observation information.

    Parameters
    ----------
    runname : str : name of current grand run
    obdf : DataFrame : observation data frame from grandsol.io.write_obslist
    orders : list : list of orders to combine
    varr_byorder : bool : (optional) return full velocity-by-order array in addition to the normal output
    usevln : bool : (optional) use zln instead of zbarn to calculate velocities

    Returns
    ---------
    vdf : DataFrame : Same as obdf with mean velocity (mnvel) and velocity uncertainty (errvel) columns added
    or
    (vdf, mnvel) : tuple : mnvel is a len(orders) X len(observations) array that contains the velocities
                           for each other before taking the mean
    """
    
    mnvel = []
    zarr = []
    for i,o in enumerate(orders):
        vdf = grandsol.io.read_vel('%s.%02d.99.vel' % (runname,o))

        if (vdf['zn'] == vdf['z0']).all():
            print "io.combine_orders: WARNING: order %d velocities are all 0.0" % o
            continue

        if usevln:
            mnvel.append(vdf['veln'].values)
            zarr.append(vdf['zn'].values)
        else:
            mnvel.append(vdf['vbarn'].values)
            zarr.append(vdf['zbarn'].values)

    mnvel = np.atleast_1d(np.vstack(mnvel))
    zarr = np.atleast_1d(np.vstack(zarr))
            
    rv = relativity.RV(z=zarr)
    bc = relativity.RV(z=vdf['z0'].values)
    prev = relativity.RV(vel=obdf['vorb'].values)

    relvel = (rv - bc) + prev
        
    vdf['mnvel'] = relvel.mean().values()
    vdf['mnvel'] -= vdf.mnvel.mean()
    vdf['errvel'] = relvel.values().std(axis=0) / np.sqrt(mnvel.shape[0])

    mdf = pd.merge(vdf, obdf, left_index=True, right_on='ind')

    if varr_byorder: return (mdf, relvel.vel)
    else: return mdf


def read_modfile(modfile):
    model = pd.read_csv(modfile, sep=' ', skipinitialspace=True,
                        names=['ind', 'order', 'pixel', 'spec', 'model', 'wav_obs', 'wav_star', 'cont',
                                'smooth_cont', 'badflag', 'tellflag', 'metflag'])

    return model
