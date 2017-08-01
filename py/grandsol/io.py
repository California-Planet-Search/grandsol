import pandas as pd
import numpy as np
import os
from scipy.constants import c
from scipy.interpolate import interp1d
import grandsol.relativity as relativity
import grandsol.rdsk as rdsk
import grandsol

def check_snr(obdf, ctslim=500, verbose=False):
    """Check SNR for each observation in an obslist DataFrame

    Check SNR for each observation in an obslist DataFrame.

    Args:
        obdf (DataFrame): obslist DataFrame
        ctslim (float): (optional) Minimum limit in counts calculated as 
            80th percentile of the mean flux over each order
    
    Returns:
        DataFrame: stripped of low SNR observations
    """
    
    datadir = os.environ['GRAND_DATADIR']

    good = []
    for i,row in obdf.iterrows():
        ob = row['obs']
        fpath = os.path.join(datadir, ob)
        d = rdsk.Dskfile(fpath)
        specarr = d.record(0)
        order_means = specarr.mean(axis=1)
        cts = np.percentile(order_means, 80)
        if cts >= ctslim:
            good.append(i)
        if verbose:
            print("{} {} {} {}".format(i, ob, fpath, cts))

    return obdf.loc[good]

def load_bc(bcfile):
    """
    Load kbcvel.ascii or abcvel.ascci and return as
    Pandas DataFrame.

    Args:
        kbcfile (string): name of bary file

    Return:
        DataFrame
    """
    
    kbc = pd.read_csv(bcfile, sep=' ', skiprows=1, skip_blank_lines=True,
                          skipinitialspace=True,
                          names=['obs', 'name', 'bc', 'jd', 'ha', 'type'])

    return kbc

def get_observations(star, thin=999):
    """
    Find observations of a given star from $GRAND_KBCVEL and 
    convert into a Pandas DataFrame.

    Args:
        star (str): star name (case insensitive)
        thin (int): limit maximum number of observations

    Returns:
       DataFrame containing the rows in $GRAND_KBCVEL corresponding 
           to the requested star
        
    """
    bcfile = os.environ['GRAND_KBCVEL']
    kbc = load_bc(bcfile)
    
    star = kbc[(kbc['name'].str.upper() == star.upper()) &
                   (kbc['type'] == 'o') & ~kbc.obs.str.startswith('rk')]

    if len(star.obs.values) > thin:
        np.random.seed(0)
        perm = np.random.permutation(len(star.obs.values))[:thin]
        star = star.iloc[perm].sort_values('jd').reset_index()

    sdf = pd.DataFrame(star)
    sdf = check_snr(sdf)
        
    return sdf

def read_obslist(infile):
    """

    Read an obslist file and translate it into a Pandas DataFrame.

    Args:
        infile (string): name of input obslist file

    Returns:
        DataFrame corresponding to input obslist.

    """
    
    df = pd.read_csv(infile, skiprows=3, sep=' ', skipinitialspace=True, 
                      names=['ind', 'obs', 'unused', 'bc', 'vorb'])

    baryfile = os.environ['GRAND_KBCVEL']
    kbc = load_bc(baryfile)
    odf = pd.merge(df, kbc, on='obs', suffixes=['','_obl'])
    
    if odf.empty:
        odf = df
        odf['jd'] = odf.ind + 15500.   # Hack to fake JD timestamps in simulated data
        
    return odf

def write_obslist(df, sysvel, datadir, outfile='obslist',
                      vorb=None, meteor=False, inst="HIRES", overwrite=False):
    """
    Write list of observations in the format that ``grand`` likes.

    Args:
        df (DataFrame): Pandas DataFrame with observation rows to be analyzed
        sysvel (float): Systemic radial velocity of the system in m/s
        outfile (string): (optional) name of output observation list
        vorb (float scalar or array): (optional) initial velocity 
            guess for each observation
        meteor (bool): Add meteor line to obslist?
        instrument (string): instrument string in obslist
        overwrite (bool): (optional) overwrite if outfile already exists?
        
    Returns:
        Pandas DataFrame that looks like the output obslist
    """
    

    df['fill'] = 0
    if 'vorb' not in df.columns or isinstance(vorb, type(None)):
        df['vorb'] = 0
    else:
        df['vorb'] = vorb.values

    odf = df.sort_values('jd').reset_index(drop=True)
    odf['ind'] = odf.index.values + 1
    
    header = 'VSYST = %.0f m/s\nRJDIR = "%s/"\nINSTRUMENT = "%s"\n' \
               % (sysvel, datadir, inst)
    if meteor:
        header += "METEOR meteor\n" 
    body = odf.to_string(index=False, header=False,
                         columns=['ind', 'obs','fill', 'bc', 'vorb'],
                         formatters=['{:03d}'.format, '{:s}'.format,
                                    '{:d}'.format, '{:.5f}'.format,
                                    '{:.5f}'.format])
    if overwrite or not os.path.isfile(outfile):
        f = open(outfile, 'w')
        print >>f, header+body
        f.close()

    return odf

def read_vel(infile):
    """
    Read velocity file (.vel) from grand output. Z is converted 
    to velocity and appended as new columns where appropriate.

    Args:
        infile (string): input file name

    Returns:
        DataFrame: Pandas DataFrame with the data contained in 
            the specified velocity file.
    
    """
    
    zdf = pd.read_csv(infile, sep=' ',
                names=['ind', 'zn', 'z0', 'zbarn', 'zsign']).set_index('ind')

    lines_find = np.loadtxt(infile.replace('.99.vel','.lines_find'))
    nlines = lines_find.shape[0]
    
    zdf['vbarn'],zdf['uvbarn'] = relativity.redshift_to_vel(zdf['zbarn'],
                                                            zdf['zsign'])

    zdf['uvbarn'] *= c / np.sqrt(nlines)
    
    zdf['veln'] = relativity.z2v(zdf['zn'])
    zdf['bc'] = relativity.z2v(zdf['z0'])
    
    return zdf

def combine_orders(runname, obdf, orders, varr_byorder=False,
                       usevln=True, get_weights=False, rv_fudge=True):
    """
    Combine velocities from multiple orders by mean and merge 
    with observation information.

    Args:
        runname (str): name of current grand run
        obdf (DataFrame): observation data frame from 
            grandsol.io.write_obslist
        orders (list): list of orders to combine
        varr_byorder (bool): (optional) return full velocity-by-order
            array in addition to the normal output
        usevln (bool): (optional) use zln instead of zbarn to 
            calculate velocities
        get_weights (bool): (optional) return the weight matrix 
            instead of the velocities
        rv_fudge (bool): (optional) apply the RV fudge factor

    Returns:
        DataFrame: Same as obdf with mean velocity (mnvel) and 
            velocity uncertainty (errvel) columns added
        or
        tuple: (vdf, mnvel) mnvel is a len(orders)xlen(observations) 
            array that contains the velocities
            for each other before taking the mean
        or
        array: weight matrix if get_weights==True)
    """
    
    mnvel = []
    zarr = []
    warr = []
    for i,o in enumerate(orders):
        fname = '%s.%02d.99.vel' % (runname,o)
        vdf = grandsol.io.read_vel(fname)

        if (vdf['zn'] == vdf['z0']).all():
            print "io.combine_orders: WARNING: order %d velocities are all 0.0"\
                     % o
            continue

        if usevln:
            mnvel.append(vdf['veln'].values)
            zarr.append(vdf['zn'].values)
        else:
            mnvel.append(vdf['vbarn'].values)
            zarr.append(vdf['zbarn'].values)

        warr.append(vdf['uvbarn'])

    mnvel = np.atleast_1d(np.vstack(mnvel))
    zarr = np.atleast_1d(np.vstack(zarr))
    warr = np.atleast_1d(np.vstack(warr))
            
    rv = relativity.RV(z=zarr)
    bc = relativity.RV(z=vdf['z0'].values)
    prev = relativity.RV(vel=obdf['vorb'].values)
    err = relativity.RV(vel=vdf['uvbarn'].values)
    
    relvel = (rv - bc) + prev

    #w = 1/warr.mean(axis=1)
    w = 1/np.std(relvel.vel, axis=1)
    w /= w.max()
    #w = w*0 + 1.0
    
    #vdf['mnvel'] = relvel.mean().values()
    #vdf['mnvel'] = np.average(relvel.values(), weights=w, axis=0)
    #weight_matrix = np.ones_like(relvel)
    #vdf['mnvel'] = np.median(relvel.values(), axis=0)
    vdf['mnvel'], weight_matrix = grandsol.utils.clipped_mean(relvel.values(),
                                                        inweights=w, sigma=4)
    vdf['absvel'], weight_matrix = grandsol.utils.clipped_mean(rv.vel,
                                                        inweights=w, sigma=4)

    if rv_fudge:
        x, y = np.genfromtxt(os.environ['GRAND']+'/py/rv_fudge_apf.csv',
                                 delimiter=',', unpack=True)
        absvel = vdf['absvel'].values
        f = interp1d(x, y, kind='linear',
                            fill_value=0.0,
                            bounds_error=False)
        corr = f(absvel)
        vdf['mnvel'] -= corr
    
    vdf['mnvel'] -= vdf['mnvel'].mean()
    vdf['errvel'] = relvel.values().std(axis=0) / np.sqrt(mnvel.shape[0])
    if (vdf['errvel'] == 0).all():
        vdf['errvel'] = err.vel

    
        
    obdf.ind = np.array(obdf.ind.values, dtype=int)
    mdf = pd.merge(vdf, obdf, left_index=True, right_on='ind')
    

    if varr_byorder: return (mdf, relvel.vel)
    elif get_weights: return weight_matrix
    else: return mdf


def read_modfile(modfile):
    """
    Read in a .mod file from a ``grand`` run


    Args:
        modfile (string): name of file

    Returns:
        Pandas DataFrame view of the .mod file
    """

    model = pd.read_csv(modfile, sep=' ', skipinitialspace=True,
                        names=['ind', 'order', 'pixel', 'spec',
                                'model', 'wav_obs', 'wav_star', 'cont',
                                'smooth_cont', 'badflag', 'tellflag',
                                'metflag'])

    return model

def read_temfile(temfile):
    """
    Read in a .mod file from a ``grand`` run

    Args:
        temfile (string): name of file

    Returns:
        Pandas DataFrame view of the .tem file
    """

    temp = pd.read_csv(temfile, sep=' ', skipinitialspace=True,
                           names=['inode', 'wav', 'temp', 'temp_prev', 'solar'])
    
    temp.drop_duplicates(subset=['temp', 'temp_prev'],
                             inplace=True, keep=False)
    
    temp.temp *= temp.solar.mean()
    temp.temp_prev *= temp.solar.mean()
    
    return temp

def read_grlsf(lsffile):
    """
    Convert ``grlsf`` output containing PSFs for all nodes 
    accross a single order into a Pandas DataFrame.

    Args:
        lsffile (any object with a read method): Name 
            of the lsffile,
            or any file buffer like object with a read() method 
            (e.g. stdout from a subprocess.Popen call).
            This must be from a multi-node ``grlsf`` call (final argument = 0)

    Returns:
        DataFrame corresponding to ``grlsf`` output with columns labeled.

    """

    df = pd.read_csv(lsffile, sep=' ', skipinitialspace=True,
                         comment="#",
                         names=['obs', 'order', 'node', 'j', 'dj', 'lsf'])

    return df

def read_truth(tfile):
    """
    Read an input "truth" file from a simulated dataset.

    Args:
        wlsfile (string): input file name, (e.g. sun.034.02)

    Returns:
        DataFrame: input data contained within `wlsfile` 
            represented as a Pandas DataFrame

    """

    df = pd.read_csv(tfile, sep=' ', skipinitialspace=True,
                         comment='#',
                         names=['waveleng', 'temp*iod',
                                'noise', 'norm_vec', 'obs_spec'])

    return df

def write_velocities(df, outfile):
    """
    Write velocities to a file

    Args:
        df (DataFrame): DataFrame output from `grandsol.io.combine_orders`
        outfile (string): Name of output file

    """

    colmap = {'bc_x': 'bc'}
    for old,new in colmap.items():
        df[new] = df[old]
        
    df.to_csv(outfile, sep=' ',
              columns=['obs', 'jd', 'mnvel', 'errvel', 'bc'],
              index=False)

