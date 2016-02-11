import pandas as pd
import os
import grandsol.relativity as relativity

def get_observations(star):
    """
    Find observations of a given star from $GRAND_KBCVEL and convert into a Pandas DataFrame.

    Parameters
    ----------
    star : str : star name (case insensitive)


    Returns
    -------
    Pandas DataFrame : containing the rows in $GRAND_KBCVEL corresponding to the requested star

    """
    
    kbc = pd.read_csv(os.environ['GRAND_KBCVEL'], sep=' ', skiprows=1, skip_blank_lines=True, skipinitialspace=True, names=['obs', 'name', 'bc', 'jd', 'ha', 'type'])
    star = kbc[(kbc['name'].str.upper() == star.upper()) & (kbc['type'] == 'o')]
    return pd.DataFrame(star)

def write_obslist(df, sysvel, outfile='obslist', vorb=0):
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
    
    f = open(outfile, 'w')

    df['fill'] = 0
    if 'vorb' not in df.columns:
        df['vorb'] = vorb

    odf = df.sort_values('jd').reset_index(drop=True)
    odf['ind'] = odf.index.values + 1
        
    header = 'VSYST = -10140 m/s\nRJDIR = "%s/"\n' % (os.environ['GRAND_DATADIR'])
    body = odf.to_string(index=False, header=False,
                         columns=['ind', 'obs','fill', 'bc', 'vorb'],
                         formatters=['{:03d}'.format, '{:s}'.format, '{:d}'.format, '{:.5f}'.format, '{:.5f}'.format])

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

    
