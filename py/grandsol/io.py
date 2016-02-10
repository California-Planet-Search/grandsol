import pandas as pd
import cpsutils
import grandsol.relativity as relativity

def get_observations(star):
    kbc = cpsutils.kbc.loadkbc(mode='kbc', quiet=True)
    star = kbc[(kbc.name.str.upper() == star.upper()) & (kbc['type'] == 'o')]
    return star

def write_obslist(df, outfile='obslist'):
    f = open(outfile, 'w')

    df['fill'] = 0
    if 'vorb' not in df.columns:
        df['vorb'] = 0

    odf = df.sort_values('jd').reset_index(drop=True)
    odf['ind'] = odf.index.values + 1
        
    header = 'VSYST = -10140 m/s\nRJDIR = "/mir3/iodspec/"\n'
    body = odf.to_string(index=False, header=False,
                         columns=['ind', 'obs','fill', 'bc', 'vorb'],
                         formatters=['{:03d}'.format, '{:s}'.format, '{:d}'.format, '{:.5f}'.format, '{:.5f}'.format])

    print >>f, header+body
    f.close()

    return odf

def read_vel(infile):
    zdf = pd.read_csv(infile, sep=' ', names=['ind', 'zn', 'z0', 'zbarn', 'zsign']).set_index('ind')

    zdf['vbarn'],zdf['uvbarn'] = relativity.redshift_to_vel(zdf['zbarn'], zdf['zsign'])
    zdf['veln'] = relativity.z2v(zdf['zn'])
    zdf['bc'] = relativity.z2v(zdf['z0'])
    
    return zdf

    
