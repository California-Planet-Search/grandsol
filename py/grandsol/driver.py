import subprocess
import os
import shutil
import grandsol
import time
import pandas as pd

def execute(cmd, cwd, verbose=True):
    """
    A thin wrapper for subprocess.Popen
    """
    os.chdir(cwd)
    o = int(cmd.split()[3])
    if verbose: print "Running command: '%s'" % cmd
    lock = subprocess.Popen(["echo 'Start Time: '`date` > order_%02d.run" % o], shell=True)

    p = subprocess.Popen(cmd.split())
    stdout, stderr = p.communicate()
    
    lock = subprocess.Popen(["echo 'End Time: '`date` >> order_%02d.run" % o], shell=True)
    time.sleep(2)
    lock = subprocess.Popen(["mv order_%02d.run order_%02d.done" % (o,o)], shell=True)
    
def run_orders(runname, obslist, ppserver=None, overwrite=False, orders=[1,2,3,4,5,6,7,8,9,10,11,12]):
    jobs = []
    for o in orders:
        cwd = os.getcwd()
        cmd = "grand %s %s %d 111111 out=%s.%02d.log vorb+ fudge- nitf=10" % (obslist, runname, o, runname, o)
        if overwrite or not os.path.isfile('order_%02d.done' % o): 
            if ppserver == None:
                execute(cmd, cwd)
            else:
                jobs.append(ppserver.submit(execute, (cmd,cwd), modules=('subprocess','os', 'time')))

    for j in jobs:
        j()
        
    ppserver.wait()
        
def run_iterations(opt, ppserver=None):
    if opt.obslist != None:
        f = open(opt.obslist, 'r')
        for l in f.readlines():
            if l.startswith('RJDIR'): datadir = l.split('=')[1].strip().replace('"','')
        f.close()
        df = pd.read_csv(opt.obslist, sep=' ', skipinitialspace=True, skiprows=2, names=['ind', 'obs', 'unused', 'bc', 'vorb'])
        df['jd'] = df['ind'] + 15000.
    else:
        df = grandsol.io.get_observations(opt.star)
        datadir = os.environ['GRAND_DATADIR']
    runname = "iGrand_" + opt.star
    rundir = os.getcwd()
    runorders = opt.orders
    
    for i in range(opt.niter):
        n = i+1
        idir = "iter%02d" % n
        if not os.path.exists(idir): os.makedirs(idir)
        os.chdir(idir)
        obfile = 'obslist_%02d' % n
        if i == 0:
            vorb = 0.0
            obdf = grandsol.io.write_obslist(df, opt.sysvel, datadir, outfile=obfile, vorb=vorb)
        else:
            vorb = (grandsol.relativity.RV( vel=vdf['mnvel'] ) + grandsol.relativity.RV( vel=vorb )).values()
            obdf = grandsol.io.write_obslist(df, opt.sysvel, datadir, outfile=obfile, vorb=vorb)

        run_orders(runname, obfile, ppserver, orders=runorders, overwrite=opt.overwrite)
        vdf, mnvel = grandsol.io.combine_orders(runname, obdf, runorders, varr_byorder=True)

        # Don't continue running orders that produce 0.0s
        goodorders = []
        for i in range(mnvel.shape[0]):
            if not (mnvel[i,:] == 0).all(): goodorders.append(runorders[i])
        runorders = goodorders
                
        os.chdir(rundir)
