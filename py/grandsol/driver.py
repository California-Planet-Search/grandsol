import subprocess
import os
import shutil
import grandsol

def execute(cmd, cwd, verbose=True):
    """
    A thin wrapper for subprocess.Popen
    """
    os.chdir(cwd)
    if verbose: print "Running command: '%s'" % cmd
    p = subprocess.Popen(cmd.split())
    stdout, stderr = p.communicate()
    
def run_orders(runname, obslist, ppserver=None, orders=[3,4,5,6,7,8,9]):
    jobs = []
    for o in orders:
        cwd = os.getcwd()
        cmd = "grand %s %s %d 111111 out=%s.%02d.log" % (obslist, runname, o, runname, o)
        if ppserver == None:
            execute(cmd, cwd)
        else:
            jobs.append(ppserver.submit(execute, (cmd,cwd), modules=('subprocess','os')))

    for j in jobs:
        j()
        
    ppserver.wait()
        
def run_iterations(opt, ppserver=None):
    df = grandsol.io.get_observations(opt.star)
    runname = "iGrand_" + opt.star
    idir = os.getcwd()
    
    for i in range(opt.niter):
        n = i+1
        idir = "iter%02d" % n
        if not os.path.exists(idir): os.makedirs(idir)
        os.chdir(idir)
        obfile = 'obslist_%02d' % n
        if i == 0:
            obdf = grandsol.io.write_obslist(df, opt.sysvel, outfile=obfile, vorb=0)
        else:
            obdf = grandsol.io.write_obslist(df, opt.sysvel, outfile=obfile, vorb=vdf['mnvel'].values)  # Need to update vorb here

        run_orders(runname, obfile, ppserver, orders=opt.orders)
        vdf = grandsol.io.combine_orders(runname, obdf, opt.orders)
        
        os.chdir(idir)
