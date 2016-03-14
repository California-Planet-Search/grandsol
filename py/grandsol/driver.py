import subprocess
import os
import shutil
import grandsol
import time
import pandas as pd
import numpy as np

def execute(cmd, cwd):
    """
    A thin wrapper for subprocess.Popen
    """
    os.chdir(cwd)
    o = int(cmd.split()[3])
    lock = subprocess.Popen(["echo 'Start Time: '`date` > order_%02d.run" % o], shell=True)

    p = subprocess.Popen(cmd.split())
    stdout, stderr = p.communicate()
    errcode = p.returncode
    
    lock = subprocess.Popen(["echo 'End Time: '`date`'\nreturn code: %s' >> order_%02d.run" % (errcode,o)], shell=True)
    time.sleep(2)
    lock = subprocess.Popen(["mv order_%02d.run order_%02d.done" % (o,o)], shell=True)

    return errcode
    
def run_orders(runname, obslist, ppserver=None, overwrite=False, fudge=True, orders=[1,2,3,4,5,6,7,8,9,10,11,12]):
    jobs = []
    good_orders = []
    for o in orders:
        cwd = os.getcwd()
        cmd = "grand %s %s %d 111111 out=%s.%02d.log vorb+ nitf=10" % (obslist, runname, o, runname, o)
        if fudge: cmd += " fudge+"
        else: cmd += " fudge-"
        if overwrite or not os.path.isfile('order_%02d.done' % o): 
            if ppserver == None:
                execute(cmd, cwd)
            else:
                print "Running command: '%s'" % cmd
                jobs.append(ppserver.submit(execute, (cmd,cwd), modules=('subprocess','os', 'time')))
        else:
            f = open('order_%02d.done' % o, 'r')
            for l in f.readlines():
                if l.startswith("return"):
                    errcode = int(l.split(":")[-1])
                    break
            f.close()
            if errcode == 0: good_orders.append(o)
                
    for o,job in zip(orders,jobs):
        e = job()
        if e == 0: good_orders.append(o)
        
    ppserver.wait()

    return good_orders
            
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

    iterdone = []
    for i in range(opt.niter):
        n = i+1
        idir = "iter%02d" % n
        if not os.path.exists(idir): os.makedirs(idir)
        os.chdir(idir)
        obfile = 'obslist_%02d' % n
        if i == 0:
            vorb = pd.Series(np.zeros_like(df['jd']))
            obdf = grandsol.io.write_obslist(df, opt.sysvel, datadir, outfile=obfile, vorb=vorb)
        else:
            vorb = vdf['mnvel']
            obdf = grandsol.io.write_obslist(df, opt.sysvel, datadir, outfile=obfile, vorb=vorb)

        runorders = run_orders(runname, obfile, ppserver, orders=runorders, overwrite=opt.overwrite, fudge=opt.fudge)
        vdf, mnvel = grandsol.io.combine_orders(runname, obdf, runorders, varr_byorder=True)

        grandsol.plotting.velplot_by_order(runname, obdf, runorders, outfile='iGrand_%s_velbyord.pdf' % opt.star)
        grandsol.plotting.velplot_by_order(runname, obdf, runorders, outfile='iGrand_%s_bcbyord.pdf' % opt.star, vsbc=True)

        iterdone.append(n)
        
        os.chdir(rundir)

        grandsol.plotting.velplot_by_iter(runname, runorders, outfile='iGrand_%s_velbyiter.pdf' % opt.star, iters=iterdone)
