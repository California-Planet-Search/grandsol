import subprocess
import os
import shutil
import grandsol
import time
import pandas as pd
import numpy as np

def execute(cmd, cwd, plotres):
    """
    A thin wrapper for subprocess.Popen

    Args:
        cmd (string): Shell command to execute
        cwd (string): Directory where the command will be executed
        plotres (bool): Make the residual plots?

    Returns:
        int: error code returned after command execution
    
    """
    os.chdir(cwd)
    o = int(cmd.split()[3])
    runname = cmd.split()[2]
    lock = subprocess.Popen(["echo 'Start Time: '`date` > order_%02d.run" % o], shell=True)

    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    errcode = p.returncode

    modfile = "%s.%02d.99.mod" % (runname, o)
    if os.path.isfile(modfile) and plotres: grandsol.plotting.plot_residuals_byobs(modfile, outfile="%s_%02d_residuals.png" % (runname, o))
    
    lock = subprocess.Popen(["echo 'End Time: '`date`'\nreturn code: %s' >> order_%02d.run" % (errcode,o)], shell=True)
    time.sleep(2)
    lock = subprocess.Popen(["mv order_%02d.run order_%02d.done" % (o,o)], shell=True)

    return errcode
    
def run_orders(runname, obslist, ppserver=None, overwrite=False, fudge=True, orders=[1,2,3,4,5,6,7,8,9,10,11,12], plotres=False):
    """
    Run ``grand`` for several orders simultaneously.

    Args:
        runname (string): name of current run
        obslist (string): name of obslist file
        ppserver (pp.Server): Parallel Python server object to send jobs
        overwrite (bool): overwrite previous run?
        fudge (bool): Apply the fudge factor?
        orders (list): list of orders to be run
        plotres (bool): make the residual plots?

    Returns:
        tuple: orders that succesfully finished, and the list of jobs sent

    """
    
    
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
                jobs.append(ppserver.submit(execute, (cmd,cwd,plotres), modules=('subprocess','os', 'time', 'grandsol')))
        else:
            f = open('order_%02d.done' % o, 'r')
            for l in f.readlines():
                if l.startswith("return"):
                    errcode = int(l.split(":")[-1])
                    break
            f.close()
            if errcode == 0:
                good_orders.append(o)
                modfile = "%s.%02d.99.mod" % (runname, o)
                if os.path.isfile(modfile) and plotres: grandsol.plotting.plot_residuals_byobs(modfile, outfile="%s_%02d_residuals.png" % (runname, o))

                
    for o,job in zip(orders,jobs):
        e = job()
        if e == 0: good_orders.append(o)
        
    ppserver.wait()            

    return good_orders, jobs
            
def run_iterations(opt, ppserver=None):
    """
    Run ``iGrand`` iterations

    Args:
        opt (argparse.ArgumentParser): command line arguments object from argparse
        ppserver (pp.Server): Parallel Python server object where jobs will be sent

    Returns:
        None

    """

    
    if opt.obslist != None:
        f = open(opt.obslist, 'r')
        for l in f.readlines():
            if l.startswith('RJDIR'): datadir = l.split('=')[1].strip().replace('"','')
            if l.startswith('VSYST'): opt.sysvel = float(l.split('=')[1].split()[0])
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

        runorders, joblist = run_orders(runname, obfile, ppserver, orders=runorders, overwrite=opt.overwrite, fudge=opt.fudge, plotres=opt.plotres)
        vdf, mnvel = grandsol.io.combine_orders(runname, obdf, runorders, varr_byorder=True)

        grandsol.plotting.velplot_by_order(runname, obdf, runorders, outfile='iGrand_%s_velbyord.pdf' % opt.star)
        grandsol.plotting.velplot_by_order(runname, obdf, runorders, outfile='iGrand_%s_bcbyord.pdf' % opt.star, vsbc=True)
            
        iterdone.append(n)
        
        os.chdir(rundir)
        
        if len(joblist) > 0 or n == opt.niter:
            grandsol.plotting.velplot_by_iter(runname, runorders, outfile='iGrand_%s_velbyiter.pdf' % opt.star, iters=iterdone)
            if opt.truth:
                grandsol.plotting.truthplot(runname, opt.truthvel, runorders, outfile='iGrand_%s_truth.pdf' % opt.star, iters=iterdone)
            if opt.phase != None:
                grandsol.plotting.phaseplot_by_iter(runname, obdf, runorders, opt.phase[1], opt.phase[0], outfile='%s_phase.pdf' % runname, iters=iterdone)
                
        if opt.plottemp: grandsol.plotting.plot_template_byiter(runname, runorders, iters=iterdone)

    if opt.plotres: grandsol.plotting.plot_resMAD_byiter(runname, obdf, runorders, iters=iterdone, outfile="%s_resMAD_byiter.pdf" % runname)

