#!/usr/bin/env python

import numpy as np
import pylab as pl
import pandas as pd
import os
import sys
import subprocess
import argparse

import grandsol


def args():
    """
    Parse command line arguments.
    """
    
    parser = argparse.ArgumentParser(description='Plot LSFs output by grand')
    parser.add_argument(metavar='runname',dest='runname',action='store',
                        help='Name of the iGrand run (e.g. iGrand_sun or iGrand_4628)', type=str)
    parser.add_argument(metavar='obs',dest='obs',action='store',
                        help='observation index from obslist (starting with 1, not 0)', type=int)
    parser.add_argument(metavar='order',dest='order',action='store',
                        help='order index (starting with 1, not 0)', type=int)
    parser.add_argument('-s','--show',dest='show',action='store_true',
                        help='display the interactive plotting window')


    opt = parser.parse_args()

    return opt
    
def main():
    opt = args()

    fig = grandsol.plotting.plot_lsf_byiter(opt.runname, opt.obs, opt.order)

    pl.savefig("%s_%03d_%02d_lsf_byiter.png" % (opt.runname, opt.obs, opt.order))
    if opt.show: pl.show()

if __name__ == '__main__':
    main()
