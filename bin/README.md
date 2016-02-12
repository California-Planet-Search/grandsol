## iGrand
Driver for parallel iterative runs of grand
```
usage: iGrand [-h] [-o] [--orders ORDERS] [--local] [--ncpus NCPUS]
              [--niter NITER]
              [stars [stars ...]] sysvel

Run iterative Grand Solution doppler analysis.

positional arguments:
  stars            Star or list of stars to analyze.
  sysvel           Systemic radial velocity to star.

optional arguments:
  -h, --help       show this help message and exit
  -o, --overwrite  Overwrite previous run [False]
  --orders ORDERS  List of spectral orders to analyze.
  --local          Run on local computer only. [True]
  --ncpus NCPUS    Number of local CPUs to utilize. [20]
  --niter NITER    Number of outer loop iterations.
```
