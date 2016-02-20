This directory contains grand solution software, which deduces stellar
parameters by fitting stellar spectra obtained with an iodine cell.

The fitting code 'grandsol/fort/grand' determines:
  intrinsic stellar spectrum (a.k.a. deconvolved template)
  stellar radial velocity for each epoch
  wavelength solutions for each epoch
  crude continuum normalization for each epoch
  line spread functions across the echelle format for each epoch

Brief descriptions of 'grandsol' subdirectories:
  'fort' - fortran source code for the fitting code 'grand.F'
  'fref' - reference files used by the fortran code 'grand'

See README files in each subdirectory for more information.

Sample environment variables and aliases for tcsh:
  setenv GRAND /Users/valenti/git/grandsol

  setenv GRAND_REFDIR $GRAND/fref

  setenv GRAND_IREFDIR $GRAND/iref
    
  setenv GRAND_DATADIR /mir3/iodspec
  
  setenv GRAND_KBCVEL /mir3/bary/kbcvel.ascii

  setenv GRANDSOL_IDL_PATH ${MINIMAL_IDL_PATH}:+$GRAND/idl

  setenv PYTHONPATH $GRAND/py:$PYTHONPATH  

  alias idlgr "env IDL_PATH=$GRANDSOL_IDL_PATH $IDL_DIR/bin/idl"

  alias grand ~/git/grandsol/fort/grand

