# Fortran code for the Grand Solution Doppler code

## Environment variables

Set GRAND_REFDIR to point at the directory containing reference files used
by 'grand'. The usual reference file directory is 'grandsol/gref'.
See 'grandsol/gref/README' for more information.


## Building `grand`

Run `make grand` in the `grandsol/fort` source directory.

Executable is written to the `grandsol/bin` directory.


## Calling syntax for `grand`
```
usage: grand obslist label m
             [nito=] [nitf=] [out=] [fudge{-,+}] [vorb+]
             [{0,1,2}{0,1,2}{0,1,2}{0,1,2}{0,1,2}{0,1,2}]
             [vel{0,1,2}] [tem{0,1,2}] [mod{0,1,2}]
             [lsf{0,1,2}] [nrm{0,1,2}] [lin{0,1,2}]

Simultaneously fit a set of observed spectra obtained with an iodine cell.
Return radial velocities, deconvolved stellar template, wavelengths,
model fit, line spread functions, and normalization vectors.

positional arguments:
  obslist               file containing list of observed spectra
  label                 root to use when generating output files
  m                     order number (starting with 1, not 0)

optional arguments:
  {0,1,2} x 6           output flags for {vel}{tem}{mod}{lsf}{nrm}{lin}
  nito=                 starting iteration number            (default: 0)
  nitf=                 final iteration number               (default: 10)
  out=                  file to contain diagnostics          (default: stdout)
  vorb{-,+}             read vorb from column 5 of obslist   (default: F)
  fudge{-,+}            apply fudge correction to model      (default: T)
  vel{0,1,2}            output flag for velocities           (default: 2)
  tem{0,1,2}            output flag for template spectrum    (default: 1)
  mod{0,1,2}            output flag for model fit            (default: 0)
  lsf{0,1,2}            output flag for line spread function (default: 1)
  nrm{0,1,2}            output flag for normalization vector (default: 0)
  lin{0,1,2}            output flag for line regions         (default: 0)

notes:
  Output flags values:
     0: do not write output files
     1: write output files for current iteration
     2: write output files for every iteration
  Maximum length of label is 59 characters
  Interpretation of relative order number m for HIRES:
      1: order 71 = 4977-5065 Angstroms
     16: order 56 = 6310-6368 Angstroms
  Output flags can be specified two different ways:
     collectively as a six-digit integer (default: 210100)
     individually as a three-letter tag and a one-digit integer (e.g., mod1)
```


## Input Files

### Observation List (e.g., `obslist`)

#### Header
| **Example** | **Default** | **Description** | **Variable** |
| :--- | :--- | :--- | :--- |
| `VSYST = -16400 m/s` | *Mandatory* | Radial velocity of center of mass | `velsys` |
| `RJDIR = "\Users\valenti\obs\"` | `".\"` | Directory containing input observations | `RJDIR` |
| `METEOR \Users\valenti\ref\meteor.dat` | `NONE` | File describing meteor location | `METEORFILE` |

#### Body
| **Column** | **Description** | **Example** | **Variable** |
| :---: | :--- | :--- | :--- |
| **0** | Observation index (start with 001, not 000) | `001` | `n` |
| **1** | Input _observation file_ | `star.001` | `FILEI` |
| **2** | Unused | `0` | `d` |
| **3** | Barycentric correction (m/s) | `29943.73508` | `barycorr` |
| **4** | Observed radial velocity around barycenter of distant system (m/s) | `-31.74` | `vorb` |

Column 4 is required if `grand` is called with the `vorb+` argument.
Otherwise `vorb` is set to zero, which is the _null hypothesis_.

See fortran subroutine `read_raw()` for more information.

### Observation File (e.g., `star.001`)

An input observation file must be in `.dsk` format,
as produced by the reduction pipeline or `wdsk_simobs.pro`.

Use an _observation list_ to specify a set of inut observation files.

See fortran subroutine `rdsk()` for more information.
This routine will need to be generalized.
Spectra are assumed to consist of 16 echelle orders, each with 4021 pixels.
Orders 1-12 (out of 0-15) are divided by a hard-coded polynomial approximation of the HIRES blaze function.

## Output Files

### Velocity File (e.g. `star.08.10.vel`)

#### Header
None.

#### Body
| **Column** | **Description** | **Example** | **Variable** |
| :---: | :--- | :--- | :--- |
| **0** | Observation index | `001` | `n` |
| **1** | RV/c for bulk shift of the entire order | `0.000073135` | `zn` |
| **2** | BC/c | `0.0000731371` | `z0` |
| **3** | mean of RV/c from spectral lines after sigma clipping | `0.0000731375` | `zbarn` |
| **4** | standard deviation of RV/c from spectral lines | `0.0000000047` | `zsign` |

Symbols:
- RV = radial velocity
- BC = barycentric correction
- c = speed of light
- z = velocity divided by the speed of light

Grand returns two estimates for radial velocity (RV).
Column 1 is the RV/c that best matches the template to an entire observed order.
Column 3 is the sigma-clipped mean of many RV/c values,
each matching the template to a different line feature in the order.
Both estimates seem to work well, so at this point we don't have a preference.

#### Formulae
```
vbarn = c * ((1+zbarn)^2    - 1) / ((1+zbarn)^2    + 1)

V = c * ((1+Z)^2 - 1) / ((1+Z)^2 + 1)

uV^2 = (dV/dZ)^2 * uZ^2

       ((1+Z)^2 + 1) * (2*(1+Z)) - ((1+Z)^2 - 1) * (2*(1+Z))
     = ----------------------------------------------------- * uZ^2
                           ((1+Z)^2 + 1)^2

       2*(1+Z) * ((1+Z)^2 + 1 - (1+Z)^2 + 1)
     = ------------------------------------- * uZ^2
                  ((1+Z)^2 + 1)^2

           4*(1+Z)
     = --------------- * uZ^2
       ((1+Z)^2 + 1)^2
```

See fortran subroutine `save_vel()` for more information.

### Model File (e.g. `star.08.99.mod`)

#### Header
None.

#### Body
| **Column** | **Description** | **Example** | **Variable** |
| :---:  | :--- | :--- | :--- |
| **0**  | Observation index (starting with 001) | `001` | `n` |
| **1**  | Order index (starting with 01) | `08` | `m` |
| **2**  | Pixel index (0001-4021 for HIRES) | `4021` | `i` |
| **3**  | Observed spectrum [arbitrary units] | `47425.8` | `pix_imn(i,m,n)` |
| **4**  | Unnormalized model spectrum | `73123.5` | `mod_imn(i,m,n)` |
| **5**  | Wavelengths [A] in the observatory (iodine) frame | `5522.335491` | `wav_imn(i,m,n)` |
| **6**  | Wavelengths [A] in the stellar frame | `5522.889612` | `wav_imn(i,m,n)/(1+zn(n))` |
| **7**  | Normalization vector such that mod*nrm matches obs | `0.642063` | `nrm_imn(i,m,n)` |
| **8**  | Smoothed mod*nrm that can be used to psuedo-normalize obs | `47803.4688` | `bar_imn(i,m,n)` |
| **9**  | Flag indicating pixel was rejected as an outlier | `1` | `rej_imn(i,m,n)` |
| **10** | Flag indicating pixel was rejected due to telluric lines | `0` | `tel_imn(i,m,n)` |
| **11** | Flag indicating pixel was rejected due to meteor | `0` | `met_imn(i,m,n)` |

Symbols:
- zn = velocity divided by the speed of light

This file has one row per observation, per order, per pixel.
Use read_grand_mod.pro to read the contents of a .mod file into IDL.
Observatory frame features include the iodine cell, telluric lines, and to some extent pixels.
Stellar frame features consist mainly of the deconvolved stellar spectrum.
Column 4 (MOD) times column 7 (NRM) is the normalized fit of the observation.

#### Sample usage in IDL

Read a .mod file into IDL:
```
IDL> read_grand_mod, 'sunsim09.08.99.mod', model
```

Examine the resulting IDL structure:
```
IDL> help, model
MODEL           STRUCT    = -> <Anonymous> Array[4021, 100]
IDL> help, /struct, model
** Structure <35004e8>, 9 tags, length=40, data length=35, refs=1:
   WOB             DOUBLE           5522.3355
   WST             DOUBLE           5522.8896
   OBS             FLOAT           47425.8
   MDL             FLOAT           73123.5
   NRM             FLOAT          0.642063
   BAR             FLOAT           47803.5
   REJ             BYTE         1
   TEL             BYTE         0
   MET             BYTE         0
```

Plot an observed spectrum and overplot the corresponding model fit:
```
IDL> m=model[*,0] & plot, m.wob, m.obs, xsty=3 & oplot, m.wob, m.mdl*m.nrm, co=255
```

Plot the fit residual for an observed spectrum. Note any divergence at the long wavelength end.
```
IDL> m=model[*,0] & plot, m.wob, m.obs-m.mdl*m.nrm, xsty=3  
```

Plot fit residuals in the observatory frame for all observations. Crop outliers at ends of the order.
```
IDL> m=model[*,0] & plot, m.wob, m.obs-m.mdl*m.nrm, ps=3, xsty=3, yr=300*[-1,1]        
IDL> for i=1,99 do begin & m=model[*,i] & oplot, m.wob, m.obs-m.mdl*m.nrm, ps=3 & endfor
```

Plot fit residuals in the stellar frame for all observations. Crop outliers at ends of the order.
```
IDL> m=model[*,0] & plot, m.wst, m.obs-m.mdl*m.nrm, ps=3, xsty=3, yr=300*[-1,1]         
IDL> for i=1,99 do begin & m=model[*,i] & oplot, m.wst, m.obs-m.mdl*m.nrm, ps=3 & endfor
```

### Stellar Template File (e.g. `star.08.99.tem`)

#### Header
None.

#### Body
| **Column** | **Description** | **Example** | **Variable** |
| :---: | :--- | :--- | :--- |
| **0** | Template node index (00001-30000) | `30000` | `s` |
| **1** | Wavelength [Angstroms] | `5986.05036` | `ww_s(s)` |
| **2** | Template spectrum at end of current internal `grand` iteration | `0.852927` | `tems(s)` |
| **3** | Template spectrum at end of previous internal `grand` iteration | `0.862997` | `tem0(s)` |
| **4** | Solar spectrum (Kurucz, Furenlid, Brault, & Testerman 1984) | `0.985750` | `sunw(ww_s(s))` |

#### Notes
At the start of each internal iteration after the first, `grand.F` reads and smoothes the template spectrum from the previous iteration. See the subroutine `smoo_tem( )` for the algorithm.

### Line spread functions (e.g., `*.lsf`)

#### Header

None.

#### Body
| **Column** | **Description** | **Example** | **Variable** |
| :---: | :--- | :--- | :--- |
| **0** | Observation index (start with 001, not 000) | `001` | `n` |
| **1** | Order index | `08` | `m` |
| **2** | LSF index for spatial variations along an echelle order (1-9)  | `9` | `q` |
| **3-15** | LSF values at 13 node | `-0.00250  0.00020  0.00242  0.01174  0.06333  0.23470  0.36176  0.25491  0.06085  0.00881  0.00301  0.00090 -0.00032` | `lsfq` |

#### Notes

Need to document the pixel offsets of the 13 nodes in the LSF used in `grand` and `grlsf`.

The number of rows in a `*.lsf` file is (number of observations) x (number of orders) x (number of LSFs along an echelle order). Since `grand` operates on one order at a time, the number of orders is 1. The number of LSFs along an echelle order is 9. Thus, the number of rows is 9 times the number of observations in the input `obslist` file.

### Other files

- .nim files:
  - col0 = 
  - col1 = 
  - etc.
- .vln files:
  - col0 = 
  - col1 = 
  - etc.
- .vsl files:
  - col0 = 
  - col1 = 
 - etc.
- .mod files:
 - col0 = 
  - col1 = 
  - etc.
- .tem files:
  - col0 = 
  - col1 = 
  - etc.
- .lines_find files:
  - col0 = 
  - col1 = 
  - etc.
- .lines_id files:
  - col0 = 
  - col1 = 
  - etc.

# Fortran code that calculates grand LSF on subpixel grid

# Building `grlsf`

Run `make grlsf` in the `grandsol/fort` source directory.

Executable is written to the `grandsol/bin` directory.

## Calling syntax for `grlsf`
```
usage: grlsf lsffile obs m k

Calculate line spread function with NJ=201 points from -10.0 to +10.0 pixels
in steps of 0.1 pixel, based on NQ=13 LSF node values in the specified *.lsf
file. The node interpolation algorithm is a copy of what is used in grand. 

positional arguments:
  lsffile               *.lsf file produced by grand
  obs                   observation index from obslist
  m                     order index (starting with 1, not 0)
  k                     node index across echelle order (1-9)

output (written to standard outout):
  line 1                NQ NJ obs m q
  line 2                LSF value at NQ=13 nodes
  line 3                LSF value at NJ=201 interpolated points

notes:
  Interpolation logic assumes NQ=13 and NJ=201.
  m=1 is HIRES order 71 (4977-5065 Angstroms)
  m=16 is HIRES order 56 (6310-6368 Angstroms)
```

## Input file

See description of the line spread function file (`*.lsf`) that `grand` outputs.

## Output (written to standard output)

#### Header
 **Row** | **Description** | **Example** | **Variable** |
| :---: | :--- | :--- | :--- |
| **1** | `NQ  NJ  obs  m i q` | `13 201   3  8 9` | `_QDIM_  _JDIM_  obs  m  q` |

Symbols:
- NQ = Number of nodes per LSF (interpolation logic assumed `NQ=13`)
- NJ = Number of finely sampled points in output LSF (interpolation logic assumed `NJ=201`)
- obs = observation index specified in argument list
- m = order index specified in argument list
- k = node index specified in argument list

#### Body
| **Row** | **Description** | **Example** | **Variable** |
| :---: | :--- | :--- | :--- |
| **2** | LSF value at nodes | `0.00605  0.01263  0.00186  0.00501  0.06210  0.23414  0.35411  0.24396  0.05303  0.00591  0.01016  0.00128 -0.00576` | `lsfq` |
| **3** | LSF value at interpolated points | `0.00000  0.00000  0.00001  0.00001  0.00002  0.00004  0.00005  0.00007  0.00010  0.00012  0.00015  0.00018  0.00021  0.00024  0.00027  0.00030  0.00033  0.00036  0.00039  0.00042  0.00045  0.00048  0.00052  0.00055  0.00058  0.00061  0.00065  0.00068  0.00072  0.00076  0.00079  0.00083  0.00087  0.00091  0.00094  0.00098  0.00102  0.00106  0.00109  0.00112  0.00115  0.00116  0.00116  0.00116  0.00114  0.00111  0.00107  0.00101  0.00095  0.00087  0.00079  0.00071  0.00063  0.00055  0.00047  0.00039  0.00031  0.00024  0.00018  0.00012  0.00008  0.00005  0.00003  0.00003  0.00004  0.00008  0.00014  0.00023  0.00035  0.00050  0.00069  0.00092  0.00121  0.00157  0.00199  0.00252  0.00314  0.00389  0.00476  0.00576  0.00690  0.00818  0.00958  0.01109  0.01270  0.01439  0.01614  0.01791  0.01970  0.02148  0.02323  0.02493  0.02656  0.02810  0.02952  0.03080  0.03191  0.03282  0.03351  0.03396  0.03416  0.03411  0.03380  0.03325  0.03247  0.03146  0.03028  0.02891  0.02739  0.02575  0.02399  0.02215  0.02024  0.01830  0.01636  0.01443  0.01255  0.01077  0.00910  0.00758  0.00622  0.00503  0.00401  0.00317  0.00248  0.00195  0.00153  0.00122  0.00099  0.00083  0.00073  0.00066  0.00063  0.00062  0.00063  0.00066  0.00069  0.00074  0.00079  0.00083  0.00088  0.00091  0.00094  0.00095  0.00095  0.00094  0.00091  0.00087  0.00082  0.00076  0.00070  0.00063  0.00055  0.00048  0.00041  0.00033  0.00027  0.00020  0.00014  0.00009  0.00004 -0.00001 -0.00005 -0.00009 -0.00013 -0.00017 -0.00021 -0.00025 -0.00029 -0.00033 -0.00037 -0.00041 -0.00044 -0.00046 -0.00048 -0.00049 -0.00049 -0.00049 -0.00048 -0.00046 -0.00043 -0.00040 -0.00037 -0.00035 -0.00032 -0.00029 -0.00026 -0.00023 -0.00020 -0.00017 -0.00014 -0.00012 -0.00009 -0.00007 -0.00005 -0.00004 -0.00002 -0.00001 -0.00001 -0.00000  0.00000` | `lsfj` |
