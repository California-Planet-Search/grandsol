# Fortran code for the Grand Solution Doppler code

## Environment variables

Set GRAND_REFDIR to point at the directory containing reference files used
by 'grand'. The usual reference file directory is 'grandsol/gref'.
See 'grandsol/gref/README' for more information.


## Building 'grand'

Run 'make' in the 'fort' source directory.


## Calling syntax
```
grand obslist label m [out=] [mod1] [211111] [fudge-] [vorb+] [nitf=5] ...
  obslist      : file specifying observations
  label        : root to use for output filenames
  m            : echelle order to analyze
  out=<STDOUT> : log file for diagnostics
  vel2         : velocity output      | 
  tem1         : template output      | 0=no output
  mod0         : model output         | 1=most recent interation
  lsf1         : lsf output           | 2=every iteration
  nrm0         : normalization output | 
  lin0         : line output          | 
  210100       : vel-tem-mod-lsf-nrm-lin output flags
  fudge+       : apply B star correction, fudge- to suppress
  vorb+        : include orbital velocities in obslist, 
  nito=0       : starting iteration number
  nitf=10      : final iteration number
grand obslist sunsim 8 out=sunsim.08.log 211111 fudge-
grand obslist2 sunsim2 8 out=sunsim2.08.log 211111 fudge- vorb+
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

Plot the fit residual an observed spectrum. Note any divergence at the ling wavelength end.
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
- .lsf files:
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

