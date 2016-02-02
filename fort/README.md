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
| **0** | Observation index | `001` | `n` |
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

