## Fortran code for the Grand Solution Doppler code

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


### input files
#### `obslist`
##### Header Description

| **Example** | **Default** | **Description** |
| :--- | :--- | :--- |
| `VSYST = -16400 m/s` | *Mandatory* | Radial velocity of center of mass |
| `RJDIR = "\Users\valenti\obs\"` | `".\"` | Directory containing observations to fit |
| `METEOR \Users\valenti\ref\meteor.dat` | `NONE` | File containing meteor location |

##### Body Description

- col0 = observation index
- col1 = file name?
- col2 = barycentric correction? (m/s)
- col3 = initial velocity guess? (m/s)

See fortran subroutine `read_raw`


### output file format
- .vel files:
  - col0 = 
  - col1 = 
  - col2 = 
  - col3 = 
  - col4 =
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

