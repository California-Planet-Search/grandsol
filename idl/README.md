## IDL Procedures for the Grand Solution

### Meteor file

##### Create a meteor file

* Set ```GRAND_RAWDIR``` environment variable to point at the directory containing raw fits files.
* Prepare an ```obslist``` file that specifies the set of observations to analyze with iGrand.
* Run the IDL procedure ```meteor_obslist, 'obslist'``` to create a ```meteor``` file.
* Add a ```METEOR meteor``` line to the header of the ```obslist``` file. 

##### Dependencies

The IDL procedure ```meteor_obslist``` requires several other IDL procedures: ```gaussbroad```, ```m_fndpks```, ```m_fords```, ```m_getarc```, ```meteor```, and ```read_grand_obslist```.

##### Format of the meteor file

| **Column** | **Description** | **Example** | **Variable** |
| :---: | :--- | :--- | :--- |
| **0** | Observation ID | `rj100.146` | `obsid` |
| **1** | row where order 8 crosses column 2000 | `289.71` | `y7at2k` |
| **2** | column where meteor crosses order 1 | `3032` | `meteor[0]` |
| **3** | column where meteor crosses order 2 | `2936` | `meteor[1]` |
| **4** | column where meteor crosses order 3 | `2837` | `meteor[2]` |
| **5** | column where meteor crosses order 4 | `2737` | `meteor[3]` |
| **6** | column where meteor crosses order 5 | `2635` | `meteor[4]` |
| **7** | column where meteor crosses order 6 | `2531` | `meteor[5]` |
| **8** | column where meteor crosses order 7 | `2424` | `meteor[6]` |
| **9** | column where meteor crosses order 8 | `2315` | `meteor[7]` |
| **10** | column where meteor crosses order 9 | `2204` | `meteor[8]` |
| **11** | column where meteor crosses order 10 | `2089` | `meteor[9]` |
| **12** | column where meteor crosses order 11 | `1972` | `meteor[10]` |
| **13** | column where meteor crosses order 12 | `1853` | `meteor[11]` |
| **14** | column where meteor crosses order 13 | `1730` | `meteor[12]` |
| **15** | column where meteor crosses order 14 | `1604` | `meteor[13]` |
| **16** | column where meteor crosses order 15 | `1474` | `meteor[14]` |
| **17** | column where meteor crosses order 16 | `1341` | `meteor[15]` |
| **18** | median amplitude of measured meteor features [DN] | `72.1` | `zmed` |
| **19** | median of absolute value of residual of meteor location fit [pixels] | `0.403` | `mar` |

### IDL procedure calling syntax and docstrings
##### fmtnum.pro
```
function fmtnum, data, fmt
;Return a string containing DATA formatted according to FMT with the
;following additional rules:
; 1) strip all leading and trailing spaces
; 2) strip all trailing zeros after the decimal point
; 3) strip trailing decimal point
```

##### locate_interval.pro
```
pro locate_interval, x, xbeg, xend, ibeg, iend, nadd=nadd, stop=stop
;
;Given a monotonic input vector and two limiting values, return two indexes
;that bracket all data values between the limiting values. Extra data points
;may be included (or discarded from) both ends of the interval.
;
;Inputs:
; x (vector) monotonically increasing or decreasing data values
; xbeg (scalar) one limiting value that defines the data interval
; xend (scalar) the other limiting value that defines the data interval
; [nadd=] (scalar) number of extra data points to include at both ends of
;   the interval (or to discard, if nadd is negative).
; [/stop] (switch) halt execution requested interval is null.
;
;Outputs:
; ibeg (scalar) lower index of data values in the requested interval.
;   return value is -1 if requested interval is null.
; iend (scalar) upper index of data values in the requested interval
;   return value is -1 if requested interval is null.
;
;Notes:
; For large data vectors, this routine is much faster than using where().
; Passing a named variable (x) is faster passing a structure element (d.x).
;
;History:
; 2008-Sep-15 Valenti  Initial coding.
```

##### meteor_obslist.pro
```
pro meteor_obslist, oblist_file, print=print, plot=plot
;Create a 'meteor' file in current working directory.
;
;Inputs:
; GRAND_RAWDIR (environment variable) directory containing raw fits files
; obslist_file (string) name of an obslist file containing a list of obsids
; [/print] (switch) print results to terminal as well as the meteor file
; [/plot] (switch) create cool diagnostic plot in terminal window
;
;Output:
; 'meteor' (disk file in current working directory) location of meteor
;   feature for each order in each observation in 'obslist' file.
```

##### rdnso.pro
```
pro rdnso, w, s, wmin, wmax
;Reads requested portion of National Solar Observatory Atlas #1, which is an
;  optical, full-disk, FTS spectrum of a relatively inactive Sun.
; w (output vector) wavelength scale for s
; s (output vector) reesidual intensity segment of atlas.
; wmin (input scalar) lowest wavelength to return in w
; wmax (input scalar) highest wavelength to return in w
;Uses the data file nso.bin which has format:
;  intarr(1073):  header containing "base wavelengths" (in Angstroms) for each
;		   atlas segment, i.e. the integer part of the first wavelength
;		   in each segment.
;  fltarr(1024):  wavelengths increments for points in atlas segment 0.  Add
;		   the base wavelength for this segment to recover the actual
;		   wavelength points.
;  intarr(1024):  disk integrated, solar, residual intensity, multiplied by
;		   a scale factor of 30000 to match range of signed integers.
;  Plus an additional {fltarr(1024), intarr(1024)} pairs for the remaining
;   1072 atlas segments.  The total file length is 6,594,658 bytes, which is
;   2*1073 + 1073 * (4*1024 + 2*1024).
;03-Aug-90 JAV	Create.
;13-Oct-90 JAV	NSO atlas spectrum returned as double precision.
;12-Aug-90 JAV	Cleaned up Paul Butler's IDL adaptation of the ANA version.
;19-Oct-94 JAV	New data file location for use on casa machines.
;10-Apr-95 JAV	Reformatted and commented to improve readability of code.
;		Extra pixels no longer added outside requested range.
;01-May-99 JAV	Made initial value of for loop 0L instead of 0 to allow for
;                more than 32K segments in a single read.
;08-Sep-22 JAV  Added /swap_if_little_endian to support Intel Macs
;2016-Jan-25 Valenti Directory containing nso.bin now specified by the
;              environment variable $GRAND_IREFDIR.
```


##### read_iodine.pro
```
pro read_iodine, cell, lab, temp, wbeg, wend, wvac, tran, utran $
               , nadd=nadd, silent=silent, noap=noap
;
;Interpolate/extrapolate iodine transmission spectra for the requested cell
;to the requested temperature.
;
;Inputs:
; cell (string) name of the iodine cell ('lick', 'aat', 'keck', 'ctio')
; lab (string) lab where iodine spectrum was obtained ('nso', 'pnnl', 'nist')
; temp (scalar) temperature (degrees C) for returned transmission spectrum
; wbeg (scalar) minimum vacuum wavelength to return
; wend (scalar) maximum vacuum wavelength to return
; [nadd=] (scalar) number of extra data points to include at both ends of
;   the wavelength interval (or to discard, if nadd is negative).
; [/noap] (switch) read the unapodized version of the spectrum
;
;Outputs:
; wvac (vector) vacuum wavelength scale for iodine cell transmission spectrum
; tran (vector) transmission spectrum for specified iodine cell
; utran (vector) uncertainty in transmission spectrum (if available)
;
;History:
; 2008-Sep-16 Valenti  Initial coding.
; 2009-Jul-21 Valenti  Added lab and temp arguments. Changed data files.
```


##### simobs_one.pro
```
pro simobs_one, star, iod, elsf, ech, velsum, ston, file $
              , verbose=verbose, seed=seed
;
;Inputs:
; star (string) - identifier for an intrinsic stellar spectrum, e.g. 'sun'
; iod (string) - identifier for an iodine cell spectrum, e.g. 'keck-nso-52.5'
; elsf (string) - identifier for an effective LSF, e.g. 'jay2011', 'laser_y'
; ech (string) - identifier for echelle spectrograph, e.g. 'hires'
; velsum (scalar) - [m/s] relativistically summed velocities
; ston (scalar) - continuum signal to noise [0 --> noiseless]
; file (string) - name of output file
; [seed=] - seed for random number generator (gives repeatable noise)
;
;Output:
; Output file is created in the current working directory.
;
;History:
; 2012-Oct-23 Valenti Adapted from make_output.pro
```


##### simobs_set.pro
```
pro simobs_set, star, iod, elsf, ech, bcmax, rvk, vsys, ston, nobs $
              , seed=seed
;
;Inputs:
; star (string) - identifier for an intrinsic stellar spectrum, e.g. 'sun'
; iod (string) - identifier for an iodine cell spectrum, e.g. 'keck-nso-52.5'
; elsf (string) - identifier for an effective LSF, e.g. 'jay2011', 'laser_y'
; ech (string) - identifier for echelle spectrograph, e.g. 'hires'
; bcmax (scalar) - [m/s] maximum barycentric correction
; rvk (scalar) - [m/s] semi-amplitude of stellar velocity perturbations
; vsys (scalar) - [m/s] radial velocity of the system
; ston (scalar) - cotinuum signal to noise ratio [value of 0 mean no noise]
; nobs (integer) - number of simulated observations to create
; [seed=] - seed for random number generator (gives repeatabable noise)
;
;Output:
; Data files are created in the current working directory.
;   - Simulated observations named <star>.###
;   - Observation list file named obslist
;
;History:
; 2012-Oct-23 Valenti Adapted from make_set.pro
```


##### wdsk_simobs.pro
```
Pro WDsk_simobs,Var,InFile,Arg3,Arg4,New=New,Old=Old,Insert=Insert,swap=swap
;General purpose routine to store variables to disk.
; Var (input variable, any type/size) variable to be stored.
; InFile (input string) name of file in which to store variable.
; Arg3 and Arg4 are Record and/or Comment in any order...
;   Record (optional input scalar) file index at which to store variable.
;     1=first storage position, 0=last storage position.
;   Comment (optional input string) comment to store with variable.
; /New (logical) forces the creation of a new file. Record must be 0 or 1.
; /Old (logical) updates existing file. Error if file doesn't exist.
; /Append (logical) allows insertion of data in the middle of file.
;   New variable must have the same length as the previous one.
;12-Jul-91 JAV	Create.
;11-May-94 JAV  Added on_error trap. 
;10-Sep-94 JAV	Added path expansion.
;20-Aug-10 JAV	Added /swap keyword to write with non-native endian-ness
```
