pro read_grand_obslist, obli, obslist
;
;Purpose:
; Read an observation list file in the current working directory.
;
;Input:
; [obslist] (string) - name of observation list file in current working
;   directory. Default name is 'obslist'.
;
;Output:
; obli (structure) - contents read from the observation list file
;  .vsys (scalar) - systemic velocity from 'VSYST =' header line
;  .obsdir (string) - observation directory from 'RJDIR =' header line
;  .nobs (scalar) - number of observations described in obslist
;  .obs (structure vector[nobs]) individual observation data
;    .file (string) - name of observation file
;    .barycorr (double) - [m/s] barycentric correction
;    .vorb (double) - [m/s] orbital velocity of companion

;Syntax
  if n_params() lt 1 then begin
    print, 'syntax: read_grand_obslist, obli [,obslist]'
    return
  endif

;Default value for optional parameters
  if ~keyword_set(obslist) then obslist = 'obslist'

;Clear return variable.
  obli = ''

;Check that observation list file exists in current working directory.
  if ~file_test(obslist) then begin
    print, 'file not found: ' + obslist
    return
  endif

;Read entire observation list file
  nline = file_lines(obslist)
  lines = strarr(nline)
  openr, unit, obslist, /get_lun
  readf, unit, lines
  free_lun, unit

;Find and parse 'VSYST' specification.
  vsyst = 0.0
  iline = where(strmid(lines, 0, 5) eq 'VSYST', count)
  if count gt 2 then message, 'multiple VSYST records in ' + obslist
  if count eq 1 then begin
    line = lines[iline[0]]
    ipos = strpos(line, '=')
    vsyst = float(strmid(line, ipos+1))
  endif

;Find and parse 'RJDIR' or 'OBSDIR' specification.
  obsdir = ''
  iline = where(strmid(lines, 0, 5) eq 'RJDIR' $
             or strmid(lines, 0, 6) eq 'OBSDIR', count)
  if count gt 2 then message, 'multiple RJDIR or OBSDIR records in ' + obslist
  if count eq 1 then begin
    line = lines[iline[0]]
    ipos1 = strpos(line, '"')
    ipos2 = strpos(line, '"', ipos1+1)
    obsdir = strmid(line, ipos1+1, ipos2-ipos1-1)
    if ~file_test(obsdir) then begin
      print, obslist + ' sets RJDIR|OBSDIR to nonexistent directory:'
      print, '  ' + obsdir
      return
    endif
  endif

;Read observation list, assumed to start immediately after RJDIR line
  iobs0 = iline[0] + 1
  words = strsplit(lines[iobs0], /extract, count=ncol)
  nobs = nline - iobs0
  data = strarr(ncol, nobs)
  for iobs=0, nobs-1 do begin
    line = lines[iobs0 + iobs]
    words = strsplit(line, /extract, count=nword)
    if nword ne ncol then message, 'inconsistent number of columns'
    data[*,iobs] = words
  endfor

;Parse data from observation list
  obsfile = reform(data[1,*])
  barycorr = double(reform(data[3,*]))
  vorb = double(reform(data[4,*]))

;Build return structure
  obs = replicate({file:'', barycorr:0d0, vorb:0d0}, nobs)
  obs.file = obsfile
  obs.barycorr = barycorr
  obs.vorb = vorb
  obli = $
    { obslist : obslist $
    , vsys    : vsyst   $
    , obsdir  : obsdir  $
    , obs     : obs     $
    }

end
