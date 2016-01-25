pro write_grand_obslist, obli
;
;Purpose:
; Write an observation list file in the current working directory.
;
;Input:
; obli (structure) - contents read from the observation list file
;  .vsys (scalar) - systemic velocity from 'VSYST =' header line
;  .obsdir (string) - observation directory from 'RJDIR =' header line
;  .obs (structure vector[nobs]) individual observation data
;    .file (string) - name of observation file
;    .barycorr (double) - [m/s] barycentric correction
;    .vorb (float) - [m/s] orbital velocity of companion

;Syntax
  if n_params() lt 1 then begin
    print, 'syntax: write_grand_obslist, obli'
    return
  endif

;Read entire observation list file
  openw, unit, obli.obslist, /get_lun

;Write header
  printf, unit, 'VSYST = ' + fmtnum(obli.vsys, '(f9.2)') + ' m/s'
  printf, unit, 'RJDIR = "' + obli.obsdir + '"'

;Write observation data.
  nobs = n_elements(obli.obs)
  for iobs=0, nobs-1 do begin
    oi = obli.obs[iobs]
    printf, unit, string(iobs+1, form='(i3.3)') + '  ' $
                , oi.file + '  0' $
                , string(oi.barycorr, form='(f14.5)') $
                , string(oi.vorb, form='(f14.5)')
  endfor

;Close file
  free_lun, unit

end
