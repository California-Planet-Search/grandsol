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

;Check syntax.
  if n_params() lt 7 then begin
    print, 'syntax: read_iodine, cell, lab, temp, wbeg, wend, wvac, tran' $
         + ' [,utran ,nadd=]'
    print, "  e.g.: read_iodine, 'keck', 'pnnl', 67.3, 5500, 5510, w, t" $
         + ", nadd=1"
    print, "  e.g.: read_iodine, 'keck', 'nso', 53.1, 3000, 9000, w, t" $
         + ", nadd=1"
    return
  endif

;Define static local variables (using a common block).
  common common_iodine, iodine_cell, iodine_lab, iodine_noap $
        , iodine_temps, iodine_wvac, iodine_tran, iodine_utran

;Default values for optional parameters.
  if n_elements(nadd) eq 0 then nadd = 0
  if n_elements(noap) eq 0 then noap = 0

;Directory where iodine spectrum is saved.
  dir = getenv('GRAND_IREFDIR') + path_sep() + 'iod' + path_sep()
  if keyword_set(noap) then dir += 'noap/'

;Restore entire iodine spectrum, unless it is already in common block.
  if n_elements(iodine_cell) eq 0 then iodine_cell = ''
  if n_elements(iodine_lab) eq 0 then iodine_lab = ''
  if n_elements(iodine_noap) eq 0 then iodine_noap = 0
  if iodine_cell ne cell or iodine_lab ne lab $
                         or iodine_noap ne noap then begin
    file = dir + 'iodine_' + cell + '_' + lab + '.sav'
    if ~keyword_set(silent) then begin
      print, 'restoring ' + file
    endif
    restore, file
    iodine_cell = cell
    iodine_lab = lab
    iodine_noap = noap
  endif

;Determine indexes of requested wavelength segment.
  locate_interval, iodine_wvac, wbeg, wend, ibeg, iend, nadd=nadd
  if iend lt ibeg then begin
    n = n_elements(iodine_wvac)
    diag, -1, 'request [' $
            + strtrim(string(wbeg, form='(f19.3)'), 2) $
            + ',' $
            + strtrim(string(wend, form='(f19.3)'), 2) $
            + '] not in allowed range [' $
            + strtrim(string(iodine_wvac[0], form='(f19.3)'), 2) $
            + ',' $
            + strtrim(string(iodine_wvac[n-1], form='(f19.3)'), 2) $
            + ']', /stop
  endif

;Handle case when transmission is available for only 1 temperature.
  ntemp = n_elements(iodine_temps)
  if ntemp eq 1 then begin
    if temp eq iodine_temps then begin
      wvac = iodine_wvac[ibeg:iend]
      tran = iodine_tran[ibeg:iend]
      if arg_present(utran) and keyword_set(iodine_utran) then begin
        utran = iodine_utran[ibeg:iend]
      endif
      return
    endif else begin
      diag, -1, 'only one spectrum at ' $
              + strtrim(fix(iodine_temps), 2) $
              + " C, so can't interpolate to " $
              + strtrim(temp, 2) + ' C', /stop
    endelse
  endif

;Prepare to interpolate in temperature.
  jtemp = 0 > value_locate(iodine_temps, temp) < (ntemp - 2)

;Extract relevant vacuum wavelengths (identical for all temperatures).
  wvac = iodine_wvac[ibeg:iend]

;Interpolate/extrapolate transmission spectra to requested temperature.
  fhi = (temp - iodine_temps[jtemp]) $
        / (iodine_temps[jtemp+1] - iodine_temps[jtemp])
  flo = 1 - fhi
  tran = flo * iodine_tran[ibeg:iend,jtemp] $
       + fhi * iodine_tran[ibeg:iend,jtemp+1]

;Calculate uncertainty in interpolated spectrum, if requested.
  if n_params() ge 6 and keyword_set(iodine_utran) then begin
    utran = sqrt( flo^2 * iodine_utran[ibeg:iend,jtemp]^2 $
                + fhi^2 * iodine_utran[ibeg:iend,jtemp+1]^2 )
  endif

end
