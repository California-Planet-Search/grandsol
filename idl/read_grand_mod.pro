pro read_grand_mod, file, fit, upix, uobs, print=print
;
;Purpose:
; Read fitted spectrum data from a .mod file produced by grand.
;
;Input:
; file (string) name of .mod file containing model fit data
; [/print] (switch) print diagnostics to screen
;
;Output:
; fit (structure[npix,nobs]) model spectrum data from a .mod file
;  .wob - wavelengths in the observatory (iodine) frame [Angstroms]
;  .wst - wavelengths in the stellar frame [Angstroms]
;  .obs - observed spectrum
;  .mdl - model spectrum (star*iodine x elsf) without normalization
;  .nrm - normalization such that mdl*nrm matches obs
;  .bar - smoothed mdl*nrm that can be used to psuedo-normalize obs
;  .rej - flag indicating pixels with dynamically rejected outliers
;  .tel - flag indicating pixels affected by telluric lines
;  .met - flag indicating pixels affected by the 'meteor' artifact
; upix (integer[npix]) one-based pixel index for fit structure
; uobs (integer[nobs]) id of observations in the fit structure

;Syntax
  if n_params() lt 2 then begin
    print, 'syntax: read_grand_mod, file, fit [upix, uobs,/print]'
    print, "  e.g.: read_grand_mod, 'sunsim.08.99.mod', fit"
    return
  endif

;Test that file exists.
  if ~file_test(file) then begin
    message, /info, 'file not found: ' + file
    return
  endif

;Determine number of rows in file.
  nline = file_lines(file)

;Allocate array
  data = dblarr(12, nline)

;Read model data that grand wrote to disk.
  openr, unit, file, /get_lun
  readf, unit, data
  free_lun, unit

;Determine unique observation, order, and pixel identifiers.
  for i=0, 2 do begin
    d = fix(reform(data[i,*]))
    u = d[uniq(d, sort(d))]
    case i of
      0: uobs = u
      1: uord = u
      2: upix = u
    endcase
  endfor
  if n_elements(uord) ne 1 then message, 'no logic to handle multiple orders'
  nuobs = n_elements(uobs)
  nupix = n_elements(upix)

;Print diagnostics
  if keyword_set(print) then begin
    print, file + ' (order ' $
         + strtrim(uord, 2) + '), ' $
         + strtrim(nuobs, 2) + ' observations, ' $
         + strtrim(nupix, 2) + ' pixels/order'
  endif

;Construct output structure.
  rec = $
    { wob : 0d0 $
    , wst : 0d0 $
    , obs : 0.0 $
    , mdl : 0.0 $
    , nrm : 0.0 $
    , bar : 0.0 $
    , rej : 0b  $
    , tel : 0b  $
    , met : 0b  $
    }
  fit = replicate(rec, nupix, nuobs)

;Populate output structure.
  ipix = value_locate(upix, round(reform(data[2,*])))
  iobs = value_locate(uobs, round(reform(data[0,*])))
  fit[ipix,iobs].wob = reform(data[5,*])
  fit[ipix,iobs].wst = reform(data[6,*])
  fit[ipix,iobs].obs = reform(data[3,*])
  fit[ipix,iobs].mdl = reform(data[4,*])
  fit[ipix,iobs].nrm = reform(data[7,*])
  fit[ipix,iobs].bar = reform(data[8,*])
  fit[ipix,iobs].rej = reform(data[9,*])
  fit[ipix,iobs].tel = reform(data[10,*])
  fit[ipix,iobs].met = reform(data[11,*])

end
