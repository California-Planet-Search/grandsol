pro read_grand_mods, mods, modinfo, upix, uobs, mlist=mlist, print=print
;
;Purpose:
; Read model data from grand for requested (default all) echelle orders.
;
;Inputs:
; Output model files (.mod) in the current directory
; [mlist=] (vector[nm]) subset of orders to read, starting with 1
; [/print] (switch) print diagnostics to screen
;
;
;Outputs:
; mods (structure[npix,nm,nobs]) model data from .mod files
;  .wob - wavelengths in the observatory (iodine) frame [Angstroms]
;  .wst - wavelengths in the stellar frame [Angstroms]
;  .obs - observed spectrum
;  .mdl - model spectrum (star*iodine x elsf) without normalization
;  .nrm - normalization such that mdl*nrm matches obs
;  .bar - smoothed mdl*nrm that can be used to psuedo-normalize obs
;  .rej - flag indicating pixels with dynamically rejected outliers
;  .tel - flag indicating pixels affected by telluric lines
;  .met - flag indicating pixels affected by the 'meteor' artifact
; modinfo (structure[nm]) auxiliary information for each order
;  .m - echelle order, starting at 1 for first order in an observation
;  .wb - minimum observatory wavelength for each order and all observations
;  .we - maximum observatory wavelength for each order and all observations
;  .file - name of file
;  .iter - last iteration number
; upix (integer[npix]) one-based pixel index for fit structure
; uobs (integer[nobs]) id of observations in the fit structure

;syntax
  if n_params() lt 2 then begin
    print, 'syntax: read_grand_mods, mods, modinfo [,upix,uobs,mlist=,/print]'
    return
  endif

;internal parameters
  sfx = 'mod'

;get information about grand output files in current directory.
  grand_files, gf
  if ~keyword_set(gf[sfx].uord) then begin
    message, /info, "no '" + sfx $
           + "' files produced by grand in current directory"
    return
  endif
  uord = gf[sfx].uord

;create complete order list, if user did not specify an order list
  if ~keyword_set(mlist) then mlist = fix(gf['mod'].uord)
  nm = n_elements(mlist)

;Allocate output info structure
  rec =          $
    { m    : 0   $
    , wb   : 0d0 $
    , we   : 0d0 $
    , file : ''  $
    , iter : ''  $
    }
  modinfo = replicate(rec, nm)

;read model data for each requested order
  for im=0, nm-1 do begin
    m = mlist[im]
    iord = where(fix(gf[sfx].uord) eq m, count)
    if count eq 0 then begin
      message, 'no model (.mod) for order ' + strtrim(m, 2)
    endif
    iter = gf[sfx].ilast[iord]
    file_mod = gf[sfx].flast[iord]
    read_grand_mod, file_mod, fit, upixm, uobsm
    if im eq 0 then begin
      upix = upixm
      uobs = uobsm
      npix = n_elements(upix)
      nobs = n_elements(uobs)
      mods = replicate(fit[0], npix, nm, nobs)
    endif else begin
      if max(upixm ne upix) eq 1 then begin
        message, 'orders have different pixels'
      endif
      if max(uobsm ne uobs) eq 1 then begin
        message, 'orders have different observations'
      endif
    endelse
    mods[*,im,*] = fit
    modinfo[im].iter = iter
    modinfo[im].file = file_mod
  endfor
  npix = n_elements(upix)

;determine first and last valid wavelength for each requested order
  for im=0, nm-1 do begin
    m = mlist[im]

;load output info structure
    modinfo[im].m  = m
    modinfo[im].wb = min(mods[0,im,*].wob)
    modinfo[im].we = max(mods[npix-1,im,*].wob)

;print diagnostics
    if keyword_set(print) then begin
      if im eq 0 then print, '        ___WAVELENGTH___'
      if im eq 0 then print, 'order     begin      end'
      print, form='(3x,i2.2,1x,2f9.2)' $
           , m, modinfo[im].wb, modinfo[im].we
    endif
  endfor

end
