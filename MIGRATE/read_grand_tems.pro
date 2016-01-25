pro read_grand_tems, tems, teminfo, mlist=mlist, print=print
;
;Purpose:
; Read template stellar spectra for requested (default all) echelle orders.
;
;Inputs:
; Output template files (.tem) in the current directory
; [mlist=] (vector[nm]) subset of orders to read, starting with 1
; [/print] (switch) print diagnostics to screen
;
;
;Outputs:
; tems (structure[ntem,nm]) template stellar spectrum from .tem files
;  .w - wavelengths for template spectrum [Angstroms]
;  .t - template spectrum
;  .tz - template spectrum shifted to solar frame
;  .sun - solar spectrum
; teminfo (structure[nm]) auxiliary information for each order
;  .m - echelle order, starting at 1 for first order in observation
;  .ib - index of first valid template point
;  .ie - index of last valid template point
;  .wb - wavelength of first valid template point
;  .we - wavelength of last valid template point
;  .file - name of file
;  .iter - last iteration number

;syntax
  if n_params() lt 2 then begin
    print, 'syntax: read_grand_tems, tems, teminfo [,mlist= ,/print]'
    return
  endif

;internal parameters
  sfx = 'tem'

;get information about grand output files in current directory.
  grand_files, gf
  if ~keyword_set(gf[sfx].uord) then begin
    message, /info, "no '" + sfx $
           + "' files produced by grand in current directory"
    return
  endif
  uord = gf[sfx].uord

;create complete order list, if user did not specify an order list
  if ~keyword_set(mlist) then mlist = fix(gf['tem'].uord)
  nm = n_elements(mlist)

;Allocate output info structure
  rec =          $
    { m    : 0   $
    , ib   : 0L  $
    , ie   : 0L  $
    , wb   : 0d0 $
    , we   : 0d0 $
    , file : ''  $
    , iter : ''  $
    }
  teminfo = replicate(rec, nm)

;read read template stellar spectrum for each requested order
  for im=0, nm-1 do begin
    m = mlist[im]
    iord = where(fix(gf[sfx].uord) eq m, count)
    if count eq 0 then begin
      message, 'no template file (.tem) for order ' + strtrim(m, 2)
    endif
    iter = gf[sfx].ilast[iord]
    file_tem = gf[sfx].flast[iord]
    read_grand_tem, file_tem, tem
    if im eq 0 then begin
      ntem = n_elements(tem)
      tems = replicate(tem[0], ntem, nm)
    endif
    tems[*,im] = tem
    teminfo[im].iter = iter
    teminfo[im].file = file_tem
  endfor

;determine first and last valid wavelength for each requested order
  for im=0, nm-1 do begin
    m = mlist[im]
    iwhr = where(tems[*,im].t ne tems[0,im].t, nwhr)
    if nwhr eq 0 then message, 'constant template for order ' + strtrim(m, 2)
    ibeg = iwhr[0]
    iwhr = where(tems[*,im].t ne tems[ntem-1,im].t, nwhr)
    iend = iwhr[nwhr-1]

;load output info structure
    teminfo[im].m  = m
    teminfo[im].ib = ibeg
    teminfo[im].ie = iend
    teminfo[im].wb = tems[ibeg,im].w
    teminfo[im].we = tems[iend,im].w

;print diagnostics
    if keyword_set(print) then begin
      if im eq 0 then print, '        ____PIXEL___   ___WAVELENGTH___'
      if im eq 0 then print, 'order   begin    end     begin      end'
      print, form='(3x,i2.2,1x,2i7,1x,2f9.2)' $
           , m, ibeg, iend, tems[ibeg,im].w, tems[iend,im].w
    endif
  endfor

end
