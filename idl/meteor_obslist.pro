pro meteor_obslist, obsid_file, print=print, plot=plot

;syntax
  if n_params() lt 1 then begin
    print, 'syntax: meteor_obslist, obslist [/print, /plot]'
    print, "syntax: meteor_obslist, 'obslist', /print"
    return
  endif

;default values for optional parameters
  if n_elements(plot) eq 0 then plot = 0

;read list of observations identifiers from an obslist file
  if keyword_set(print) then print, 'reading ' + obsid_file
  read_grand_obslist, obli, obslist
  obsids = obli.obs.file
  nobs = n_elements(obsids)

;open output file
  fmt = '(a-10,f7.2,16i5,f6.1,f6.3)'
  file = 'meteor'
  if keyword_set(print) then print, 'writing ' + file
  openw, unit, file, /get_lun
  for iobs=0L, nobs-1 do begin
    obsid = obsids[iobs]
    meteor, obsid, meteor, zmed, mar, y7at2k, plot=plot, errmsg=errmsg
    if keyword_set(meteor) then begin
      if keyword_set(print) then begin
        print, form=fmt, obsid, y7at2k, meteor, zmed, mar
      endif
      printf, unit, form=fmt, obsid, y7at2k, meteor, zmed, mar
    endif else begin
      print, obsid + ' ' + errmsg
      printf, unit, obsid + ' ' + errmsg
    endelse
  endfor
  free_lun, unit

end
