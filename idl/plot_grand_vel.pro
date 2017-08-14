pro plot_grand_vel, oblist, errbar=errbar, nord=nord

;syntax
  if n_params() lt 0 then begin
    print, 'syntax: plot_grand_vel [,obslist ,/errbar]'
    return
  endif

;constants
  c = 2.99792458d8
  sfx = 'vel'

;default value for optional parameters
  if ~keyword_set(obslist) then obslist = 'obslist'

;get current working directory
  cd, curr=cwd
  dir = file_basename(cwd)

;get information about grand output files in current directory.
  grand_files, gf
  if ~keyword_set(gf[sfx].uord) then begin
    message, /info, "no '" + sfx $
           + "' files produced by grand in current directory"
    return
  endif

;read velocities for all orders
  if ~keyword_set(nord) then nord = n_elements(gf[sfx].uord)
  for iord=0, nord-1 do begin
    file_vel = gf[sfx].flast[iord]
    read_grand_vel, file_vel, vel
    if iord eq 0 then begin
      nobs = n_elements(vel)
      vels = replicate(vel[0], nobs, nord)
    endif
    vels[*,iord] = vel
  endfor

;calculate mean of all orders
  if nord gt 1 then begin
    iord_beg = -1
    mnvel = vel
    sdvel = vel
    ntag = n_elements(tag_names(vel))
    for iobs=0, nobs-1 do begin
      for itag=0, ntag-1 do begin
        mnvel[iobs].(itag) = mean(vels[iobs,*].(itag), /double)
        sdvel[iobs].(itag) = stddev(vels[iobs,*].(itag), /double)
      endfor
    endfor
  endif else begin
    iord_beg = 0
  endelse

;loop through orders
  for iord=iord_beg, nord-1 do begin
    ord = strtrim(fix(gf[sfx].olast[iord]), 2)
    iter = gf[sfx].ilast

;data for current order (or for mean of all orders)
    if iord eq -1 then begin
      file_vel = ''
      vel = mnvel
    endif else begin
      file_vel = gf[sfx].flast[iord]
      vel = vels[*,iord]
    endelse

;read Z values determined by the grand solution.
    dvbarn = vel.vbarn - vel.v0
    dveln = vel.veln - vel.v0

;calculate mean
    mn_vbarn = mean(dvbarn)
    mn_veln = mean(dveln)

;calculate standard deviation
    sd_vbarn = stddev(dvbarn)
    sd_veln = stddev(dveln)

;plot setup
    if iord eq iord_beg then begin
      thi = 5
      chars = 1.3
      psfile = 'grand_vel.ps'
      set_plot, 'ps'
      !p.font = 0
      device, /iso, bits=8, file=psfile, /color
      device, xsize=10.0, ysize=7.5, /inch

;x plot limits
      xmin = min(vels.v0, max=xmax)
      xrlo = xmin - 0.05 * (xmax-xmin)
      xrhi = xmax + 0.05 * (xmax-xmin)
    endif

;y plot limits for each echelle order
    ymin = min([dvbarn, dveln], max=ymax)
    yrlo = ymin - 0.10 * (ymax-ymin)
    yrhi = ymax + 0.10 * (ymax-ymin)

;plot window
    plot, [xrlo,xrhi], [yrlo,yrhi], ps=3, /nodata $
        , xr=[xrlo,xrhi]/1e3, xsty=1, yr=[yrlo,yrhi], ysty=1 $
        , xtit='Barycentric Velocity  (km/s)' $
        , ytit='Measured Orbital Velocity  (m/s)' $
        , xmarg=[9,1], ymarg=[3.3,3.0] $
        , chars=chars, xthi=thi, ythi=thi
    xcr = !x.crange
    dxcr = xcr[1] - xcr[0]
    ycr = !y.crange
    dycr = ycr[1] - ycr[0]

;velocities
    oplot, xcr, mn_vbarn+[0,0], thi=thi, lin=1, co=c24(2)
    oplot, xcr, mn_veln+[0,0], thi=thi, lin=1, co=c24(4)
    psym, 'circle', thi=0.67*thi
    oplot, vel.v0/1d3, dvbarn, ps=8, symsiz=1.5, col=c24(2)
    psym, 'square', thi=0.67*thi
    oplot, vel.v0/1d3, dveln, ps=8, symsiz=1.5, col=c24(4)

;uncertainties
    if keyword_set(errbar) then begin
      if iord eq -1 then begin
        for iobs=0, nobs-1 do begin
          oplot, vel[iobs].v0/1d3+[0,0] $
               , dvbarn[iobs]+[-1,1]*sdvel[iobs].vbarn $
               , thi=thi, col=c24(2)
          oplot, vel[iobs].v0/1d3+[0,0] $
               , dveln[iobs]+[-1,1]*sdvel[iobs].veln $
               , thi=thi, col=c24(4)
        endfor
      endif
    endif

;annotations
    if iord eq -1 then begin
      xyouts, xcr[0]+0.02*dxcr, ycr[1]+0.015*dycr, size=chars $
            , 'Mean, ' + strtrim(nord, 2) + ' orders, ' $
            + dir
    endif else begin
      xyouts, xcr[0]+0.02*dxcr, ycr[1]+0.015*dycr, size=chars $
            , 'Order ' + ord + ',  ' + dir + '!e !n/!e !n' + file_vel
    endelse
    xyouts, xcr[0]+0.98*dxcr, ycr[1]+0.050*dycr, size=chars, co=c24(4) $
          , '!ms!x veln!e !n=!e !n' $
          + strtrim(string(sd_veln, form='(f19.3)'), 2) $
          + ' m/s', align=1
    xyouts, xcr[0]+0.98*dxcr, ycr[1]+0.015*dycr, size=chars, co=c24(2) $
          , '!ms!x vbarn!e !n=!e !n' $
          + strtrim(string(sd_vbarn, form='(f19.3)'), 2) $
          + ' m/s', align=1

;end of loop through roders.
  endfor
  
;plot cleanup
  device, /close
  set_plot, 'x'
  !p.font = -1
  print, 'writing ' + psfile

end
