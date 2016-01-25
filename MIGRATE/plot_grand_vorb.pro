pro plot_grand_vorb, oblist, errbar=errbar, nord=nord $
                   , vorbi=vorbi, vorbm=vorbm

;syntax
  if n_params() lt 0 then begin
    print, 'syntax: plot_grand_vel [,obslist ,/errbar ,vorbi= ,vorbm=]'
    return
  endif

;constants
  c = 2.99792458d8
  sfx = 'vel'

;default value for optional parameters
  if ~keyword_set(obslist) then begin
    obslists = file_search('obslist*', count=nobslist)
    if nobslist eq 0 then message, 'could not identify obslist file'
    obslist = obslists[nobslist-1]
  endif

;default value for optional parameters
  if obslist eq 'obslist' then begin
    obslist0 = obslist
  endif else begin
    obslist0 = file_search('obslist', count=nobslist)
    if nobslist eq 0 then message, 'could not identify obslist0 file'
    if nobslist ge 2 then message, 'multiple obslist0 files'
    obslist0 = obslist0[0]
  endelse

;read initial velocity
  read_grand_obslist, obli0, obslist0
  vorb0 = obli0.obs.vorb

;get current working directory
  cd, curr=cwd
  dir = file_basename(cwd)

;read obslist
  read_grand_obslist, obli, obslist
  bc = obli.obs.barycorr
  vorb = obli.obs.vorb

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

;calculate median absolute error
    mad_vbarn = median(abs(dvbarn-vorb))
    mad_veln = median(abs(dveln-vorb))

;plot setup
    if iord eq iord_beg then begin
      thi = 5
      chars = 1.3
      psfile = 'grand_vorb.ps'
      set_plot, 'ps'
      !p.font = 0
      device, /iso, bits=8, file=psfile, /color
      device, xsize=10.0, ysize=7.5, /inch

;x plot limits
      xmin = min(vorb0, max=xmax)
      xrlo = xmin - 0.05 * (xmax-xmin)
      xrhi = xmax + 0.05 * (xmax-xmin)
    endif

;y plot limits for each echelle order
    ymin = min([dvbarn, dveln], max=ymax)
    yrlo = ymin - 0.10 * (ymax-ymin)
    yrhi = ymax + 0.10 * (ymax-ymin)

;use same plot limits for x and y
    xrlo = xrlo < yrlo
    xrhi = xrhi > yrhi
    yrlo = xrlo
    yrhi = xrhi

;plot window
    plot, [xrlo,xrhi], [yrlo,yrhi], ps=3, /nodata $
        , xr=[xrlo,xrhi], xsty=1, yr=[yrlo,yrhi], ysty=1 $
        , xtit='Input Orbital Velocity  (m/s)' $
        , ytit='Measured Orbital Velocity  (m/s)' $
        , xmarg=[9,1], ymarg=[3.3,3.0] $
        , chars=chars, xthi=thi, ythi=thi
    xcr = !x.crange
    dxcr = xcr[1] - xcr[0]
    ycr = !y.crange
    dycr = ycr[1] - ycr[0]

;velocities
;   oplot, xcr, mn_vbarn+[0,0], thi=thi, lin=1, co=c24(2)
;   oplot, xcr, mn_veln+[0,0], thi=thi, lin=1, co=c24(4)
    ran = [min(vorb0),max(vorb0)]
    oplot, ran, ran, thi=hti, lin=1
    psym, 'circle', thi=0.67*thi
    oplot, vorb0, vorb+dvbarn, ps=8, symsiz=1.5, col=c24(2)
    psym, 'square', thi=0.67*thi
    oplot, vorb0, vorb+dveln, ps=8, symsiz=1.5, col=c24(4)

;uncertainties
    if keyword_set(errbar) then begin
      if iord eq -1 then begin
        for iobs=0, nobs-1 do begin
          oplot, vorb0[iobs]+[0,0] $
               , dvbarn[iobs]+[-1,1]*sdvel[iobs].vbarn $
               , thi=thi, col=c24(2)
          oplot, vorb0[iobs]+[0,0] $
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
          , 'median |veln-vorb|!e !n=!e !n' $
          + strtrim(string(mad_veln, form='(f19.1)'), 2) $
          + ' m/s', align=1
    xyouts, xcr[0]+0.98*dxcr, ycr[1]+0.015*dycr, size=chars, co=c24(2) $
          , 'median |vbarn-vorb|!e !n=!e !n' $
          + strtrim(string(mad_vbarn, form='(f19.1)'), 2) $
          + ' m/s', align=1

;end of loop through roders.
  endfor
  
;plot cleanup
  device, /close
  set_plot, 'x'
  !p.font = -1
  print, 'writing ' + psfile

;load return argument
  vorbi = vorb
  vorbm = dvbarn

end
