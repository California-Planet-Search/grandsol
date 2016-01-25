pro plot_grand_iterate, nord=nord, yr=yr

;syntax
  if n_params() lt 0 then begin
    print, 'syntax: plot_grand_iterate [,nord= ,yr=]'
    return
  endif

;default value for optional parameters
  if ~keyword_set(yr) then yr = [0,0]

;constants
  c = 2.99792458d8

;read velocities from observation lists
  files = file_search('obslist*', count=nfile)
  for ifile=0, nfile-1 do begin
    read_grand_obslist, obli, files[ifile]
    if ifile eq 0 then begin
      nobs = n_elements(obli.obs)
      niter = nfile - 1
      vorbs = dblarr(nobs, niter)
      vtrue = obli.obs.vorb
      isort = sort(vtrue)
      vtrue = vtrue[isort]
      bc = obli.obs.barycorr
    endif else begin
      vorbs[*,ifile-1] = obli.obs[isort].vorb
    endelse
  endfor

;errors
  sd = fltarr(niter)
  sdstr = strarr(niter)
  for iiter=0, niter-1 do begin
    sd[iiter] = stddev(vorbs[*,iiter] - vtrue)
    sdstr[iiter] = 'Iter ' + strtrim(iiter, 2) + ':  !ms!x=' $
                 + strtrim(string(sd[iiter], form='(f9.2)'), 2) + ' m/s'
  endfor

;plot setup
  thi = 5
  chars = 1.3
  psfile = 'iterate_grand.ps'
  set_plot, 'ps'
  !p.font = 0
  device, /iso, bits=8, file=psfile, /color
  device, xsize=10.0, ysize=7.5, /inch, /land

;plot window
  plot, vtrue, vorbs[*,1]-vtrue, /nodata $
      , xsty=3, yr=yr, ysty=3 $
      , xtit='Input Orbital Velocity  (m/s)' $
      , ytit='Measured - Input Orbital Velocity  (m/s)' $
      , xmarg=[8,2], ymarg=[3.7,1] $
      , chars=chars, xthi=thi, ythi=thi
  xcr = !x.crange
  dxcr = xcr[1] - xcr[0]
  ycr = !y.crange
  dycr = ycr[1] - ycr[0]

;plot sequence of orbital velocities.
  clist = c24([2,7,10,4,11])
  for iiter=1, niter-1 do begin
    col = clist[(iiter-1) mod 5]
    oplot, vtrue, vorbs[*,iiter]-vtrue, col=col, thi=thi
    xyouts, xcr[0]+0.65*dxcr, ycr[0]+(0.96-0.04*iiter)*dycr $
          , size=chars, sdstr[iiter], col=col
  endfor

;plot window
  psym, 'circle', thi=thi
  plot, vtrue, vorbs[*,niter-1]-vtrue, syms=2, ps=8 $
      , xsty=3 $
      , xtit='Input Orbital Velocity  (m/s)' $
      , ytit='Measured - Input Orbital Velocity  (m/s)' $
      , xmarg=[9,2], ymarg=[3.4,2.2] $
      , chars=chars, xthi=thi, ythi=thi
  xcr = !x.crange
  dxcr = xcr[1] - xcr[0]
  ycr = !y.crange
  dycr = ycr[1] - ycr[0]
  xyouts, xcr[0]+0.5*dxcr, ycr[1]+0.018*dycr, align=0.5 $
        , size=chars, sdstr[niter-1]

;plot window
  psym, 'circle', thi=thi
  plot, bc/1e3, vorbs[*,niter-1]-vtrue, syms=2, ps=8 $
      , xsty=3 $
      , xtit='Barycentric Correction  (km/s)' $
      , ytit='Measured - Input Orbital Velocity  (m/s)' $
      , xmarg=[9,2], ymarg=[3.4,2.2] $
      , chars=chars, xthi=thi, ythi=thi
  xcr = !x.crange
  dxcr = xcr[1] - xcr[0]
  ycr = !y.crange
  dycr = ycr[1] - ycr[0]
  xyouts, xcr[0]+0.5*dxcr, ycr[1]+0.018*dycr, align=0.5 $
        , size=chars, sdstr[niter-1]

;plot cleanup
  device, /close
  set_plot, 'x'
  !p.font = -1
  print, 'writing ' + psfile

;save
  if 1 then begin
    dumpfile = 'dump_vels.sav'
    print, 'saving ' + dumpfile
    save, file=dumpfile, vtrue, vorbs, bc, sd, sdstr
  endif

end
