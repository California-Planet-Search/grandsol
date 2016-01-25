pro plot_simvel, oblist, ps=ps

;Syntax
  if n_params() lt 0 then begin
    print, 'syntax: plot_simvel [,obslist]'
    return
  endif

;Constants
  c = 2.99792458d8

;Default value for optional parameters
  if ~keyword_set(obslist) then obslist = 'obslist'

;Read observation list file.
  read_grand_obslist, obli, obslist

;Read information from simulated observations.
  read_siminfo, obli.obsdir, obli.obs.file, info
  vinp = info.velsum

;Find velocity output files.
  fpatt = '*.[0-9][0-9].vel'
  files = file_search(fpatt, count=nfile)
  if nfile eq 0 then begin
    print, 'no velocity saves files in current directory'
    return
  endif

;Figure out the 
  file_vel = files[nfile-1]	;assumes final velocities in last file

stop

;Read velocities that grand wrote to disk.
  data = dblarr(5, nobs)
  openr, unit, file_vel, /get_lun
  readf, unit, data
  free_lun, unit
  zn = reform(data[1,*])	;measured bulk shift of entire order
  z0 = reform(data[2,*])	;input barycentric shift
  zbarn = reform(data[3,*])	;mean shift for lines after sigma clipping
  zsign = reform(data[4,*])	;std dev of line shifts after sigma clipping

;calculate velocities
  v0 = c * ((1+z0)^2 - 1) / ((1+z0)^2 + 1)
  veln = c * ((1+zn)^2 - 1) / ((1+zn)^2 + 1) - v0
  vbarn = c * ((1+zbarn)^2 - 1) / ((1+zbarn)^2 + 1) - v0
  vinp = c * ((1+zplus1-1)^2 - 1) / ((1+zplus1-1)^2 + 1) - v0

;fit velocities
  coef = poly_fit(vinp, vbarn, 1, double)

;plot measured velocity versus input velocity
  thi = 1
  chars = 1.8
  if keyword_set(ps) then begin
    ps_open, 'plot_simvel', /ps, /encap, /color
    device, xsize=6.5, ysize=6.5, /inch
    thi = 3
    chars = 1.3
  endif
  vmax = max(abs(vinp)) > max(abs(vbarn)) > max(abs(veln))
  vr = [-vmax, vmax]
  psym, 'circle'
  plot, vinp, vbarn, ps=8, /nodata $
      , xr=vr, xsty=3, yr=vr, ysty=3 $
      , xtit='Input Orbital Velocity  (m/s)' $
      , ytit='Measured Orbital Velocity  (m/s)' $
      , xmarg=[8,1], ymarg=[3.5,0.5] $
      , chars=chars, xthi=thi, ythi=thi
  oplot, vr, vr, thi=2*thi, lin=2, co=c24(10)
  oplot, vr, poly(vr, coef), thi=2*thi, co=c24(2)
  oplot, vinp, vbarn, ps=8, thi=0.67*thi, symsiz=1.5

  xcr = !x.crange
  dxcr = xcr[1] - xcr[0]
  ycr = !y.crange
  dycr = ycr[1] - ycr[0]
  xyouts, xcr[0]+0.1*dxcr, ycr[0]+0.9*dycr, size=1.2*chars, co=c24(2) $
        , 'Slope = ' + strtrim(string(coef[1], form='(f19.3)'), 2)

  if keyword_set(ps) then ps_close

end
