pro plot_grand_tems, mlist=mlist, print=print

;syntax
  if n_params() lt 0 then begin
    print, 'syntax: plot_grand_tems [,mlist= ,print]'
    return
  endif

;default value for optional parameters
  if n_elements(mlist) eq 0 then mlist = 0

;get current working directory
  cd, curr=cwd
  dir = file_basename(cwd)

;get information about grand output files in current directory.
  grand_files, gf

;read template data
  read_grand_tems, tems, ti, mlist=mlist, print=keyword_set(print)
  nm = n_elements(mlist)

;plot setup
  thi = 5
  chars = 1.3
  psfile = 'grand_tems.ps'
  set_plot, 'ps'
  !p.font = 0
  device, /iso, bits=8, file=psfile, /color
  device, xsize=10.0, ysize=7.5, /inch
  clist = c24([4,2,10])

;first plot all requested orders, if user requested more than one
  if nm ge 2 then begin
    xmin = min(ti.wb)
    xmax = max(ti.we)
    dx = xmax - xmin
    xr = [xmin-0.05*dx, xmax+0.05*dx]
    plot, tems.w, tems.t, /nodata $
        , xr=xr, xsty=1, yr=[0,1.05*max(tems.t)], ysty=1 $
        , xtit='Template Wavelength  (' + string(197b) + ')' $
        , ytit='Template Intensity' $
        , xmarg=[7,1], ymarg=[3.3,1.8] $
        , chars=chars, xthi=thi, ythi=thi
    for im=0, nm-1 do begin
      ib = ti[im].ib
      ie = ti[im].ie
      oplot, tems[ib:ie,im].w $
           , tems[ib:ie,im].t $
           , co=clist[im mod 2]
    endfor
  endif

;label the orders
  xcr = !x.crange
  dxcr = xcr[1] - xcr[0]
  ycr = !y.crange
  dycr = ycr[1] - ycr[0]
  for im=0, nm-1 do begin
    ib = ti[im].ib
    ie = ti[im].ie
    xyouts, 0.5*(tems[ib,im].w + tems[ie,im].w) $
          , ycr[0]+0.03*dycr, align=0.5, size=chars $
          , co=clist[im mod 2], string(mlist[im], form='(i2.2)')
  endfor
  imin = min(ti.iter, max=imax)
  if imin eq '99' then imin = min(gf['iter'].last)
  if imax eq '99' then imax = max(gf['iter'].last)
  istr = strtrim(imin, 2)
  if imin ne imax then istr += ' to ' + strtrim(imax, 2)
  xyouts, xcr[0]+0.02*dxcr, ycr[1]+0.015*dycr $
        , dir + ', ' + istr + ' iterations'

;plot each requested order on a separate page
  dm = [-1, 1, 0]				;sequence of order plots
  for im=0, nm-1 do begin
    m = mlist[im]

;determine plot range
    xmin = ti[im].wb
    xmax = ti[im].we
    dx = xmax - xmin
    xr = [xmin-0.05*dx, xmax+0.05*dx]
    iplt = where(tems[*,im].w ge xr[0] and tems[*,im].w le xr[1])

;create plot window for current order
    plot, tems[iplt,im].w, tems[iplt,im].t, /nodata $
        , xr=xr, xsty=1, yr=[0,1.05*max(tems.t)], ysty=1 $
        , xtit='Template Wavelength  (' + string(197b) + ')' $
        , ytit='Template Intensity' $
        , xmarg=[7,1], ymarg=[3.3,1.8] $
        , chars=chars, xthi=thi, ythi=thi

;label the orders
    xcr = !x.crange
    dxcr = xcr[1] - xcr[0]
    ycr = !y.crange
    dycr = ycr[1] - ycr[0]
    iter = ti[im].iter
    if iter eq '99' then iter = gf['iter'].last[im]
    xyouts, xcr[0]+0.02*dxcr, ycr[1]+0.015*dycr $
          , 'Order ' + string(m, form='(i2.2)') $
          + ',  ' + dir + '!e !n/!e !n' + ti[im].file $
          + ' (' + iter + ' iterations)'

;plot current order
    for i=0, 2 do begin
      jm = im + dm[i]
      if jm lt 0 or jm ge nm then continue
      oplot, tems[iplt,jm].w, tems[iplt,jm].t $
           , co=clist[i mod 3]
      xyouts, 0.5*((xr[0] > ti[jm].wb) + (ti[jm].we < xr[1])) $
            , ycr[0]+0.03*dycr, align=0.5, size=chars $
            , co=clist[i mod 3], string(mlist[jm], form='(i2.2)')
    endfor

;end of loop through roders.
  endfor
  
;plot cleanup
  device, /close
  set_plot, 'x'
  !p.font = -1
  print, 'writing ' + psfile

end
