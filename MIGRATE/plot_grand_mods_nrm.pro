pro plot_grand_mods_nrm, mlist=mlist, print=print $
                       , mods=mods, modinfo=mi

;syntax
  if n_params() lt 0 then begin
    print, 'syntax: plot_grand_mod [,mlist=, mods=, modinfo= ,/print]'
    return
  endif

;default value for optional parameters
  if n_elements(mlist) eq 0 then mlist = 0

;get current working directory
  cd, curr=cwd
  dir = file_basename(cwd)

;get information about grand output files in current directory.
  grand_files, gf

;get locations of nodes for normalization vector
  if max(gf.keys() eq 'par') eq 0 then message, 'no .par files'

;read model fit data [expensive], if not passed in by argument.
  if ~keyword_set(mods) or ~keyword_set(mi) then begin
    read_grand_mods, mods, mi, mlist=mlist, print=keyword_set(print)
  endif
  if ~keyword_set(mlist) then mlist = mi.m
  sz = size(mods)
  npix = sz[1]
  nm   = sz[2]
  nobs = sz[3]
  xpix = 1 + findgen(npix)
  yobs = 1 + findgen(nobs)

;calculate mean over all observations
  mean_nrm = total(mods.nrm, 3, /double) / nobs

;;
;; PLOT MEAN NORMALIZATION VECTOR FOR ALL ORDERS TOGETHER
;;

;plot setup
  xsav = !x
  ysav = !y
  psav = !p
  !p.font = 0
  !p.multi = 0
  psfile = 'grand_mods_nrm.ps'
  set_plot, 'ps'
  device, /iso, bits=8, file=psfile, /color
  device, xsize=10.0, ysize=7.5, /inch
  rainbow_colors, clist, nm
  psym, 'circle', /fill
  chars = 1.3
  thi = 5

;make first plot window that will show all orders
  xr = [0, npix]
  yr = minmax(mean_nrm)
  ytit = 'Normalization,!e !n averaged over ' $
       + strtrim(nobs, 2) + ' observations'
  plot, xr, yr, /nodata $
      , xr=xr, xsty=3, yr=yr, ysty=1 $
      , xtit='Pixels' $
      , ytit=ytit $
      , xmarg=[7,1], ymarg=[3.3,1.8] $
      , chars=chars, xthi=thi, ythi=thi
  xcr = !x.crange
  dxcr = xcr[1] - xcr[0]
  ycr = !y.crange
  dycr = ycr[1] - ycr[0]

;plot the normalization vector for each order, averaged over all observations
  for im=0, nm-1 do begin
    m = mlist[im]
    oplot, xpix, mean_nrm[*,im], thi=thi, co=clist[im]
    xyouts, xcr[0]+(im+2)/(nm+2.0)*dxcr, ycr[0]+0.93*dycr $
          , align=0.5, size=chars, col=clist[im], string(m, form='(i2.2)')
  endfor
  xyouts, xcr[0]+1/(nm+2.0)*dxcr, ycr[0]+0.93*dycr $
          , align=0.5, size=chars, 'Ord'

;annotate with title that includes the iteration number
  imin = min(mi.iter, max=imax)
  if imin eq '99' then imin = min(gf['iter'].last)
  if imax eq '99' then imax = max(gf['iter'].last)
  istr = strtrim(imin, 2)
  if imin ne imax then istr += ' to ' + strtrim(imax, 2)
  xyouts, xcr[0]+0.02*dxcr, ycr[1]+0.015*dycr, chars=chars $
        , dir + ', ' + istr + ' iterations'

;loop through requested orders
  for im=0, nm-1 do begin
    m = mi[im].m
    nrmi = reform(mods[*,im,*].nrm)

;read parameters for current order
    iwhr = where(gf['par'].uord eq m, nwhr)
    if nwhr ne 1 then message, 'error finding .par file for current order'
    file_par = gf['par'].flast[iwhr]
    read_grand_par, file_par[0], par

;get barycentric range
    wob1 = mods[floor(npix/2.0)  ,im,nobs/2].wob
    wob2 = mods[floor(npix/2.0)+1,im,nobs/2].wob
    wmid = 0.5 * (wob2 + wob1)
    dwob = wob2 - wob1
    dbc = par.maxbc - par.minbc
    dbcpix = (wmid / dwob) * (dbc / 2.9979246e8)

;get node locations for normalization vector
    xnrm = par.xnrm
    nnrm = n_elements(xnrm)
    nnode = 10

;;
;; PLOT MEAN NORMALIZATION VECTOR FOR ALL ORDERS TOGETHER
;;

;load color map with gray scale and then red, green, and blue
    tvlct, [bindgen(253),0b,0b,255b] $
         , [bindgen(253),0b,255b,0b] $
         , [bindgen(253),255b,0b,0b]

;display image of normalization for current order
    med = median(nrmi)
    zhw = 5.0 * median(abs(nrmi-med))
    !x.margin = [7,2]
    !y.margin = [3,1.8]
    !p.charsize = chars
    !p.multi = 0
    display, nrmi, xpix, yobs $
           , min=med-zhw, max=med+zhw $
           , xtit='Pixel' $
           , ytit='Observation Number'

;computed plot limits
    xcr = !x.crange
    dxcr = xcr[1] - xcr[0]
    ycr = !y.crange
    dycr = ycr[1] - ycr[0]

;oplot a grid angled with the stellar frame
    xstep = 100
    for x=xstep, npix, xstep do begin
      col = median(nrmi[x-1,*]) lt med ? 252b : 1b
;     oplot, [x,x], ycr[0]+dycr*[0.03,0.97], lin=1, thi=1, col=col
      oplot, x+[0.5,-0.5]*dbcpix*0.94, ycr[0]+dycr*[0.03,0.97] $
           , lin=1, thi=1, col=col
;     oplot, x-[0.5,-0.5]*dbcpix*0.94, ycr[0]+dycr*[0.03,0.97] $
;          , lin=1, thi=1, col=col
    endfor

;make a 2D plot for the middle observation
    y = nrmi[*,nobs/2]
    ymed = median(y)
    yscale = 0.4 * dycr / max(y-ymed)
    yy = ycr[0] + 0.5 * dycr + yscale * (y - ymed)
    oplot, xpix, yy, thi=thi, co=255b
    oplot, !x.crange, ycr[0]+0.5*dycr+[0,0], lin=1, thi=2, col=255b

;plot 1% indicator
    hw1p = 0.01 * ymed
    x1p = xstep * (round(0.5 * npix / xstep) + 0.5) + [0,0]
    y1p = interpol(yy, xpix, x1p) + yscale * hw1p * [-1,1]
    oplot, x1p, y1p, thi=thi, co=255b
    oplot, x1p[0]+[-1,1]*0.006*dxcr, y1p[0]+[0,0], thi=thi, co=255b
    oplot, x1p[0]+[-1,1]*0.006*dxcr, y1p[1]+[0,0], thi=thi, co=255b

;label the current order
    imin = min(mi.iter, max=imax)
    if imin eq '99' then imin = min(gf['iter'].last)
    if imax eq '99' then imax = max(gf['iter'].last)
    istr = strtrim(imin, 2)
    xyouts, xcr[0]+0.02*dxcr, ycr[1]+0.015*dycr, chars=chars $
          , 'nrmi,  order ' + string(m, form='(i2.2)') $
          + ',  ' + dir + '!e !n/!e !n' + mi[im].file $
          + ' (' + istr + ' iterations)'

;;
;; PLOT LEFT EDGE OF CURRENT ORDER
;;

;extract segment to plot
    no = 5
    jo = round((nobs - 1) * findgen(no) / (no - 1))
    ib = 0
    ie = ceil(0.5*(xnrm[nnode-1]+xnrm[nnode])) - 1
    x = xpix[ib:ie,jo]
    y = nrmi[ib:ie,jo]
    ymin = min(y, max=ymax)
    yrlo = ymin - 0.10 * (ymax-ymin)
    yrhi = ymax + 0.10 * (ymax-ymin)

;plot both ends of the current order
    !p.multi = [0,2,1]
    xr = [0,max(x)]
    yr = [yrlo,yrhi]
    plot, xr, yr, /nodata $
        , xr=xr, /xsty, yr=yr, /ysty $
        , xtit='Pixel' $
        , ytit='Normalization Vector' $
        , chars=chars, xthi=thi, ythi=thi $
        , xmarg=[9,1], ymarg=[3,1.8]

;computed plot limits
    xcr = !x.crange
    dxcr = xcr[1] - xcr[0]
    ycr = !y.crange
    dycr = ycr[1] - ycr[0]

;plot data
    rainbow_colors, clist, no
    for io=0, no-1 do begin
      oplot, x, y[*,io], thi=thi, psym=3, col=clist[io]
      xnode = xnrm[0:nnode-1]
      oplot, xnode, nrmi[xnode-1,jo[io]], psym=8 $
           , symsiz=0.5, col=clist[io], /noclip
    endfor

;label the current order
    xyouts, xcr[0]+0.02*dxcr, ycr[1]+0.015*dycr, chars=chars $
          , 'nrmi,  order ' + string(m, form='(i2.2)') $
          + ',  ' + dir + '!e !n/!e !n' + mi[im].file $
          + ' (' + istr + ' iterations)'

;;
;; PLOT RIGHT EDGE OF CURRENT ORDER
;;

;extract segment to plot
    ib = ceil(0.5*(xnrm[nnrm-nnode-1]+xnrm[nnrm-nnode])) - 1
    ie = npix-1
    x = xpix[ib:ie,jo]
    y = nrmi[ib:ie,jo]
    ymin = min(y, max=ymax)
    yrlo = ymin - 0.10 * (ymax-ymin)
    yrhi = ymax + 0.10 * (ymax-ymin)

;plot both ends of the current order
    xr = [min(x),npix]
    yr = [yrlo,yrhi]
    plot, xr, yr, /nodata $
        , xr=xr, /xsty, yr=yr, /ysty $
        , xtit='Pixel' $
        , ytit='Normalization Vector' $
        , chars=chars, xthi=thi, ythi=thi $
        , xmarg=[9,1], ymarg=[3,1.8]

;computed plot limits
    xcr = !x.crange
    dxcr = xcr[1] - xcr[0]
    ycr = !y.crange
    dycr = ycr[1] - ycr[0]

;plot data
    rainbow_colors, clist, no
    for io=0, no-1 do begin
      oplot, x, y[*,io], thi=thi, psym=3, col=clist[io]
      xnode = xnrm[nnrm-nnode:nnrm-1]
      oplot, xnode, nrmi[xnode-1,jo[io]], psym=8 $
           , symsiz=0.5, col=clist[io], /noclip
      xyouts, xcr[0]+(io+2.0)/(no+2)*dxcr, ycr[0]+0.94*dycr $
            , align=0.5, size=chars, col=clist[io] $
            , strtrim(round(yobs[jo[io]]), 2)
    endfor
    xyouts, xcr[0]+1.0/(no+2)*dxcr, ycr[0]+0.94*dycr $
          , align=0.5, size=chars, 'Obs'

;end of loop through roders.
  endfor

;plot cleanup
  device, /close
  set_plot, 'x'
  !p.font = -1
  print, 'writing ' + psfile
  !x = xsav
  !y = ysav
  !p = psav

end
