pro plot_jump_phase, imsum

  if n_params() lt 1 then begin
    print, 'syntax: flat_jump_phase, imsum'
    print, "  e.g.: fits_read, 'j550225.fits', im, hd, exten=2"
    print, '        imsum = total(float(im), 1)'
    print, '        flat_jump_phase, imsum'
    return
  endif

;internal parameters.
  ;per = 1024 / 15.0
  ;dx = 39.37
  per = 1024.0
  dx = 654.

  xr = [-6, 6]
  yr = [0.985, 1.03]

;remove all but the smallest scale features
  y = imsum / median(imsum, 7)

;construct phased response
  nx = n_elements(y)
  x = findgen(nx) - dx
  dpix = (x + per) mod per
  iwhr = where(dpix ge per/2, nwhr)
  if nwhr gt 0 then dpix[iwhr] -= per

;construct median response per phased pixel
  nmed = xr[1] - xr[0] + 1
  xmed = xr[0] + findgen(nmed); + 0.5
  ymed = fltarr(nmed)
  for imed=0, nmed-1 do begin
    cen = xmed[imed]
    iwhr = where(abs(dpix - xmed[imed]) le 0.5, nwhr)
    ymed[imed] = median(y[iwhr])
  endfor

;calculate mean of medians
  ymedsum = total(ymed) / nmed
;  ymedsum = 1.0

;plot parameters
  thi = 5
  chars = 1.8
  xtit = '(Column - ' $
       + strtrim(string(dx, form='(f9.2)'), 2) $
       + ') mod ' $
       + strtrim(string(per, form='(f9.2)'), 2)

;open postscript file
  psfile = 'plot_jump_phase.eps'
  set_plot, 'ps'
  !p.font = 0
  device, /iso, bits=8, file=psfile, /color, /encap
  device, xsize=8.0, ysize=6.0, xoff=0.25, yoff=2.5, /inch, /port

;plot window
  plot, xr, yr, /nodata $
      , xr=xr, /xsty, yr=yr, /ysty $
      , xthi=thi, ythi=thi $
      , xticks=8 $
      , xtit=xtit $
      , ytit='Flat Field Response' $
      , chars=chars $
      , xmarg=[9,1], ymarg=[3.2,0.5]

;annotation coordinates
  xcr = !x.crange
  dxcr = xcr[1] - xcr[0]
  ycr = !y.crange
  dycr = ycr[1] - ycr[0]

;plot data
; oplot, xcr, [1,1], thi=thi, co=c24(2)
  oplot, xcr, 1.0+[0,0], thi=thi, lin=2, co=c24(10)
  oplot, [0,0], !y.crange, thi=thi, co=c24(4)
  oplot, dpix, y, ps=7, thi=thi
  oplot, xmed, ymed, psy=10, thi=2*thi, co=!p.background
  oplot, xmed, ymed, psy=10, thi=thi, co=c24(10)

;annotations
  xyouts, xcr[0]+0.96*dxcr, ycr[0]+0.06*dycr, siz=chars $
        , align=1, co=c24(10) $
        , 'Integral!e !n=!e !n' + strtrim(string(ymedsum, form='(f19.6)'), 2)

;close postscript file.
  device, /close
  set_plot, 'x'
  !p.font = -1
  print, psfile + ' closed'

end
