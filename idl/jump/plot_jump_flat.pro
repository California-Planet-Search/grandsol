pro plot_jump_flat, imsum

  if n_params() lt 1 then begin
    print, 'syntax: plot_jump_flat, imsum'
    print, "  e.g.: fits_read, 'j550225.fits', im, hd, exten=2"
    print, '        imsum = total(float(im), 1)'
    print, '        plot_jump_flat, imsum'
    return
  endif

;get jump locations
  jump, xjmp

;remove all but the smallest scale features
  nx = n_elements(imsum)
  x = indgen(nx)
  c = median(imsum, 7)
  y = imsum / c

;plot parameters
  thi = 5
  chars = 1.8
  xr = [0, 4000]
  yr1 = [2.3e6, 4e6]
  yr2 = [0.98499, 1.02501]
  xtit = 'Pixel'

;open postscript file
  psfile = 'plot_jump_flat.eps'
  set_plot, 'ps'
  !p.font = 0
  device, /iso, bits=8, file=psfile, /color, /encap
  device, xsize=8.0, ysize=6.0, xoff=0.25, yoff=2.5, /inch, /port

;plot setup
  !p.multi = [0, 1, 2]
  !y.omargin = [2.9, 0]

;plot window
  plot, xr, yr1, /nodata $
      , xr=xr, /xsty, yr=yr1, /ysty $
      , xthi=thi, ythi=thi $
      , xtickn=replicate(' ', 19) $
      , ytickn=['2.5','3.0','3.5','4.0'] + string(215b) + '10!a6!n' $
      , ytit='Mashed Flat (ADU)' $
      , chars=chars $
      , xmarg=[10.5,2.2], ymarg=[0.0,0.5]

;annotation coordinates
  xcr = !x.crange
  dxcr = xcr[1] - xcr[0]
  ycr = !y.crange
  dycr = ycr[1] - ycr[0]

;plot data
  oplot, x, imsum, thi=thi
  oplot, x, c, thi=0.2*thi, co=c24(2)

;plot window
  plot, xr, yr2, /nodata $
      , xr=xr, /xsty, yr=yr2, /ysty $
      , xthi=thi, ythi=thi $
      , xtit=xtit, yminor=5 $
      , ytit='Normalized Flat' $
      , chars=chars $
      , xmarg=[10.5,2.2], ymarg=[0.0,0.5]

;annotation coordinates
  xcr = !x.crange
  dxcr = xcr[1] - xcr[0]
  ycr = !y.crange
  dycr = ycr[1] - ycr[0]

;plot data
  oplot, x, y, thi=thi
  oplot, xcr, [1,1], thi=thi, lin=2, col=c24(2)

;indicate anomaly spacing
  arrow, /data, 1950, 1.02, 1700, 1.02, hsize=-0.4, thi=thi, col=c24(4)
  arrow, /data, 2430, 1.02, 2680, 1.02, hsize=-0.4, thi=thi, col=c24(4)
  xyouts, 0.5*(1700+2680), 1.0203, siz=0.8*chars, col=c24(4) $
        , align=0.5, '1024'
  xyouts, 0.5*(1700+2680), 1.0167, siz=0.8*chars, col=c24(4) $
        , align=0.5, 'Pixels'

;indicate 15 jumps per 1024 pixels
  iw = where(xjmp ge 2730 and xjmp le 3690, nw)
  for i=0, nw-1 do begin
    oplot, xjmp[iw[i]]+[0,0], [1.009,1.0115], thi=thi, col=c24(4)
  endfor
  xyouts, total(xjmp[iw])/nw, 1.013, siz=0.8*chars, align=0.55 $
        , ' 15 per 1024', col=c24(4)

;plot cleanup
  !p.multi = 0
  !y.omargin = 0

;close postscript file.
  device, /close
  set_plot, 'x'
  !p.font = -1
  print, psfile + ' closed'

end
