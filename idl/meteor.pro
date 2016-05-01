pro meteor, obsid, meteor, zmed, mar, y7at2k $
          , plot=plot, print=print, trace=trace, errmsg=errmsg

;syntax
  if n_params() lt 1 then begin
    print, 'syntax: meteor, obsid' $
         + ' [,meteor ,zmed ,mar ,y7at2k ,/plot ,/trace]'
    print, "  e.g.: meteor, 'rj76.106'"
    return
  endif

;default values of optional parameters
  if n_elements(plot) eq 0 then plot = 0
  trc = keyword_set(trace)

;clear return variables.
  meteor = 0
  zmed = 0
  mar = 0
  y7at2k = 0
  errmsg = ''

;internal parameters
  xnomb = [ 3072, 2977, 2879, 2780, 2678, 2574, 2469, 2361 $
          , 2250, 2136, 2020, 1899, 1775, 1648, 1518, 1384 ]
  xnomt = [ 3021, 2924, 2826, 2724, 2613, 2525, 2420, 2313 $
          , 2202, 2087, 1970, 1848, 1723, 1594, 1465, 1331 ]
  nnomb = n_elements(xnomb)
  if n_elements(xnomt) ne nnomb then begin
    message, 'nominal crossings above and below disagree'
  endif

;location of raw images
  rawdir = getenv('GRAND_RAWDIR') + path_sep()

;parse obsid
  words = strsplit(obsid, '.', count=nword, /extract)
  if nword ne 2 then message, 'obsid must contain one period'
  if strmid(words[0], 0, 1) ne 'r' then message, "obsid must begin with 'r'"
  run = strmid(words[0], 1)
  obsnum = long(words[1])

;read data form fits file
  file = rawdir + run + string(obsnum, form='(i4.4)') + '.fits'
  if ~file_test(file) then begin
    message, obsnum + 'has no raw image: ' + file
  endif
  im = mrdfits(file, 2, hd2, /dscale, /silent)
  im = rotate(temporary(im), 3)
  hd0 = headfits(file)

;subtract bias.
  im -= median(im[*,5:13])
  sz = size(im)
  nx = sz[1]
  ny = sz[2]

;check CCD binning
  if ny ne 713 then begin
    errmsg = ' CCD image is not binned'
    print, obsid + ' ' + errmsg
    return
  endif

;trim image
  im = im[0:4020,24:706]
  sz = size(im)
  nx = sz[1]
  ny = sz[2]

;skipping nonlinearity correction.
; im = nonlinear(im)

;median filter image, since all we are doing is locating orders and the meteor
  im = median(im, 3)

;find orders
  swid = 32
  m_fords, im, swid, orc, trace=trc
  nord = n_elements(orc[0,*])
  if nord ne 16 then message, 'expected 16 orders, got ' + strtrim(nord, 2)
  if nord ne nnomb then message, 'need nominal meteor location for each order'
  y7at2k = poly(2000, orc[*,7])

;extend orders
  ncoef = n_elements(orc[*,0])
  xorc = dblarr(ncoef, nord+2)
  xorc[*,0] = 2 * orc[*,0] - orc[*,1]
  xorc[*,1:nord] = orc
  xorc[*,nord+1] = 2 * orc[*,nord-1] - orc[*,nord-2]

;extract background vector below each order
  hwsp = 6
  hwbg = 2
  sb = dblarr(nx,nord)
  x = dindgen(nx)
  ybb = dblarr(nx,nord)
  ybt = dblarr(nx,nord)
  for iord=0, nord-1 do begin
    m_getarc, im, xorc, iord+1, 2*hwbg, arc, dypix=-hwsp-hwbg, yb=yb, yt=yt
    sb[*,iord] = median(arc, 3)
    ybb[*,iord] = yb
    ybt[*,iord] = yt
  endfor

;extract background vector above each order
  st = dblarr(nx,nord)
  ytb = dblarr(nx,nord)
  ytt = dblarr(nx,nord)
  for iord=0, nord-1 do begin
    m_getarc, im, xorc, iord+1, 2*hwbg, arc, dypix=hwsp+hwbg, yb=yb, yt=yt
    st[*,iord] = median(arc, 3)
    ytb[*,iord] = yb
    ytt[*,iord] = yt
  endfor

;calculate below order minus above order difference spectrum
  hwg = 5
  diff = sb - st
  diff -= median(diff)
  for iord=0, nord-1 do diff[*,iord] = gaussbroad(x, diff[*,iord], hwg)

;determine approximate offset of meteor using central order
  zthr = (10 * median(abs(diff))) < 15
  hwtols = [50, 100, 150, 200]
  ihwtol = 0
  repeat begin
    if ihwtol ge n_elements(hwtols) then begin
      message, 'no meteor in ' + obsid + ' even in largest search window'
    endif
    hwtol = hwtols[ihwtol]
    iord = nord / 2
    ibb = round(xnomb[iord] - hwtol)
    ieb = round(xnomb[iord] + hwtol)
    zb = max(diff[ibb:ieb,iord], imax)
    ibt = round(xnomt[iord] - hwtol)
    iet = round(xnomt[iord] + hwtol)
    zt = -min(diff[ibt:iet,iord], imin)
    ihwtol++
  endrep until zb gt zthr or zt gt zthr
  xshiftb = ibb + imax - xnomb[iord]
  xshiftt = ibt + imin - xnomt[iord]
  xshift = zb gt zt ? xshiftb : xshiftt
  if trc then begin
    print, 'meteor shift: xb=' + strtrim(round(xshiftb), 2) $
         + ' pix, zb=' + strtrim(round(zb), 2)
    print, 'meteor shift: xt=' + strtrim(round(xshiftt), 2) $
         + ' pix, zt=' + strtrim(round(zt), 2)
  endif

;set actual meteor location BELOW to maximum in search window
  hwtol = 30
  xmetb = fltarr(nord)
  ymetb = fltarr(nord)
  zmetb = fltarr(nord)
  for iord=0, nord-1 do begin
    ib = round(xnomb[iord] - hwtol) + xshift
    ie = round(xnomb[iord] + hwtol) + xshift
    zmetb[iord] = max(diff[ib:ie,iord], imax)
    ymetb[iord] = 0.5 * (ybb[ib+imax,iord] + ybt[ib+imax,iord])
    xmetb[iord] = ib + imax
  endfor

;set actual meteor location ABOVE to minimum in search window
  hwtol = 30
  xmett = fltarr(nord)
  ymett = fltarr(nord)
  zmett = fltarr(nord)
  for iord=0, nord-1 do begin
    ib = round(xnomt[iord] - hwtol) + xshift
    ie = round(xnomt[iord] + hwtol) + xshift
    zmett[iord] = min(diff[ib:ie,iord], imin)
    ymett[iord] = 0.5 * (ytb[ib+imin,iord] + ytt[ib+imin,iord])
    xmett[iord] = ib + imin
  endfor

;combine measurements. discard weak peaks
  xmet = [xmetb, xmett]
  ymet = [ymetb, ymett]
  zmet = [zmetb, zmett]
  igd = where(abs(zmet) gt zthr, ngd)
  if ngd lt 5 then message, obsid + ' has too few points'
  xmet = xmet[igd]
  ymet = ymet[igd]
  zmet = zmet[igd]
  isort = sort(xmet)
  xmet = xmet[isort]
  ymet = ymet[isort]
  zmet = zmet[isort]
  zmed = median(abs(zmet))

;robust quadratic fit to meteor location
  cmet = poly_fit(xmet, ymet, 1, /double)
  resid = ymet - poly(xmet, cmet)
  thresh = 5 * median(abs(resid))
  igd = where(abs(resid) lt thresh, ngd)
  cmet = poly_fit(xmet[igd], ymet[igd], 2, /double)
  resid = ymet - poly(xmet, cmet)
  thresh = 5 * median(abs(resid))
  igd = where(abs(resid) lt thresh, ngd)
  cmet = poly_fit(xmet[igd], ymet[igd], 2, yfit=fmet, /double)
  resid = ymet - poly(xmet, cmet)
  mar = median(abs(resid))

;calculate where meteor crosses each order
  meteor = intarr(nord)
  for iord=0, nord-1 do begin
    coef = orc[*,iord]
    coef[0:2] = coef[0:2] - cmet
    roots = fz_roots(coef, /double)
    iwhr = where(real_part(roots) gt 0 $
             and real_part(roots) lt nx-1 $
             and imaginary(roots) eq 0, nwhr)
    if nwhr ne 1 then message, strtrim(nwhr, 2) + ' meteor intersections'
    meteor[iord] = round(real_part(roots[iwhr]))
  endfor

;;
;; DIAGNOSTICS
;;

;print results
  if keyword_set(print) then begin
    fmt = '(a-10,f7.2,16i5,f6.1,f6.3)'
    print, form=fmt, obsid, y7at2k, meteor, zmed, mar
  endif

;calculate mode
  if plot ge 1 then begin
    h = histogram(im, min=0)
    junk = max(h, mode)

;display image
    dz = 20
    window, 0, xsize=1000, ysize=800
    display, im, min=mode-dz, max=mode+dz

;overplot orders ridgelines and background boundaries
    cosp = c24(4)
    cob = c24(3)
    cot = c24(2)
    oplot, x, poly(x, xorc[*,0]), co=cosp
    for iord=0, nord-1 do begin
      oplot, x, poly(x, orc[*,iord]), co=cosp
      oplot, x, ybb[*,iord], co=cob
      oplot, x, ybt[*,iord], co=cob
      oplot, x, ytb[*,iord], co=cot
      oplot, x, ytt[*,iord], co=cot
    endfor
    oplot, x, poly(x, xorc[*,nord+1]), co=cosp

;overplot meteor location
    psym, 'circle'
    oplot, xmet, ymet, ps=8, co=c24(11)
    psym, 'circle', /fill
    oplot, x, poly(x, cmet), lin=2, co=c24(11)
    if plot gt 1 then if get_kbrd(1) eq 'q' then retall
  endif

;plot extracted background spectra
  if plot ge 2 then begin
    ib = 1100
    ie = 3300
    yr = minmax(diff[ib:ie,*])
    plot, x[ib:ie], diff[ib:ie,0], /nodata $
        , /xsty, yr=yr, /ysty
    for iord=0, nord-1 do oplot, x[ib:ie], diff[ib:ie,iord]
    oplot, !x.crange, [0,0], col=c24(4)

;overplot meteor peaks and threshold
    oplot, !x.crange, [0,0]+zthr, lin=2, co=c24(2)
    oplot, !x.crange, [0,0]-zthr, lin=2, co=c24(2)
    ipos = where(zmet gt 0, npos, comp=ineg, ncomp=nneg)
    if npos gt 0 then oplot, xmet[ipos], zmet[ipos], ps=-8, co=c24(5)
    if nneg gt 0 then oplot, xmet[ineg], zmet[ineg], ps=-8, co=c24(5)
  endif

end
