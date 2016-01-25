pro plot_simobs, ord, pixb, pixe

;syntax
  if n_params() lt 3 then begin
    print, 'syntax: plot_simobs, ord, pixb, pixe'
    print, '  e.g.: plot_simobs, 6, 2450, 2550'
    print, '  e.g.: plot_simobs, 1, 1250, 1400'
    return
  endif

;build requested pixel scale
  npix = pixe - pixb + 1
  pix = pixb + findgen(npix)

;simulated observations in current working directory
  files = file_search('*.[0-9][0-9][0-9]', count=nfile)

;restore wavelength solution used to construct first simulated observation
  rdsk_struct, files[0], dsk
  head = dsk.r2.data
  ihead = where(strmid(head, 0, 8) eq 'echwave:', count)
  if count eq 0 then begin
    print, 'no echwave record in header of ' + files[0]
    return
  endif
  file_wave = strtrim(strmid(head[ihead[0]], 8), 2)
  restore, file_wave
  wave = wave[pixb:pixe,ord]
  wmin = min(wave, max=wmax)
  if wmax eq 0 then begin
    print, 'no data for order ' + strtrim(ord, 2)
    return
  endif
  dpix = pix - 0.5 * (pixb + pixe)
  wcoef = poly_fit(dpix, wave, 4, yfit=yfit, /double)
  if max(abs(wave-yfit)) gt 1e-6 then message, 'error fitting wavelengths'

;read each simulated observation and associated value of z+1 shift
  specs = fltarr(npix, nfile)
  zplus1 = dblarr(nfile)
  pixsh = fltarr(nfile)
  for ifile=0, nfile-1 do begin
    file = files[ifile]
    rdsk_struct, file, dsk
    specs[*,ifile] = dsk.r1.data[pixb:pixe,ord]

    head = dsk.r2.data
    ihead = where(strmid(head, 0, 4) eq 'z+1:', count)
    if count eq 0 then message, 'no z+1 record in header of ' + file
    zplus1[ifile] = double(strmid(head[ihead[0]], 4))
    ihead = where(strmid(head, 0, 6) eq 'pixsh:', count)
    if count eq 0 then message, 'no pixsh record in header of ' + file
    pixsh[ifile] = float(strmid(head[ihead[0]], 6))
  endfor

;read iodine spectrum
  ihead = where(strmid(head, 0, 7) eq 'iodine:', count)
  if count eq 0 then message, 'no iodine record in header of ' + file
  iodine = strtrim(strmid(head[ihead[0]], 7), 2)
  words = strtrim(strsplit(iodine, ',', /extract, count=nword), 2)
  if nword ne 3 then message, 'iodine info has ' + strtrim(nword, 2) + ' words'
  cell = strmid(words[0], 5)
  lab = strmid(words[1], 4)
  temp = float(strmid(words[2], 5))
  read_iodine, 'keck', lab, temp, wmin, wmax, wiod, tiod, nadd=10
  siod = gaussbroad(wiod, tiod, 0.25*(wmin+wmax)/80000)

;read solar spectrum
  rdnso, wsun, ssun, wmin, wmax

;plot result
  !p.multi = [0, 1, 2]
  for loop=0, 1 do begin
    xtit = loop ? 'Stellar Wavelengths' : 'Observatory Wavelengths'
    plot, wave, specs[*,0], /nodata $
        , xsty=3, ytickf='(i6)' $
        , xtit=xtit, ytit='Signal (ADU)', chars=1.5

    for ifile=0, nfile-1 do begin
      x = poly(dpix + pixsh[ifile], wcoef)
      if loop eq 1 then x /= zplus1[ifile]
      line = max(ifile eq [0, round(nfile/2), nfile-1])
      psym = 3 * (1 - line)
      color = c24(line ? 2 + round(2 * float(ifile) / nfile) : 1)
      oplot, x, specs[*,ifile], psym=psym, thi=2, col=color

    endfor

;overplot iodine spectrum
    if loop eq 0 then begin
      scale = max(specs[*,0]) / max(siod)
      oplot, wiod, scale*siod, thi=4, col=!p.background
      oplot, wiod, scale*siod, thi=2, col=c24(6)
    endif

;overplot solar spectrum
    if loop eq 1 then begin
      scale = max(specs[*,0]) / max(ssun)
      oplot, wsun, scale*ssun, thi=4, col=!p.background
      oplot, wsun, scale*ssun, thi=2, col=c24(5)
    endif
  endfor
  !p.multi = 0

;Diagnostics
  c = 2.9979246d8			;m/s
  disp = (wave[npix-1]-wave[0]) / (npix - 1)
  print, strtrim(string(min((zplus1-1)*c), form='(f19.2)'), 2) $
       + ' < z*c < ' $
       + strtrim(string(max((zplus1-1)*c), form='(f19.2)'), 2) $
       + ' m/s'
  print, strtrim(string(min(pixsh), form='(f19.2)'), 2) $
       + ' < pixsh < ' $
       + strtrim(string(max(pixsh), form='(f19.2)'), 2) $
       + ' pixels'
  print, 'dispersion = ' $
       + strtrim(string(disp, form='(f19.5)'), 2) $
       + ' A/pixel = ' $
       + strtrim(string(c*disp/mean(wave)/1e3, form='(f19.2)'), 2) $
       + ' km/s/pixel'

end
