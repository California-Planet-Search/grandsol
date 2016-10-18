pro simobs_one, star, iod, elsf, ech, velsum, ston, file $
              , verbose=verbose, seed=seed, dump=dump
;
;Inputs:
; star (string) - identifier for an intrinsic stellar spectrum, e.g. 'sun'
; iod (string) - identifier for an iodine cell spectrum, e.g. 'keck-nso-52.5'
; elsf (string) - identifier for an effective LSF, e.g. 'jay2011', 'laser_y'
; ech (string) - identifier for echelle spectrograph, e.g. 'hires'
; velsum (scalar) - [m/s] relativistically summed velocities
; ston (scalar) - continuum signal to noise [0 --> noiseless]
; file (string) - name of output file
; [seed=] - seed for random number generator (gives repeatable noise)
; [dump=] - 0 or undefined: only write simulated observations to disk
;           1: also write pixel level diagnostics for each order to disk
;           2: also write supsampled pixel diagnostics to disk
;
;Output:
; Output file is created in the current working directory.
;
;History:
; 2012-Oct-23 Valenti Adapted from make_output.pro

;Syntax
  if n_params() lt 7 then begin
    print, 'syntax: simobs_one, star, iod, elsf, ech, velsum, noise, file'
    print, '                 [, seed=, dump=]'
    print, "  e.g.: simobs_one, 'sun', 'keck-nso-52.5', 'jay2012', 'hires'" $
         + ", 27345.6, 0, 'sun.023'"
    return
  endif

;Handle optional keyword.
  if keyword_set(seed) then begin
    seed_entry = seed
    seed_exit = seed
  endif

;Internal parameters.
  osamp = 10					;grid points per pixel

;Internal constants.
  c = 2.99792458d8                                      ;m/s

;Location of reference files.
  dir_base = getenv('GRAND_IREFDIR') + path_sep()

;;
;; INPUTS
;;

;Radial velocity factor
  beta = velsum / c
  zplus1 = sqrt((1+beta)/(1-beta))

;
; Echelle Spectrum Characteristics
;

;Restore wavelength array for echelle spectrograph
  dir_ech = dir_base + 'ech' + path_sep()
  file_wave = dir_ech + ech + '.wave'
  if keyword_set(verbose) then print, 'restoring ' + file_wave
  restore, file_wave
  sz = size(wave)
  npix = sz[1]
  nord = sz[2]

;Restore continuum normalization for echelle spectrograph
  file_norm = dir_ech + ech + '.norm'
  if keyword_set(verbose) then print, 'restoring ' + file_norm
  restore, file_norm
  sz = size(norm)
  if sz[1] ne npix then message, 'wave and norm have different pixels'
  if sz[2] ne nord then message, 'wave and norm have different orders'

;
; Intrinsic Stellar Spectrum (ISS)
;

;Solar atlas.
  if star eq 'sun' then begin
    file_star = 'Kurucz et al. (1984) Solar Flux Atlas from 296 to 1300 nm'
    if keyword_set(verbose) then print, 'reading NSO solar atlas'
    iwave = where(wave gt 0)
    wave_min = min(wave[iwave], max=wave_max)
    rdnso, wstar, sstar, wave_min-10, wave_max+10

;Spectra stored in IDL save files.
  endif else begin
    dir_star = dir_base + 'star' + path_sep()
    file_star = dir_star + star + '.star'
    if keyword_set(verbose) then print, 'restoring ' + file_star
    restore, file_star
  endelse

;Wavelength interval.
  wstar_min = min(wstar, max=wstar_max)

;
; Line Spread Function (LSF)
;

;Restore effective LSF.
  dir_elsf = dir_base + 'elsf' + path_sep()
  file_elsf = dir_elsf + elsf + '.elsf'
  if keyword_set(verbose) then print, 'restoring ' + file_elsf
  restore, file_elsf
  xlsf_in = xlsf
  ylsf_in = ylsf

;Interpolate LSF onto properly sampled grid.
  xmin = min(xlsf, max=xmax)
  nlsf = round((xmax - xmin) * osamp) + 1
  xlsf = xmin + dindgen(nlsf) / double(osamp)
  splc = spl_init(xlsf_in, ylsf_in, /double)
  ylsf = spl_interp(xlsf_in, ylsf_in, splc, xlsf, /double)

;Renormalize LSF
  ylsf /= total(ylsf, /double)

;Deterimine half-widths on both sides
  iz = where(xlsf eq 0, nz)
  if nz eq 0 then message, 'xlsf must have a point at zero'
  hwlo = iz[0]
  hwhi = n_elements(xlsf) - hwlo - 1

;Create pixel indexes for computational grid
  nx = hwlo + ((npix - 1) * osamp + 1) + hwhi
  x = (dindgen(nx) - hwlo) / double(osamp)

;;
;; CONSTRUCT SIMULATED OBSERVATION FOR EACH ORDER
;;

;Initialize output spectrum
  spec = fltarr(npix, nord)

;Generate a random shoft for all orders.
  pixsh = 2d0 * (randomu(seed_exit, 1))[0] - 1d0

;Loop through echelle orders.
  for iord=0, nord-1 do begin

;
; Generate Oversampled Wavelengths for Extended Order
;

;Wavelengths for current order.
    word = wave[*,iord]

;Skip echelle orders with no wavelengths
    if max(word) eq 0 then begin
      if keyword_set(verbose) then begin
        print, 'skipping order ' + strtrim(iord+1, 2) + ' (missing wavelengths)'
      endif
      continue
    endif

;Fit wavelength solution with polynomial. 'pix' starts at index 1.
    pixb = 20
    pixe = 4000
    pix = pixb + dindgen(pixe - pixb + 1)
    wcoef = poly_fit(pix, word[pixb:pixe], 4, yfit=yfit, /double)
    if max(abs(word[pixb:pixe]-yfit)) gt 1e-6 then begin
      message, 'wavelengths are not a 4th order polynomial'
    endif

;Extend wavelengths beyond order ends and introduce random pixel shift.
    w = poly(x+pixsh, wcoef)
    wmin = min(w, max=wmax)
    wminz = wmin / zplus1
    wmaxz = wmax / zplus1

;
; Extract Stellar Spectrum for Extended Order
;

;Check whether stellar spectrum spans extended order.
    if wminz lt wstar_min or wmaxz gt wstar_max then begin
      if keyword_set(verbose) then begin
        print, 'skipping order ' + strtrim(iord+1, 2) $
             + ' (missing stellar spectrum)'
      endif
      continue
    endif

;interpolate onto computational grid
    yiss = interpol(sstar, wstar*zplus1, w)

;
; Iodine Transmission for Extended Order
;

;Parse the iodine specification.
    words = strsplit(iod, '-', /extract, count=count)
    if count ne 3 then message, 'iod must have two dashes: ' + iod
    iodcell = words[0]
    iodlab = words[1]
    iodtemp = float(words[2])

;Read iodine transmission spectrum
    iod_param = 'cell=' + iodcell + ', lab=' + iodlab $
              + ', temp=' + string(iodtemp, form='(f4.1)') + ' C'
    read_iodine, iodcell, iodlab, iodtemp, wmin, wmax, wiod, tiod, nadd=1

;Interpolate onto computational grid
    splc = spl_init(wiod, tiod, /double)
    yiod = spl_interp(wiod, tiod, splc, w, /double)

;
; Composite Spectrum for Extended Order
;

    if n_elements(yiss) ne n_elements(yiod) then begin
      message, 'yiss and yiod have different lengths'
    endif
    ycom = yiss * yiod

;
; Convolve Extended Order with LSF.
;

    y = convol(ycom, ylsf)

;
; Sample at Original Pixels.
;

    jsim = hwlo + osamp * lindgen(npix)
    ssim = y[jsim]

;
; Add noise
;

   if keyword_set(ston) and max(ssim) gt 0 then begin
     denom = sqrt(ston^2.0 * ssim / max(ssim))
     sigma = 1.0 / denom
     noise = sigma * randomn(seed_exit, npix)
     ssim += noise
   endif else begin
     noise = replicate(0.0, npix)
   endelse

;
; Normalize Continuum of Current Order. Save into Full Array.
;

    spec[*,iord] = ssim * norm[*,iord]

;Write auxilliary pixel data for current order to file.
    file_aux = file + '.' + string(iord+1, '(i2.2)')
    if keyword_set(dump) && dump ge 1 then begin
      cmd = "simobs_one,'" + star + "','" + iod + "','" + elsf + "','" $
          + ech + "'," + strtrim(string(velsum,'(f19.2)'),2) + "," $
          + strtrim(ston,2) + ",'" + file + "'"
;     if keyword_set(seed_entry) then cmd += ",seed=" + strtrim(seed_entry,2)
      if keyword_set(dump) then cmd += ",dump=" + strtrim(dump,2)
      rx = (dindgen(npix)-2010) / 1990d0
      coef = reform(poly_fit(rx, w[jsim], 4, yfit=wfit, /doub))
      if max(abs(1 - w[jsim]/wfit)) gt 1e-10 then begin
        message, 'wavelengths are not quadratic'
      endif
      openw, u, file_aux, /get
      printf, u, '# ' + cmd
      printf, u, '#'
      printf, u, '# oversampling: ' + strtrim(osamp,2) + ' points/pixel'
      printf, u, '# input wavelengths: ' + file_wave
      printf, u, '# pixel shift: ' + strtrim(pixsh,2) + ' pixels'
      printf, u, '# iodine spectrum: ' + iod_param
      printf, u, '# stelar spectrum: ' + file_star
      printf, u, '# stellar velocity shift: ' + strtrim(velsum,2) + ' m/s'
      printf, u, '# z+1: ' + strtrim(string(zplus1, form='(f20.11)'),2)
      printf, u, '# effective LSF: ' + file_elsf
      printf, u, '# signal/noise: ' + strtrim(ston,2)
      printf, u, '# order normalization: ' + file_norm
      printf, u, '#'
      printf, u, '# waveleng = polyval(rx,coef)'
      printf, u, '# rx = np.arange(-2010,2011) / 1990.0
      printf, u, '# coef = [' $
               + strjoin(strtrim(string(coef,'(g20.10)'),2),',') + ']'
      printf, u, '# obs_spec = (temp*iod + noise) * norm_vec'
      printf, u, '#'
      printf, u, '# waveleng temp*iod     noise norm_vec obs_spec'
      printf, u, '# -------- -------- --------- -------- --------'
      for ipix=0, npix-1 do begin
        j = jsim[ipix]
        printf, u, form='(f10.5,f9.6,f10.6,2f9.2)' $
              , w[j], ssim[ipix]-noise[ipix], noise[ipix] $
              , norm[ipix,iord], spec[ipix,iord]
      endfor
      free_lun, u
    endif else begin
      if file_test(file_aux) then print,'WARNING - ' + file_aux + ' is stale'
    endelse

;Write auxilliary supersampled data for current order to file.
    file_aux = file + '.' + string(iord+1, '(i2.2)') + 'x'
    if keyword_set(dump) && dump ge 2 then begin
      openw, u, file_aux, /get
      printf, u, '# convolved = np.convolve(stellar*iodine, elsf)'
      printf, u, '#'
      printf, u, '# pixel wavelength stellar   iodine  star*iod convolved'
      printf, u, '# ----- ---------- -------- -------- -------- --------'
      for ix=0, nx-1 do begin
        printf, u, form='(f7.2,f11.5,4f9.6)' $
              , x[ix], w[ix], yiss[ix], yiod[ix], ycom[ix], y[ix]
      endfor
      free_lun, u
    endif else begin
      if file_test(file_aux) then print,'WARNING - ' + file_aux + ' is stale'
    endelse

;Done looping through echelle orders
  endfor

;;
;; OUTPUT
;;

;
; Write Simulated Observation to .dsk File
;

  if keyword_set(verbose) then print, 'writing ' + file
  wdsk_simobs, spec, file, /new, /swap

;
; Construct Header
;

  head = $
    [ 'star: ' + file_star $
    , 'iodine: ' + iod_param $
    , 'elsf: ' + file_elsf $
    , 'echwave: ' + file_wave $
    , 'echnorm: ' + file_norm $
    , 'velsum: ' + strtrim(velsum, 2) $
    , 'ston: ' + strtrim(ston, 2) $
    , 'z+1: ' + strtrim(string(zplus1, form='(f20.11)'), 2) $
    , 'pixsh: ' + strtrim(pixsh, 2) $
    , 'lsfhwlo: ' + strtrim(hwlo, 2) $
    , 'lsfhwhi: ' + strtrim(hwhi, 2) $
    , 'osamp: ' + strtrim(osamp, 2) $
    , 'idlpro: ' + 'simobs_one.pro' ]

;
; Write Header to .dsk File
;

  wdsk_simobs, head, file, /swap

;
; Write seeds to .dsk File
;

  if keyword_set(seed_entry) then begin
    wdsk_simobs, seed_entry, file, /swap
    wdsk_simobs, seed_exit, file, /swap

;
; Update Seed Passed by Caller
;

    seed = seed_exit
  endif

end
