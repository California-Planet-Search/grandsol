pro simobs_set, star, iod, elsf, ech, bcmax, rvk, vsys, ston, nobs $
              , seed=seed
;
;Inputs:
; star (string) - identifier for an intrinsic stellar spectrum, e.g. 'sun'
; iod (string) - identifier for an iodine cell spectrum, e.g. 'keck-nso-52.5'
; elsf (string) - identifier for an effective LSF, e.g. 'jay2011', 'laser_y'
; ech (string) - identifier for echelle spectrograph, e.g. 'hires'
; bcmax (scalar) - [m/s] maximum barycentric correction
; rvk (scalar) - [m/s] semi-amplitude of stellar velocity perturbations
; vsys (scalar) - [m/s] radial velocity of the system
; ston (scalar) - cotinuum signal to noise ratio [value of 0 mean no noise]
; nobs (integer) - number of simulated observations to create
; [seed=] - seed for random number generator (gives repeatabable noise)
;
;Output:
; Data files are created in the current working directory.
;   - Simulated observations named <star>.###
;   - Observation list file named obslist
;
;History:
; 2012-Oct-23 Valenti Adapted from make_set.pro

;Syntax
  if n_params() lt 9 then begin
    print, 'syntax: simobs_set, star, iod, elsf, ech, bcmax, rvk, vsys, ston, nobs'
    print, '                 [, seed=]'
    print, "  e.g.: simobs_set, 'sun', 'keck-nso-52.5', 'jay2011', 'hires'" $
         + ', 3e4, 3, 0, 200, 100, seed=200'
    print, "  e.g.: simobs_set, 'sun', 'keck-nso-52.5', 'jay2012', 'hires'" $
         + ', 3e4, 300, 0, 0, 100, seed=200'
    return
  endif

;Internal constants
  c = 2.99792458d8					;m/s

;Handle optional keyword.
  if keyword_set(seed) then begin
    seed_entry = seed
    seed_exit = seed
  endif

;;
;; Generate Velocities
;;

;Generate uniformly distributed random stellar velocity perturbations
  velorb = [0.0, rvk * (2 * randomu(seed_exit, nobs-1, /double) - 1)]

;Generate uniformly distributed random barycentric corrections
  barycorr = [0.0, bcmax * (2 * randomu(seed_exit, nobs-1, /double) - 1)]

;Calculate beta = v/c
  velear = -barycorr					;m/s
  betsys = vsys   / c
  betorb = velorb / c
  betear = velear / c

;Calculate relativistic velocity sum.
  betsum = (betsys+betorb+betear + betsys*betorb*betear) $
         / (1 + betsys*betorb + betsys*betear + betorb*betear)
  velsum = betsum * c

;Sort velocities into ascending order.
  isort = sort(velsum)
  velsum = velsum[isort]
  velorb = velorb[isort]
  barycorr = barycorr[isort]

;;
;; Create Observations
;;

;Create new observation list file.
  openw, unit, 'obslist', /get_lun

;Write header.
  cd, curr=cwd
  printf, unit, 'VSYST = ' + fmtnum(vsys, '(f20.10)') + ' m/s'
  printf, unit, 'RJDIR = "' + cwd + '/"'

;Loop through observations
  for iobs=0, nobs-1 do begin
    file = star + '.' + string(iobs+1, form='(i3.3)')

;Generate observations and associated files
    simobs_one, star, iod, elsf, ech, velsum[iobs], ston, file $
              , seed=seed_exit

;Update observation list
    printf, unit, form='(i3.3,2x,a,i3,2f14.5)' $
          , iobs+1, file, 0, barycorr[iobs], velorb[iobs]

;Done looping through observations
  endfor

;Close observation list
  free_lun, unit

;Report final seed value.
  if keyword_set(seed) then seed = seed_exit

;Write simulation parameters.
  openw, unit, 'params.txt', /get_lun
  printf, unit, "star = '" + star + "'"
  printf, unit, "iod = '"  + iod + "'"
  printf, unit, "elsf = '" + elsf + "'"
  printf, unit, "ech = '"  + ech + "'"
  printf, unit, "bcmax = " + strtrim(string(bcmax, form='(f9.2)'), 2) + ' m/s'
  printf, unit, "rvk = " + strtrim(string(rvk, form='(f9.2)'), 2) + ' m/s'
  printf, unit, "vsys = " + strtrim(string(rvk, form='(f9.2)'), 2) + ' m/s'
  printf, unit, "ston = " + strtrim(ston, 2)
  printf, unit, "nobs = " + strtrim(nobs, 2)
  if keyword_set(seed_entry) then begin
    printf, unit, "seed_entry = " + strtrim(seed_entry, 2)
  endif
  free_lun, unit

end
