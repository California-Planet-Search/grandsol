pro jump, xjmp, print=print

  if n_params() lt 1 then begin
    print, 'syntax: jump, xjmp [,/print]'
    return
  endif

;Locations of features in blue and green detectors (red is different).
  xjmp = [   39,  107,  175,  244,  312,  380,  448,  517,  585,  653 $
         ,  722,  790,  858,  926,  995, 1063, 1131, 1199, 1268, 1336 $
         , 1404, 1472, 1541, 1609, 1677, 1746, 1814, 1882, 1950, 2019 $
         , 2087, 2155, 2223, 2292, 2360, 2428, 2496, 2565, 2633, 2701 $
         , 2770, 2838, 2906, 2974, 3043, 3111, 3179, 3247, 3316, 3384 $
         , 3452, 3520, 3589, 3657, 3725, 3794, 3862, 3930, 3998 ]
  njmp= n_elements(xjmp)
  dxjmp = xjmp[1:njmp-1] - xjmp[0:njmp-2]

;Pattern repeats every 1024 pixels (15 intervals).
  if keyword_set(print) then begin
    print, form='(15i5)', xjmp
    print, form='(15i5)', dxjmp
  endif

end
