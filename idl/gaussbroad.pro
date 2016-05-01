function gaussbroad,w,s,hwhm
;Smooths a spectrum by convolution with a gaussian of specified hwhm.
; w (input vector) wavelength scale of spectrum to be smoothed
; s (input vector) spectrum to be smoothed
; hwhm (input scalar) half width at half maximum of smoothing gaussian.
;Returns a vector containing the gaussian-smoothed spectrum.
;Edit History:
;  -Dec-90 GB,GM Rewrote with fourier convolution algorithm.
;  -Jul-91 AL	Translated from ANA to IDL.
;22-Sep-91 JAV	Relaxed constant dispersion check; vectorized, 50% faster.
;05-Jul-92 JAV	Converted to function, handle nonpositive hwhm.

if n_params() lt 3 then begin
  print,'syntax: snew=gaussbroad(w,sold,hwhm)'
  retall
endif

;Warn user if hwhm is negative.
  if hwhm lt 0.0 then $
    message,/info,'Warning! Forcing negative smoothing width to zero.'

;Return input argument if half-width is nonpositive.
  if hwhm le 0.0 then return,s			;true: no broadening

;Calculate (uniform) dispersion.
  nw = n_elements(w)				;# points in spectrum
  dw = (w(nw-1) - w(0)) / (nw-1)		;wavelength change per pixel

;Make smoothing gaussian; extend to 4 sigma.
;Note: 4.0 / sqrt(2.0*alog(2.0)) = 3.3972872 and sqrt(alog(2.0))=0.83255461
;  sqrt(alog(2.0)/pi)=0.46971864 (*1.0000632 to correct for >4 sigma wings)
  if(hwhm gt 5*(w(nw-1) - w(0))) then return, replicate(total(s)/nw, nw)
  nhalf = long(3.3972872d0 * hwhm/dw)		;# points in half gaussian
  ng = 2L * nhalf + 1L				;# points in gaussian (odd!)
  wg = dW * (findgen(nG) - (ng-1L)/2.d0)	;wavelength scale of gaussian
  xg = (0.83255461d0 / hwhm) * wg 		;convenient absisca
  gpro = (0.46974832d0 * dw / hwhm) * exp(-xg*xg);unit area gaussian w/ FWHM
  gpro=gpro/total(gpro)

;Pad spectrum ends to minimize impact of Fourier ringing.
  npad = nhalf + 2				;# pad pixels on each end
  spad = [replicate(s(0),npad),s,replicate(s(nw-1),npad)]

;Convolve and trim.
  sout = convol(spad,gpro)			;convolve with gaussian
  sout = sout(npad:npad+nw-1)			;trim to original data/length
  return,sout					;return broadened spectrum.

End
