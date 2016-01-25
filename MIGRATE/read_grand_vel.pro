pro read_grand_vel, file, vel
;
;Purpose:
; Read velocity data from a .vel file produced by grand.
;
;Input:
; file (string) name of .vel file containing velocity data
;
;Output:
; vel (structure) velocity data from a .vel file
;
;Notes:
;
;  vbarn = c * ((1+zbarn)^2    - 1) / ((1+zbarn)^2    + 1)
;
;  V = c * ((1+Z)^2 - 1) / ((1+Z)^2 + 1)
;
;  uV^2 = (dV/dZ)^2 * uZ^2
;
;         ((1+Z)^2 + 1) * (2*(1+Z)) - ((1+Z)^2 - 1) * (2*(1+Z))
;       = ----------------------------------------------------- * uZ^2
;                             ((1+Z)^2 + 1)^2
;
;         2*(1+Z) * ((1+Z)^2 + 1 - (1+Z)^2 + 1)
;       = ------------------------------------- * uZ^2
;                    ((1+Z)^2 + 1)^2
;
;             4*(1+Z)
;       = --------------- * uZ^2
;         ((1+Z)^2 + 1)^2


;Syntax
  if n_params() lt 2 then begin
    print, 'syntax: read_grand_vel, file, vel'
    print, "  e.g.: read_grand_vel, 'sunsim.08.10.vel', vel"
    return
  endif

;Test that file exists.
  if ~file_test(file) then begin
    message, /info, 'file not found: ' + file
    return
  endif

;Constants
  c = 2.99792458d8

;Determine number of lines in file.
  nobs = file_lines(file)

;Allocate array
  data = dblarr(5, nobs)

;Read velocities that grand wrote to disk.
  openr, unit, file, /get_lun
  readf, unit, data
  free_lun, unit

;Extract individual quantities from array.
  zn    = reform(data[1,*])	;measured bulk shift of entire order
  z0    = reform(data[2,*])	;input barycentric shift
  zbarn = reform(data[3,*])	;mean shift for lines after sigma clipping
  zsign = reform(data[4,*])	;std dev of line shifts after sigma clipping

;Calculate velocities.
  veln  = c * ((1+zn)^2       - 1) / ((1+zn)^2       + 1)
  v0    = c * ((1+z0)^2       - 1) / ((1+z0)^2       + 1)
  vbarn = c * ((1+zbarn)^2    - 1) / ((1+zbarn)^2    + 1)

;Propagate uncertainties from zsign to vbarn.
  uvbarn = zsign^2 * 4 * (1+zbarn) / ((1+zbarn)^2 + 1)^2

;Construct output structure.
  rec = $
    { v0    : 0d0 $
    , veln  : 0d0 $
    , vbarn : 0d0 $
    , uvbarn: 0d0 $
    , z0    : 0d0 $
    , zn    : 0d0 $
    , zbarn : 0d0 $
    , zsign : 0d0 $
    }
  vel = replicate(rec, nobs)

;Populate output structure.
  vel.v0     = v0
  vel.veln   = veln
  vel.vbarn  = vbarn
  vel.uvbarn = uvbarn
  vel.z0     = z0
  vel.zn     = zn
  vel.zbarn  = zbarn
  vel.zsign  = zsign

end
