pro read_grand_tem, file, tem, print=print
;
;Purpose:
; Read template stellar spectrum from a .tem file produced by grand.
;
;Input:
; file (string) name of .tem file containing a stellar template spectrum
; [/print] (switch) print diagnostics to screen
;
;Output:
; tem (structure[ntem]) template stellar spectrum from a .tem file
;  .w - wavelengths for template spectrum [Angstroms]
;  .t - template spectrum
;  .tz - template spectrum shifted to solar frame
;  .sun - solar spectrum

;Syntax
  if n_params() lt 2 then begin
    print, 'syntax: read_grand_tem, file, tem [,/print]'
    print, "  e.g.: read_grand_tem, 'sunsim.08.99.tem', tem"
    return
  endif

;Test that file exists.
  if ~file_test(file) then begin
    message, /info, 'file not found: ' + file
    return
  endif

;Determine number of rows in file.
  ntem = file_lines(file)

;Allocate array
  data = dblarr(5, ntem)

;Read stellar template spectrum that grand wrote to disk.
  openr, unit, file, /get_lun
  readf, unit, data
  free_lun, unit

;Print diagnostics
  if keyword_set(print) then begin
    wmin = min(data[1,*], max=wmax)
    print, file + ', ' $
         + strtrim(ntem, 2) + ' nodes, wmin=' $
         + strtrim(string(wmin, form='(f9.2)'), 2) + '-' $
         + strtrim(string(wmax, form='(f9.2)'), 2) + ' Angstroms'
  endif

;Construct output structure.
  rec = $
    { w   : 0d0 $
    , t   : 0.0 $
    , tz  : 0.0 $
    , sun : 0.0 $
    }
  tem = replicate(rec, ntem)

;Populate output structure.
  tem.w   = reform(data[1,*])
  tem.t   = reform(data[2,*])
  tem.tz  = reform(data[3,*])
  tem.sun = reform(data[4,*])

end
