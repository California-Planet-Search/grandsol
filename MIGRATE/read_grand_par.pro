pro read_grand_par, file, par, print=print
;
;Purpose:
; Read analysis parameters from a .par file produced by grand.
;
;Input:
; file (string) name of .par file containing analysis parameters
; [/print] (switch) print diagnostics to screen
;
;Output:
; par (structure) analysis parameters from a .par file

;Syntax
  if n_params() lt 2 then begin
    print, 'syntax: read_grand_par, file, par [,/print]'
    print, "  e.g.: read_grand_par, 'sunsim.08.par', par"
    return
  endif

;Test that file exists.
  if ~file_test(file) then begin
    message, /info, 'file not found: ' + file
    return
  endif

;Read parameter information from .par file
  nline = file_lines(file)
  lines = strarr(nline)
  openr, unit, file, /get_lun
  readf, unit, lines
  free_lun, unit

;Loop through lines, parsing parameter information.
  par = ''
  for iline=0, nline-1 do begin
    line = strtrim(lines[iline], 2)
    if line eq '' or strmid(line, 0, 1) eq '#' then continue

    words = strtrim(strsplit(line, ':', count=nword, /extract), 2)
    if nword ne 2 then begin
      message, /info, 'error parsing line ' + strtrim(iline, 2)
      message, line
    endif
    tag = words[0]
    value = words[1]

    if valid_num(value) then begin
      value = strpos(value, '.') ge 0 ? double(value) : long(value)
    endif

    if ~keyword_set(par) then begin
      par = create_struct(tag, value)
    endif else begin
      if max(strupcase(tag) eq tag_names(par)) eq 1 then continue
      par = create_struct(par, tag, value)
    endelse
  endfor

;add node locations for normalization vector
  par = create_struct(par $
    , 'xnrm', 1 + 20 * lindgen(200) $
    )

end
