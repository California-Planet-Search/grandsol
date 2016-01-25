function fmtnum, data, fmt
;Return a string containing DATA formatted according to FMT with the
;following additional rules:
; 1) strip all leading and trailing spaces
; 2) strip all trailing zeros after the decimal point
; 3) strip trailing decimal point

if n_params() lt 2 then begin
  print, 'syntax: string = fmtnum(data, fmt)'
  retall
endif

  zero = (byte('0'))(0)			;ASCII for zero

;Convert to string.
  str = string(data, form=fmt)		;use requested format
  str = strtrim(str, 2)			;strip leading/trailing spaces

;Return if integer or exponential format was used.
  iexp = strpos(strlowcase(str), 'e')	;look for exponent
  idec = strpos(str, '.')		;look for decimal point
  if iexp gt 0 or idec lt 0 then begin
    return, str
  endif

;Strip trailing zeros after decimal point.
  byt = byte(str)			;get byte array equivalent
  len = max(where(byt ne zero)) + 1	;length desired after stripping
  str = strmid(str, 0, len)		;strip trailing zeros

;Strip trailing decimal point.
  idec = strpos(str, '.')		;look for decimal point
  len = strlen(str)			;string length
  if idec+1 eq len then begin
    str = strmid(str, 0, len-1)		;strip trailing decimal point
  endif
  return, str

end
