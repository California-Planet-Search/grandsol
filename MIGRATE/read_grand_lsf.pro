pro read_grand_lsf, xq, lsfq, xj, lsfj, lsffile, obs, ord, k
;
;Purpose:
; Read LSF node values from a .lsf file produced by grand.
; Calculate supersampled LSF from node values.
;
;Input:
; lsffile (string) name of .lsf file containing node values
; obs (integer) observation number [1-100]
; ord (integer) echelle order number [8]
; k (integer) region number within an order [1-9]
;
;Output:
; xq (vector[qdim]) pixel offsets for lsfq
; lsfq (vector[qdim]) lsf node values
; xj (vector[jdim]) pixel offsets for lsfj
; lsfj (vector[jdim]) supersampled lsf
;
;Notes:
; Calls fortran program.

;Syntax
  if n_params() lt 8 then begin
    print, 'syntax: read_grand_lsf, xq, lsfq, xj, lsfj, lsffile, obs, ord, k'
    print, "  e.g.: read_grand_lsf, xq, lsfq, xj, lsfj, 'sunsim.lsf', 50, 8, 5"
    return
  endif

;Test that file exists.
  if ~file_test(lsffile) then begin
    message, /info, 'file not found: ' + lsffile
    return
  endif

;Call fortran program.
  obsstr = strtrim(obs)
  ordstr = strtrim(ord)
  kstr = strtrim(k)
  spawn, /noshell, ['./makelsf', lsffile, obsstr, ordstr, kstr], stdout

;Extract data from fortran output.
  a = intarr(5)
  reads, stdout[0], a
  qdim = a[0]
  jdim = a[1]
  lsfq = fltarr(qdim)
  reads, stdout[1], lsfq
  lsfj = fltarr(jdim)
  reads, stdout[2], lsfj

;generate pixel offsets
  osamp = 10
  xq = findgen(qdim) - 0.5*(qdim-1)
  xj = (findgen(jdim) - 0.5*(jdim-1)) / osamp

end
