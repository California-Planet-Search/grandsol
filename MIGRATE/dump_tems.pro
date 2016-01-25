pro dump_tems, dir

  if ~keyword_set(dir) then dir = '.'

  fpatt = dir + path_sep() + '*.tem'
  files = file_search(fpatt, count=nfile)
  print, strtrim(nfile,2) + ' tem files found'
  if nfile eq 0 then return

  for ifile=0, nfile-1 do begin
    file = files[ifile]
    basename = file_basename(file)
    words = strsplit(basename, '.', count=nword, /extract)
    word0 = words[0]
    len0 = strlen(word0)
    iterstr = strmid(word0, len0-2)
    iter = fix(iterstr)

    read_grand_tem, file, tem
    if ifile eq 0 then begin
      nw = n_elements(tem)
      tems = replicate(tem[0], nw, nfile)
    endif
    tems[*,ifile] = tem
  endfor

  file = 'dump_tems.sav'
  print, 'saving ' + file
  save, file=file, tems

end
