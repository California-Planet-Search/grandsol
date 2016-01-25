pro grand_files_iter, gf, print=print
;
;Notes:
; n_elements(gf) returns the number of suffixes
; sfxs = gf.keys() gets list of suffixes with file information
; if gf.haskey('vel') test if info exists for a particular suffix 
; gf['iter'].ord lists the orders that have .par files
; gf['iter'].iter lists the last iteration for each order with a .par file

;Syntax
  if n_params() lt 1 then begin
    print, 'syntax: grand_files_iter, gf [,/print]'
    return
  endif

;List of file suffixes
  sfxs = ['cnr', 'lsf', 'mod', 'nim', 'par', 'tem', 'vel', 'vln', 'vsl', 'vsn']
  nsfx = n_elements(sfxs)

;Initialize storage.
  gf = hash()

;Loop through suffixes.
  for isfx=0, nsfx-1 do begin
    sfx = sfxs[isfx]

;Find files with current suffix.
    fpatt = '*.' + sfx
    files = file_search(fpatt, count=nfile)

;Parse filenames that have order number, but no iteration number.
    if sfx eq 'par' then begin
      regex = '(.*)' $			;all characters before what follows
            + '([0-9][0-9]\.)' $	;two digit outer loop number, e.g. '05.'
            + '([0-9][0-9]\.)' $	;two digit order number, e.g. '08.'
            + '(' + sfx + ')'		;suffix without a period, e.g. 'cnr'
      comps = stregex(files, regex, /subex, /extr)
      ids = reform(comps[1,*])
      loops = strmid(reform(comps[2,*]), 0, 2)
      ords = strmid(reform(comps[3,*]), 0, 2)
      its = replicate('', n_elements(ords))
    
;Parse using regular expressions.
    endif else begin
      regex = '(.*)' $			;all characters before what follows
            + '([0-9][0-9]\.)' $	;two digit outer loop number, e.g. '05.'
            + '([0-9][0-9]\.)' $	;two digit order number, e.g. '08.'
            + '([0-9][0-9]\.)' $	;iteration number, e.g. '10.'
            + '(' + sfx + ')'		;suffix without a period, e.g. 'cnr'
      comps = stregex(files, regex, /subex, /extr)
      ids = reform(comps[1,*])
      loops = strmid(reform(comps[2,*]), 0, 2)
      ords = strmid(reform(comps[3,*]), 0, 2)
      its = strmid(reform(comps[4,*]), 0, 2)
    endelse

;Get list of unique id, outer loop number, order number, iteration number.
    uids   = ids[uniq(ids, sort(ids))]
    uloops = loops[uniq(loops, sort(loops))]
    uords  = ords[uniq(ords, sort(ords))]
    uits   = its[uniq(its, sort(its))]
    nuid   = n_elements(uids)
    nuloop = n_elements(uloops)
    nuord  = n_elements(uords)
    nuit   = n_elements(uits)

;Put filenames into an nuord-by-nuit grid.
    fgrid = strarr(nuloop, nuord, nuit)
    ogrid = strarr(nuloop, nuord, nuit)
    igrid = strarr(nuloop, nuord, nuit)
    for i=0, nfile-1 do begin
      iuloop = where(uloops eq loops[i])
      iuord = where(uords eq ords[i])
      iuit = where(uits eq its[i])
      if fgrid[iuloop,iuord,iuit] ne '' then begin
        print, 'WARNING - ignoring multiple files with same ORD, IT, and SFX!!'
      endif
      fgrid[iuloop,iuord,iuit] = files[i]
      ogrid[iuloop,iuord,iuit] = ords[i]
      igrid[iuloop,iuord,iuit] = its[i]
    endfor

;Create lists with last iteration info for each order.
    flast = strarr(nuord)
    olast = strarr(nuord)
    ilast = strarr(nuord)
    for iuord=0, nuord-1 do begin
      max = max(igrid[*,iuord,*], imax)
      flast[iuord] = fgrid[*,iuord,imax]
      olast[iuord] = ogrid[*,iuord,imax]
      ilast[iuord] = max
    endfor

;Print diagnostics.
    if keyword_set(print) then begin
      if keyword_set(uits) then begin
        print, sfx + ': ' $
          + strtrim(nuid,2)   + ' id(' + strjoin(uids,',')   + '), ' $
          + strtrim(nuloop,2) + ' loop(' + strjoin(uloops,',') + '), ' $
          + strtrim(nuord,2)  + ' ord(' + strjoin(uords,',')  + '), ' $
          + strtrim(nuit,2)   + ' it(' + strjoin(uits,',')   + ')'
      endif else begin
        print, sfx + ': ' $
          + strtrim(nuid,2)   + ' id(' + strjoin(uids,',')   + '), ' $
          + strtrim(nuloop,2) + ' loop(' + strjoin(uloops,',') + '), ' $
          + strtrim(nuord,2)  + ' ord(' + strjoin(uords,',')  + ')'
      endelse
    endif

;Save results
    fi = $
      { uid: uids $
      , uloop: uloops $
      , uord: uords $
      , uit: uits $
      , ilast: ilast $
      , olast: olast $
      , flast: flast $
      , fgrid: fgrid $
      }
    gf[sfx] = fi

;Done looping through suffixes.
  endfor

;Get complete list of orders.
  ords = gf['par'].uord
  nord = n_elements(ords)

;Determine the last explicitly listed iteration for each order
  last = replicate(-1, nord)
  keys = gf.keys()
  nkey = n_elements(keys)
  for ikey=0, nkey-1 do begin
    key = keys[ikey]
    if max(gf[key].uord ne ords) eq 1 then continue
    if ~keyword_set(gf[key].ilast) then continue
    ilast = fix(gf[key].ilast)
    iwhr = where(ilast ne 99, nwhr)
    if nwhr gt 0 then last[iwhr] = last[iwhr] > ilast[iwhr]
  endfor
  gf['iter'] = $
    { ord  : ords $
    , last : string(last, form='(i2.2)') $
    }
  
end
