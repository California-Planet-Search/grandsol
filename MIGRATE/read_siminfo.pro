pro read_siminfo, obsdir, obsfiles, info
;
;Purpose:
;  Read information from simulated observation files.
;
;Input:
;  obsdir (string) - directory containing simulated observations
;  obsfiles (string vector[nobs]) - names of observation files
;
;Output:
; info (structure) - information from observation files
;  .vsys (scalar) - systemic velocity from 'VSYST =' header line
;  .obsdir (string) - observation directory from 'RJDIR =' header line
;  .nobs (scalar) - number of observations described in obslist
;  .obsfile (string vector[nobs]) - names of observation files
;  .barycorr (double vector[nobs]) - [m/s] barycentric corrections
;  .vorb (float vector[nobs]) - [m/s] orbital velocities of companion

;Syntax
  if n_params() lt 3 then begin
    print, 'syntax: read_siminfo, obsdir, obsfiles, info'
    return
  endif

;Check that observation directory exists.
  if ~file_test(obsdir, /dir) then begin
    print, 'directory not found: ' + obsdir
    info = ''
    return
  endif

;Initialize return structure.
  nobs = n_elements(obsfiles)
  info = $
    { obsfile: strarr(nobs) $
    , star: '' $
    , iodine: '' $
    , elsf: '' $
    , echwave: '' $
    , echnorm: '' $
    , velsum: dblarr(nobs) $
    , zplus1: dblarr(nobs) $
    , pixsh: fltarr(nobs) $
    , ston: 0 $
    , lsfhwlo: 0 $
    , lsfhwhi: 0 $
    , osamp: 0 $
    , idlpro: '' $
    , seed_entry: lonarr(36, nobs) $
    , seed_exit: lonarr(36, nobs) $
    }
  infotags = strlowcase(tag_names(info))

;Loop through simulated observation files.
  for iobs=0, nobs-1 do begin

;Check that specified observation file exists.
    file = obsfiles[iobs]
    path_file = obsdir + path_sep() + file
    if ~file_test(path_file) then begin
      info = ''
      print, 'simulated observation not found: ' + file
      return
    endif
    info.obsfile[iobs] = file

;Read entire contents of .dsk file.
    rdsk_struct, path_file, dsk
    ndsk = n_tags(dsk)

;Verify that simulation header is present
    if ndsk lt 2 then begin
;     message, 'no simulation header in ' + file
      continue
    endif
    head = dsk.r2.data

;Parse each line assuming 'tag : value' format.
    nrec = n_elements(head)
    tags = strarr(nrec)
    values = strarr(nrec)
    for irec=0, nrec-1 do begin
      rec = head[irec]
      words = strsplit(rec, ':', /extract, count=nword)
      if nword ne 2 then message, 'bad tag:value syntax in header: ' + rec
      tags[irec] = strtrim(words[0], 2)
      values[irec] = strtrim(words[1], 2)
      if tags[irec] eq 'z+1' then tags[irec] = 'zplus1'
    endfor

;Loop through tags. Convert type. Check uniqueness. Save value.
    for irec=0, nrec-1 do begin
      tag = tags[irec]
      val = values[irec]
      case tag of
        'star'   : info.star = val
        'iodine' : info.iodine = val
        'elsf'   : info.elsf = val
        'echwave': info.echwave = val
        'echnorm': info.echnorm = val
        'velsum' : info.velsum[iobs] = double(val)
        'zplus1' : info.zplus1[iobs] = double(val)
        'pixsh'  : info.pixsh[iobs] = float(val)
        'ston'   : info.ston = round(float(val))
        'lsfhwlo': info.lsfhwlo = fix(val)
        'lsfhwhi': info.lsfhwhi = fix(val)
        'osamp'  : info.osamp = fix(val)
        'idlpro' : info.idlpro = val
        else: message, 'unrecognized tag in header: ' + tag
      endcase

;Check that scalar values remain identical.
      iinfo = (where(infotags eq tag, count))[0]
      if iobs gt 0 and n_elements(info.(iinfo)) eq 1 then begin
        if info.(iinfo) ne prev.(iinfo) then begin
          message, "'" + tag + "' is not identical for all observations"
        endif
      endif
    endfor

;Save random number seeds before and after each simulated observation.
    if ndsk ge 3 then begin
      seed = dsk.r3.data
      case n_elements(seed) of
        1: info.seed_entry[0,iobs] = seed
        36: info.seed_entry[*,iobs] = seed
        else: message, 'seed_entry does not have 1 or 36 elements'
      endcase
    endif
    if ndsk ge 4 then begin
      seed = dsk.r4.data
      case n_elements(seed) of
        1: info.seed_exit[0,iobs] = seed
        36: info.seed_exit[*,iobs] = seed
        else: message, 'seed_exit does not have 1 or 36 elements'
      endcase
    endif

;Done looping through observations.
    prev = info
  endfor

end
