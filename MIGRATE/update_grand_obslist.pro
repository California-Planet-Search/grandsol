pro update_grand_obslist, nord=nord, useveln=useveln

;syntax
  if n_params() lt 0 then begin
    print, 'syntax: update_grand_obslist [,nord= ,/useveln]'
    return
  endif

;constants
  c = 2.99792458d8

;determine previous obslist
  obslists = file_search('obslist*', count=nobslist)
  if nobslist eq 0 then message, 'no previous obslist file'
  obslist_prev = obslists[nobslist-1]
  iiter_prev = fix(strmid(obslist_prev, 7))

;read previous obslist file into structure
  read_grand_obslist, obli_prev, obslist_prev

;construct name for new obslist
  if iiter_prev gt 99 then message, 'maximum iternation is 99'
  obslist_new = 'obslist' + string(iiter_prev+1, form='(i2.2)')

;initialize next obslist structure
  obli_next = obli_prev
  obli_next.obslist = obslist_new
  obli_next.obs.vorb = 0

;handle case when creating obsfile01 for first iternation
  if iiter_prev le 0 then begin
    print, 'writing ' + obli_next.obslist
    write_grand_obslist, obli_next
    return
  endif

;get information about grand output files in current directory.
  sfx = 'vel'
  grand_files_iter, gf
  if ~keyword_set(gf[sfx].uord) then begin
    message, /info, "no '" + sfx $
           + "' files produced by grand in current directory"
    return
  endif

;read velocities for all orders
  if ~keyword_set(nord) then nord = n_elements(gf[sfx].uord)
  for iord=0, nord-1 do begin
    file_vel = gf[sfx].flast[iord]
    read_grand_vel, file_vel, vel
    if iord eq 0 then begin
      nobs = n_elements(vel)
      vels = replicate(vel[0], nobs, nord)
    endif
    vels[*,iord] = vel
  endfor

;calculate mean of all orders
  if nord gt 1 then begin
    iord_beg = -1
    mnvel = vel
    ntag = n_elements(tag_names(vel))
    for iobs=0, nobs-1 do begin
      for itag=0, ntag-1 do begin
        mnvel[iobs].(itag) = mean(vels[iobs,*].(itag), /double)
      endfor
    endfor
  endif else begin
    iord_beg = 0
  endelse

;loop through orders
  for iord=iord_beg, nord-1 do begin
    ord = strtrim(fix(gf[sfx].olast[iord]), 2)
    iter = gf[sfx].ilast

;data for current order (or for mean of all orders)
    vel = iord eq -1 ? mnvel : vels[*,iord]

;calculate orbital velocity, using special relativity
    betaprev = -vel.v0 / c
    betacurr = keyword_set(useveln) ? vel.veln / c : vel.vbarn / c
    betanew  = (betaprev+betacurr) / (1 + betaprev*betacurr)
    vorb = betanew * c

stop

;end of loop through orders.
  endfor

;add current and previous orbital velocities
  beta1 = vorb / c
  beta2 = obli_prev.obs.vorb / c
  beta  = (beta1+beta2) / (1 + beta1*beta2)
  vorb_next = beta * c

;update orbital velocity in obslist structure
  obli_next.obs.vorb = vorb_next

;write obslist
  print, 'writing ' + obli_next.obslist
  write_grand_obslist, obli_next

end
