pro iterate_grand, id, olist, niter, useveln=useveln, fudge=fudge

;syntax
  if n_params() lt 3 then begin
    print, 'syntax: iterate_grand, id, olist, niter [,/useveln ,/fudge]'
    print, "  e.g.: iterate_grand, 'sunsim', [6,7,8], 10"
    return
  endif

;internal parameters
  maxwait = 7200			;seconds

;computed qunatities
  mlist = string(olist, form='(i2.2)')
  nmlist = n_elements(mlist)
  pids = strarr(nmlist)

;begin main loop
  for iiter=1, niter-1 do begin
    iter = string(iiter, form='(i2.2)')
    obslist = 'obslist' + iter
    label = id + iter
    flags = '211111'
    fudge = keyword_set(fudge) ? 'fudge+' : 'fudge-'
    vorb = 'vorb+'

;update observation list
    update_grand_obslist, useveln=keyword_set(useveln)

;loop through orders
    for imlist=0, nmlist-1 do begin
      m = mlist[imlist]
      out = label + '.' + m + '.log'

;construct grand solution command with arguments
      cmd = [ 'grand' $
            , obslist $
            , label $
            , m $
            , 'out='+out $
            , flags $
            , fudge $
            , vorb $
          ]
      fullcmd = strjoin(cmd, ' ')

;call grand solution for current iteration.
      spawn, fullcmd+' &', stdout, stderr
      if keyword_set(stderr) or n_elements(stdout) ne 1 $
                             or strmid(stdout, 0, 4) ne '[1] ' then begin
        print, fullcmd
        print, 'stdout: ' + stdout
        print, 'stderr: ' + stderr
        stop
      endif
      pids[imlist] = strtrim(strmid(stdout[0], 4), 2)
      print, pids[imlist] + ': ' + fullcmd

;end order loop
    endfor

;wait until spawned jobs for all orders are done
    time0 = systime(/sec)
    cmd = ['ps', '-p', pids]
    repeat begin
      done = 0
      wait, 60
      dt = long(systime(/sec)-time0)
      spawn, /nosh, cmd, stdout, stderr
      if strpos(stdout[0], 'PID TTY') lt 0 then begin
        message, 'error getting status of spawned processes for ' + label
      endif
      if n_elements(stdout) gt 1 then begin
        running = strmid(stdout[1:*], 0, 5)
        print, strtrim(dt, 2) $
             + 's: still waiting for pids ' $
             + strjoin(running, ' ')
      endif else done = 1
      if dt gt maxwait then begin
        print, 'maxtime (' + strtrim(maxwait, 2) $
             + ' exceeded waiting for ' + label
        stop
      endif
    endrep until done

;end main loop
  endfor

;one final observation list (with velocities from final iteration)
  update_grand_obslist

end
