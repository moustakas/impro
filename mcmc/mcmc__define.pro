;+
; NAME:
;  MCMC (IDL Class File)
;
; PURPOSE:
;  Calculate a Monte Carlo Markov Chain. 
;
; METHODS:
;  For detailed info for each method, use
;  IDL> doc_method, 'mcmc:methodname'
;
;  mcmc::run(step_func, like_func, nstep, parguess, like=, /log, seed=,
;            printstep=, file=)
;  mcmc::run, step_func, like_func, nstep, parguess, pars, like, /log, seed=,
;        printstep=, file=
;      Run a markov chain, optionally outputting to a file.  Both a functional
;      and procedural interface are provided.
;  mcmc::read_trials(file):  Read a set of trials from an output file.
;  mcmc::extract_stats(trialp, burning): Extract stats for each parameter from
;     a trials parameter array with the specified burn in.
;
; MODIFICATION HISTORY:
;  Created: 2005: Erin Sheldon, NYU.
;  Incorporated Hogg's more general method, and a bug fix thanks to Morad
;    Masjedi.  2006-10-14 Erin Sheldon, NYU
;
;-

FUNCTION mcmc::init
  return,1
END 

;docstart::mcmc::run
; NAME:
;  MCMC::RUN()
;  MCMC::RUN
;
; PURPOSE:
;  Run a Monte Carlo Markov Chain.  Both a functional and procedural interface
;  are provided.
;
; CALLING SEQUENCE:
;  trial_pars = mcmc->run(step_func, like_func, nstep, parguess, like=, 
;                        /log, seed=, file=)
;    OR
;  mcmc->run, step_func, like_func, nstep, parguess, trial_pars, like=, 
;        /log, seed=, file= 
;
; INPUTS:
;  step_func: A function that takes a step in the chain.
;  like_func: A function that computes the likelihood of an input set of
;             parameters. 
;  nstep: The number of steps to take in the Markov Chain.
;  parguess:  A starting guess at the parameters.
;
; KEYWORD PARAMETERS:
;  /log: "like_func" will calculate the log likelihood.  This is usually a
;        better way to go because of precision issues with exponentials in the
;        likelihood function.
;  seed: A starting seed.  If not input one is generated using IDL's builtin
;        seed generator from randomu.
;  printstep:  Print out a progress report every printstep trials.
;  file: Write the results to an output file rather than returning as a
;        variable. The value of trial_pars will be -1 in this case. The
;        results are read using the MCMC::READ_TRIALS() function.
;
; OUTPUTS:
;  trial_pars: An array [npars,ntrials] containing the trial steps through the
;              parameter space.  Will be -1 if file= is sent.
;
; OPTIONAL OUTPUTS:
;  like: The likelihood for each trial.
;
;
; EXAMPLES:
; 
; ; Fit a constant to n measurements.
;
; ; log likelihood.  
; function const_like, p
;   common mcmc_test_block, x, y, ivar, psigma, npars
;   chi2 = total( (y-p[0])^2*ivar )
;   return,-0.5*chi2
; end 
;
; function const_step, seed, pars
;   common mcmc_test_block, x, y, ivar, psigma, npars
;   return, pars + psigma*randomn(seed)
; end 
; function const_truepars, sigma=sigma
;   sigma = 1.0
;   return, 10.0
; end 
;
; ; For this test, just fake the data.
; pro const_setup
;   common mcmc_test_block, x, y, ivar, psigma, npars
;   ny = 10
;   pars = mcmc_test_const_truepars(sigma=sigma)
;   y = replicate(pars, ny)
;   yerr = replicate(sigma, ny)
;
;   ivar = 1.0/yerr^2
;   y[*] = y[*] + yerr*randomn(seed, ny)
;   psigma = yerr[0]
;   npars = 1
; end 
;
; IDL> mcmc=obj_new('mcmc')
; IDL> const_setup
; IDL> nstep = 10000
; IDL> parguess = 5.0 ; initial guess.  
; IDL> trials = mcmc->run('const_step', 'const_like', nstep, parguess, /log)
;
; ; Or you can write them to a file.  Better when it will take a long time
; ; or might crash.  This way you keep a record of what happened.
; IDL> mcmc->run, 'const_step', 'const_like', nstep, parguess, file=file, /log
; IDL> trials = mcmc->read_trials(file)
;
; MODIFICATION HISTORY:  
;  Created: 2005: Erin Sheldon, NYU.
;  Incorporated Hogg's more general method, and a bug fix thanks to Morad
;    Masjedi.  2006-10-14 Erin Sheldon, NYU
;
;docend::mcmc::run

; function version
function mcmc::run, step_func, like_func, nstep, parguess, like=like, log=log, seed=seed, file=file, printstep=printstep
  self->run, step_func, like_func, nstep, parguess, pars, like=like, log=log, seed=seed, file=file, printstep=printstep
  return, pars
end 
; procedure version
pro mcmc::run, step_func, like_func, nstep, parguess, pars, like=like, log=log, seed=seed, file=file, printstep=printstep

  if n_elements(seed) eq 0 then begin 
      seed = long( randomu(seed2)*10000 )
      tmp = randomu(seed)
  endif 

  npars   = n_elements(parguess)
  oldpars = float(parguess)

  ;; Should we write this to a file?
  if n_elements(file) ne 0 then begin 
      openw, lun, file, /get_lun

      ;; longs
      writeu, lun, long(nstep)
      writeu, lun, long(npars)

      pars = -1
  endif else begin 
      pars    = replicate(parguess[0],npars,nstep)
      like    = fltarr(nstep)
  endelse 

  for ii=0L,nstep-1L do begin

      self->_step, seed, oldpars, oldlike, step_func, like_func, newpars, newlike,$
        log=log

      if n_elements(file) ne 0 then begin 
          writeu, lun, float(newpars)
      endif else begin 
          pars[*,ii] = newpars
          like[ii] = newlike
      endelse 
      
      oldpars = newpars
      oldlike = newlike
      
      if n_elements(printstep) ne 0 then begin 
          if (ii mod  printstep) eq 0 then begin 
              print, $
                strjoin(strarr(21)+string(byte(8)),''), $
                'MCMC: ',100L*ii/nstep,' percent', $
                format= '($,A,A,I2.2,A)'
          endif 
      endif 
  endfor

  if n_elements(file) ne 0 then begin
      flush, lun
      free_lun, lun
  endif 

  if n_elements(printstep) ne 0 then begin 
      print, strjoin(strarr(21)+string(byte(8)),'')+'MCMC: done      '
  endif 
  return

end 

pro mcmc::_step, seed, pars, like, step_func, like_func, newpars, newlike, log=log

  if (not keyword_set(like)) then like= call_function(like_func,pars)

  newpars= call_function(step_func,seed,pars)
  newlike= call_function(like_func,newpars)

  if keyword_set(log) then likeratio= newlike-like $
  else likeratio= newlike/like

  randomnumber= randomu(seed)

  if keyword_set(log) then randomnumber= alog(randomnumber)
  if NOT ((newlike GT like) OR $
          (randomnumber LT likeratio)) then begin
      newpars= pars
      newlike= like
  endif

  return
end

;docstart::mcmc::read_trials
; NAME:
;  MCMC::READ_TRIALS()
;
; PURPOSE:
;  Read an output file from the MCMC ::run procedure.
;
; CALLING SEQUENCE:
;  trial_pars = mcmc->read_trials(file)
;
; INPUTS:
;  file: A file created by the MCMC::RUN procedure with file=file input.
;
; OUTPUTS:
;  trial_pars: An array [npars,ntrials] containing the trial steps through the
;              parameter space. 
;
;
;
; MODIFICATION HISTORY:  
;  Created: 2005: Erin Sheldon, NYU.
;
;docend::mcmc::read_trials


FUNCTION mcmc::read_trials, file
  openr, lun, file, /get_lun
  ntrial=0L
  npar=0L

  readu, lun, ntrial
  readu, lun, npar

  pars = fltarr(npar, ntrial)
  readu, lun, pars
  free_lun, lun

  return, pars

END 


;docstart::mcmc::extract_stats
; NAME:
;  MCMC::EXTRACT_STATS()
;
; PURPOSE:
;  Extract statistics from an MCMC chain with burn in.
;
; CALLING SEQUENCE:
;  stats = mcmc->extract_status(trial_pars, burnin)
;
; INPUTS:
;  trial_pars: The outputs of the MCMC chain.
;  burnin: Number of points from the beginning to skip for burn in.
;
; OUTPUTS:
;  Mean and standard deviation for each parameter.  A [npar,2] array.
;
; MODIFICATION HISTORY:  
;  Created: 2005: Erin Sheldon, NYU.
;
;docend::mcmc::extract_stats


FUNCTION mcmc::extract_stats, trialp, like, burnin

    nn = (n_elements(trialp) gt 0) + (n_elements(like) gt 0) + (n_elements(burnin) gt 0)
    IF nn ne 3 then begin
        on_error, 2
        print,'-Syntax: parms = mcmc->extract_stats(trialp, like, burnin)'
        print
        message,'Halting'
    ENDIF 

    sst = size(trialp,/struct)

    if sst.n_dimensions ne 2 then begin 
        on_error,2
        message,'trialp must be a [npar,ntrial] array'
    endif 


    npar = sst.dimensions[0]
    ntrial = sst.dimensions[1]

    parms = dblarr(npar, 2)

    w=lindgen(ntrial)
    w = w[burnin:ntrial-1]
    nuse = n_elements(w)

    maxl = max(like[w], best)
    FOR i=0L, npar-1 DO BEGIN 

        value = trialp[i,w[best]] 

        ; calc sdev around the best fit
        ; Numerically-stable "two-pass" formula, which offers less
        ; round-off error. Page 613, Numerical Recipes in C.
        resid = trialp[i,w] - value
        var = (TOTAL(resid^2, /double) - $
               (TOTAL(resid, /double)^2)/nuse)/(nuse-1.0)

        parms[i,0] = value
        parms[i,1] = sqrt(var)
    ENDFOR 

    return, parms
END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Some plotting routines
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mcmc::_pmulti_value, nplot
    if nplot eq 1 then return, 0
    if nplot eq 2 then return,[0,0,2]
    if nplot eq 3 or nplot eq 4 then return,[0,2,2]
    if nplot eq 5 or nplot eq 6 then return,[0,3,2]
    if nplot ge 7 and nplot le 9 then return,[0,3,3]
    if nplot ge 10 and nplot le 12 then return,[0,4,3]
    if nplot ge 13 then return,[0,4,4]
end


; all 1d marginalizations
PRO mcmc::plot1d, trials, like, burnin, names=names, nsmooth=nsmooth, _extra=_extra

    nn = (n_elements(trials) gt 0) + (n_elements(like) gt 0) + (n_elements(burnin) gt 0)
    if nn ne 3 then begin
        print,'-Syntax: mcmc->plot1d, trials, like, burnin, names=, nsmooth=, _extra='
        on_error, 2
        message,'Halting'
    endif

    if n_elements(nsmooth) eq 0 then nsmooth=2
    sz=size(trials, /dim)
    npar = sz[0]
    ntrial = sz[1]

    xmold = !x.margin
    ymold = !y.margin
    pmold = !p.multi
    !p.multi = self->_pmulti_value(npar)
    if npar gt 2 then begin
        !x.margin = [4.,1.]
        !y.margin = [1.,0.5]
    endif

    nnames=n_elements(names)
    if nnames ne 0 then begin
        if nnames ne npar then message,'size of names must be number of parameters'
    endif else begin
        names = 'p'+strtrim(lindgen(npar),2)
;       names = 'p'+ntostr(lindgen(npar))
    endelse

    ; get the stats
    parms = self->extract_stats(trials, like, burnin)

    w=lindgen(ntrial)
    w = w[burnin:ntrial-1]
    nsig = 4.0
    for i=0L, npar-1 do begin
        mn = parms[i,0]
        sdev = parms[i,1]
        bin = sdev/7.0

        xrange = mn + nsig*sdev*[-1.0,1.0]
        plothist, trials[i,w], xhist, yhist, bin=bin, min=xrange[0], max=xrange[1], $
            /noplot
        yhist = smooth(float(yhist),nsmooth)
        yhist = yhist/max(yhist)
        plot, xhist, yhist, psym=10, $
            yrange=[0,1.2], ystyle=1, xtitle=names[i], $
            ytickname=replicate(' ',20)
        oplot, [mn,mn],[0,10], color=djs_icolor('blue')
        oplot, [mn+sdev,mn+sdev],[0,10], color=djs_icolor('blue'), line=2
        oplot, [mn-sdev,mn-sdev],[0,10], color=djs_icolor('blue'), line=2
;       pplot, xhist, yhist, psym=10, aspect=1, $
;           yrange=[0,1.2], ystyle=1, xtitle=names[i], $
;           ytickname=replicate(' ',20)
;       oplot, [mn,mn],[0,10], color=!blue
;       oplot, [mn+sdev,mn+sdev],[0,10], color=!blue, line=2
;       oplot, [mn-sdev,mn-sdev],[0,10], color=!blue, line=2

    endfor

    !x.margin=xmold
    !y.margin=ymold
    !p.multi=0

END 

; all 2d marginalizations
PRO mcmc::plot2d, trials, like, burnin, nbin=nbin, nsmooth=nsmooth, names=names, _extra=_extra

    nn = (n_elements(trials) gt 0) + (n_elements(like) gt 0) + (n_elements(burnin) gt 0)
    if nn ne 3 then begin
        print,'-Syntax: mcmc->plot1d, trials, like, burnin, nbin=, nsmooth=, names=, _extra='
        on_error, 2
        message,'Halting'
    endif

    sz=size(trials, /dim)
    npar = sz[0]
    if npar eq 1 then message,'There must be at least two parameters for 2-d plots'
    nplots = npar-1

    if n_elements(nsmooth) eq 0 then nsmooth=2

    xmold = !x.margin
    ymold = !y.margin
    pmold = !p.multi
    !p.multi = self->_pmulti_value(nplots)
    if nplots gt 2 then begin
        !x.margin = [4.,1.]
        !y.margin = [1.,0.5]
    endif

    ntrial = sz[1]

    if n_elements(nbin) eq 0 then nbin = 40
    nx=nbin
    ny=nbin
    
    nnames=n_elements(names)
    if nnames ne 0 then begin
        if nnames ne npar then message,'size of names must be number of parameters'
    endif else begin
        names = 'p'+strtrim(lindgen(npar),2)
;       names = 'p'+ntostr(lindgen(npar))
    endelse

    ; get the stats
    parms = self->extract_stats(trials, like, burnin)

    w=lindgen(ntrial)
    w = w[burnin:ntrial-1]
    nsig = 5.0
    for i=0L, nplots-1 do begin
        p1 = i
        p2 = i+1

        mn1 = parms[p1,0]
        sdev1 = parms[p1,1]
        bin1 = sdev1/10.0
        mn2 = parms[p2,0]
        sdev2 = parms[p2,1]
        bin2 = sdev2/10.0

        ; get region within nsig sigma
        ;range1 = mn1 + nsig*sdev1*[-1.0,1.0]; > min(trials[p1,w]) < max(trials[p1,w])
        ;range2 = mn2 + nsig*sdev2*[-1.0,1.0]; > min(trials[p2,w]) < max(trials[p2,w])
        ;hist=histogram_2d(trials[p1,w], trials[p2,w], nx, ny,$
        ;    xrange=range1,yrange=range2)
        hist=histogram_2d(trials[p1,w], trials[p2,w], nx, ny)

        hist.map = smooth(hist.map, nsmooth)

        Lrel = hist.map/max(hist.map)
        ;; convert chi^2 levels to likelihood ratio levels
        Llevels = exp( -0.5*!siglevels2 )
        Llevels = Llevels[sort(Llevels)]

        contour, Lrel, hist.xbins, hist.ybins, xtitle=names[p1], ytitle=names[p2], $
            levels=Llevels
        oplot, [mn1], [mn2], psym=7, symsize=2, color=!red

        ;oplot, mn1+sdev1*[-1,-1], [-1000, 1000], color=!blue, line=2
        ;oplot, mn1+sdev1*[1,1], [-1000, 1000], color=!blue, line=2

    endfor

    !x.margin=xmold
    !y.margin=ymold
    !p.multi=pmold

END 





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; cleanup
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO mcmc__define
  struct = { $
             mcmc, $
             mcmc_dummy_var: 0 $
           }
END 
