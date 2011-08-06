;+
; NAME:
;   MONTE_LOG12OH_KK04
;
; PURPOSE:
;   Compute oxygen abundances using the KK04 calibration, with
;   proper Monte Carlo of the errors.
;
; INPUTS: 
;   r23, o32 - if you don't know what these are you probably shouldn't
;     be using this routine [NOBJ]
;
; OPTIONAL INPUTS: 
;   nmonte - number of Monte Carlo realizations; note that if
;     NMONTE=0 then LOG12OH_AVG as well as the errors will not be
;     computed (default 250)
;   maxiter - maximum number of iterations (default 50)
;   seed - random number seed
;   line - emission-line structure containing the fluxes and errors of
;     the [OII], [OIII] 4959+5007, and H-beta emission lines from which
;     R23 and O32 were derived [NOBJ]:
;
;     IDL> help, line, /str
;     ** Structure <1544078>, 3 tags, length=24, data length=24, refs=1:
;        OII         FLOAT     Array[2]
;        OIII        FLOAT     Array[2]
;        HBETA       FLOAT     Array[2]
;
; KEYWORD PARAMETERS: 
;   log - indicates that R23 and O32 are in the log, i.e., log(R23),
;     log(O32) 
;   debug - render a QAplot and wait for a keystroke
;
; OUTPUTS: 
;   result - output data structure (mostly self-explanatory)
;
; OPTIONAL OUTPUTS:
;   result_monte - RESULT, for each Monte Carlo realization 
;
; COMMENTS:
;   Note that LINE is only required if NMONTE>0.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Feb 06, NYU, written
;   jm08mar21nyu - allow NMONTE=0
;   jm09sep25ucsd - added SEED optional inpu
;   jm10jul31ucsd - near-complete rewrite of uncertainties and
;     convergence testing; slower, but correct!
;
; Copyright (C) 2008-2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function log12oh_kk04, x, y, maxiter=maxiter, silent=silent
; driver routine to find the metallicity and log U iteratively 
    
    light = 2.99792458D10 ; speed of light [cm/s]

    nobj = n_elements(x)
    res = replicate({log12oh_upper: 0.0, log12oh_lower: 0.0, $
      log12oh_avg: 0.0, logu_upper: 0.0, logu_lower: 0.0, $
      logu_avg: 0.0, converge_upper: 1, converge_lower: 1},nobj)

    for ii = 0L, nobj-1L do begin
; upper
       logoh_upper = 9.0
       logq_upper = 7.5
       iter = 0
       quit = 0
       while (quit eq 0) and (iter lt maxiter) do begin
          logoh_iter = logoh_upper
          logq_iter = logq_upper
       
          logq_upper = (32.81D - 1.153*y[ii]^2 + logoh_upper*$
            (-3.396 - 0.025*y[ii] + 0.1444*y[ii]^2)) / $
            (4.603D - 0.3119*y[ii] - 0.163*y[ii]^2 + logoh_upper*$
            (-0.48 + 0.0271*y[ii] + 0.02037*y[ii]^2))
          logoh_upper = 9.72D - 0.777*x[ii] - 0.951*x[ii]^2 - $
            0.072*x[ii]^3 - 0.811*x[ii]^4 - $
            logq_upper*(0.0737 - 0.0713*x[ii] - 0.141*x[ii]^2 + $
            0.0373*x[ii]^3 - 0.058*x[ii]^4)

          if (abs(logoh_upper-logoh_iter) lt 0.01) and $
            (abs(logq_upper-logq_iter) lt 0.01) then quit = 1
;         print, iter, logoh_upper, logoh_iter, logoh_upper-logoh_iter, $
;           logq_upper, logq_iter, logq_upper-logq_iter
          iter++
       endwhile
       res[ii].log12oh_upper = logoh_upper
       res[ii].logu_upper = logq_upper - alog10(light)
       if (iter eq maxiter) then res[ii].converge_upper = 0
       if (keyword_set(silent) eq 0) then if (res[ii].converge_upper eq 0) then $
         splog, 'KK04/upper on object '+strtrim(ii,2)+' did not converge!'

; lower       
       logoh_lower = 8.0
       logq_lower = 7.5
       iter = 0 & quit = 0
       while (quit eq 0) and (iter lt maxiter) do begin
          logoh_iter = logoh_lower
          logq_iter = logq_lower

          logq_lower = (32.81D - 1.153D*y[ii]^2 + logoh_lower*$
            (-3.396D - 0.025D*y[ii] + 0.1444D*y[ii]^2)) / $
            (4.603D - 0.3119D*y[ii] - 0.163D*y[ii]^2 + logoh_lower*$
            (-0.48D + 0.0271D*y[ii] + 0.02037D*y[ii]^2))
          logoh_lower = 9.40D + 4.65D*x[ii] - 3.17D*x[ii]^2 - $
            logq_lower*(0.272D + 0.547D*x[ii] - 0.513D*x[ii]^2)

          if (abs(logoh_lower-logoh_iter) lt 0.01) and $
            (abs(logq_lower-logq_iter) lt 0.01) then quit = 1
;         print, iter, logoh_lower, logoh_iter, logoh_lower-logoh_iter, $
;           logq_lower, logq_iter, logq_lower-logq_iter
          iter++
       endwhile
       res[ii].log12oh_lower = logoh_lower
       res[ii].logu_lower = logq_lower - alog10(light)
       if (iter eq maxiter) then res[ii].converge_lower = 0
       if (keyword_set(silent) eq 0) then if (res[ii].converge_lower eq 0) then $
         splog, 'KK04/lower on object '+strtrim(ii,2)+' did not converge!'
    endfor 
    res.log12oh_avg = (res.log12oh_upper+res.log12oh_lower)/2.0
    res.logu_avg = (res.logu_upper+res.logu_lower)/2.0
    
return, res
end
    
function monte_log12oh_kk04, allr231, allo321, log=log, nmonte=nmonte, $
  line=line, result_monte=allresult_monte, seed=seed, maxiter=maxiter, $
  debug=debug
    
    if (n_elements(nmonte) eq 0) then nmonte = 250
    if (n_elements(seed) eq 0) then delvarx, seed
    if (n_elements(maxiter) eq 0) then maxiter = 50
    
    nallobj = n_elements(allr231)
    
    allresult = {$
;     nmonte:                 0,$
;     frac:              -999.0,$
      log12oh_avg:       -999.0,$
      log12oh_avg_err:   -999.0,$
      log12oh_upper:     -999.0,$
      log12oh_upper_err: -999.0,$
      log12oh_lower:     -999.0,$
      log12oh_lower_err: -999.0,$
      logu_avg:          -999.0,$
      logu_avg_err:      -999.0,$
      logu_upper:        -999.0,$
      logu_upper_err:    -999.0,$
      logu_lower:        -999.0,$
      logu_lower_err:    -999.0,$
      converge_upper:         1,$
      converge_lower:         1}
    allresult = replicate(allresult,nallobj)

    if (nmonte gt 0L) and arg_present(allresult_monte) then begin
       allresult_monte = {log12oh_upper: fltarr(nmonte), $
         log12oh_lower: fltarr(nmonte), logq_upper: fltarr(nmonte), $
         logq_lower: fltarr(nmonte)}
       allresult_monte = replicate(allresult_monte,nallobj)
    endif

; check for good measurements    
    these = where((allr231 gt -900.0) and (allo321 gt -900.0),nobj)
    if (nobj eq 0L) then return, allresult else begin
       r231 = allr231[these]
       o321 = allo321[these]
       result = allresult[these]
       if (nmonte gt 0L) and arg_present(allresult_monte) then $
         result_monte = allresult_monte[these]
    endelse

    if keyword_set(log) then begin
       r23 = 10.0^r231
       o32 = 10.0^o321
    endif else begin
       r23 = r231
       o32 = o321
    endelse

; compute the abundances on the lower and upper branches
    x = alog10(r23)
    y = alog10(o32)

    res = log12oh_kk04(x,y,maxiter=maxiter)
    result = im_struct_assign(res,result,/nozero)

; compute the average metallicity and ionization parameter
    result.log12oh_avg = (result.log12oh_upper+result.log12oh_lower)/2.0
    result.logu_avg = (result.logu_upper+result.logu_lower)/2.0
    
; use Monte Carlo to figure out the uncertainties; use the linear
; uncertainties on the individual line-fluxes to construct the Monte
; Carlo distributions of R23 and O32, since the two quantities are
; correlated; see IM_ABUNDANCE() for an example of what the LINE
; structure should look like
    if (nmonte gt 0L) then begin
       if (n_elements(line) eq 0L) then message, $
         'Monte Carlo uncertainties require LINE'
       
       oiimonte = rebin(reform(line.oii[0],1,nobj),nmonte,nobj)+randomn(seed,nmonte,nobj)*$
         rebin(reform(line.oii[1],1,nobj),nmonte,nobj)
       oiiimonte = rebin(reform(line.oiii[0],1,nobj),nmonte,nobj)+$
         randomn(seed,nmonte,nobj)*rebin(reform(line.oiii[1],1,nobj),nmonte,nobj)
       hbmonte = rebin(reform(line.hbeta[0],1,nobj),nmonte,nobj)+randomn(seed,nmonte,nobj)*$
         rebin(reform(line.hbeta[1],1,nobj),nmonte,nobj)

       r23monte = (oiimonte + oiiimonte) / hbmonte
       o32monte = oiiimonte / oiimonte

; old, wrong code!       
;      r23monte = rebin(reform(r23,1,nobj),nmonte,nobj)+randomn(seed,nmonte,nobj)*$
;        rebin(reform(r23err,1,nobj),nmonte,nobj)
;      o32monte = rebin(reform(o32,1,nobj),nmonte,nobj)+randomn(seed,nmonte,nobj)*$
;        rebin(reform(o32err,1,nobj),nmonte,nobj)

       for ii = 0L, nobj-1L do begin
; enforce positivity
          montepos = where((r23monte[*,ii] gt 0.0) and $
            (o32monte[*,ii] gt 0.0),nmontepos)
;         if (nmontepos ne nmonte) then message, 'Problem here'
          if (nmontepos gt 1L) then begin

             xmonte = alog10(r23monte[montepos,ii])
             ymonte = alog10(o32monte[montepos,ii])

             res = log12oh_kk04(xmonte,ymonte,maxiter=maxiter,/silent)
             conv = where(res.converge_upper and res.converge_lower,$
               nconv,comp=fail,ncomp=nfail)
             if (nconv eq 0) then begin
                splog, 'Monte Carlo errors on object '+$
                  strtrim(ii,2)+' did not converge!'
             endif else begin
                result[ii].log12oh_upper_err = djsig(res[conv].log12oh_upper)
                result[ii].log12oh_lower_err = djsig(res[conv].log12oh_lower)
                result[ii].logu_upper_err = djsig(res[conv].logu_upper)
                result[ii].logu_lower_err = djsig(res[conv].logu_lower)

                good = where((res[conv].log12oh_lower lt res[conv].log12oh_upper),ngood)
                if (ngood eq 0L) then begin ; average abundance error not defined
                endif else begin
                   result[ii].log12oh_avg_err = djsig([res[conv[good]].log12oh_lower,$
                     res[conv[good]].log12oh_upper])
                   result[ii].logu_avg_err = djsig([res[conv[good]].logu_lower,$
                     res[conv[good]].log12oh_upper])
                endelse
             endelse
          endif

; debugging plot
          if keyword_set(debug) then begin
             djs_plot, alog10(r23monte[*,ii]), res.log12oh_upper, ps=6, xsty=3, ysty=3, $
               yr=[7,9.5], xrange=[-0.5,1.0], color='orange'
             djs_oplot, alog10(r23monte[*,ii]), res.log12oh_lower, ps=6, color='orange'               
             oploterror, x[ii], result[ii].log12oh_upper, result[ii].log12oh_upper_err, $
               psym=7, symsize=2.0, color=djs_icolor('cyan'), errcolor=djs_icolor('cyan'), $
               thick=3, errthick=3
             oploterror, x[ii], result[ii].log12oh_lower, result[ii].log12oh_lower_err, $
               psym=7, symsize=2.0, color=djs_icolor('cyan'), errcolor=djs_icolor('cyan'), $
               thick=3, errthick=3
             oploterror, x[ii], result[ii].log12oh_avg, result[ii].log12oh_avg_err, $
               psym=7, symsize=2.0, color=djs_icolor('cyan'), errcolor=djs_icolor('cyan'), $
               thick=3, errthick=3
; models          
             model_logq = alog10([1.2D8,8D7,4D7,2D7,1D7]) ; alog10(4E7)
             model_logr23 = findgen((1.1-(-0.5))/0.001+1)*0.001-0.5
             linestyle = [0,1,2,3,4,5]
             for iq = 0L, n_elements(model_logq)-1L do begin
                model_logoh_upper = 9.72D - 0.777*model_logr23 - 0.951*model_logr23^2 - 0.072*model_logr23^3 - $
                  0.811*model_logr23^4 - model_logq[iq]*(0.0737 - 0.0713*model_logr23 - 0.141*model_logr23^2 + $
                  0.0373*model_logr23^3 - 0.058*model_logr23^4)
                model_logoh_lower = 9.40D + 4.65D*model_logr23 - 3.17D*model_logr23^2 - model_logq[iq]*$
                  (0.272D + 0.547D*model_logr23 - 0.513D*model_logr23^2)
                model_good1 = where((model_logoh_upper gt model_logoh_lower))
                model_good2 = where((model_logoh_lower[model_good1] gt 7.5))
                djs_oplot, model_logr23[model_good1], model_logoh_upper[model_good1], linestyle=linestyle[iq]
                djs_oplot, model_logr23[model_good1], model_logoh_lower[model_good1], linestyle=linestyle[iq]
             endfor
             cc = get_kbrd(1)
          endif
       endfor ; close object loop
    endif 
       
    allresult[these] = result
    if (nmonte gt 0L) and arg_present(allresult_monte) then $
      allresult_monte[these] = result_monte

return, allresult
end

