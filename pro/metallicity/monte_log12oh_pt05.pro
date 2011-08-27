;+
; NAME:
;   MONTE_LOG12OH_PT05
;
; PURPOSE:
;   Compute oxygen abundances using the PT05 calibration, with
;   proper Monte Carlo of the errors.
;
; INPUTS: 
;   r23, p - if you don't know what these are you probably shouldn't
;     be using this routine 
;
; OPTIONAL INPUTS: 
;   nmonte - number of Monte Carlo realizations; note that if
;     NMONTE=0 then LOG12OH_AVG as well as the errors will not be
;     computed  
;
; KEYWORD PARAMETERS: 
;   log - indicates that R23 is in the log, i.e., log(R23); note that
;     P should *always* be linear
;
; OUTPUTS: 
;   result - output data structure (mostly self-explanatory)
;
; OPTIONAL OUTPUTS:
;   result_monte - RESULT, for each Monte Carlo realization 
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Feb 06, NYU, written
;   jm08mar21nyu - allow NMONTE=0
;   jm10jul31ucsd - better Monte Carlo errors
;
; Copyright (C) 2008, 2010, John Moustakas
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

function monte_log12oh_pt05, allr231, allp, log=log, nmonte=nmonte, $
  line=line, result_monte=allresult_monte

    if (n_elements(nmonte) eq 0L) then nmonte = 250
    nallobj = n_elements(allr231)
    
    allresult = {$
;     nmonte:                0L, $
;     frac:              -999.0, $
      log12oh_avg:       -999.0, $
      log12oh_avg_err:   -999.0, $
      log12oh_upper:     -999.0, $
      log12oh_upper_err: -999.0, $
      log12oh_lower:     -999.0, $
      log12oh_lower_err: -999.0}
    allresult = replicate(allresult,nallobj)

    if (nmonte gt 0L) and arg_present(allresult_monte) then begin
       allresult_monte = {log12oh_upper: fltarr(nmonte), log12oh_lower: fltarr(nmonte)}
       allresult_monte = replicate(allresult_monte,nallobj)
    endif

; check for good measurements    
    these = where((allr231 gt -900.0) and (allp gt -900.0),nobj)
    if (nobj eq 0L) then return, allresult else begin
       r231 = allr231[these]
       p = allp[these]
       result = allresult[these]
       if (nmonte gt 0L) and arg_present(allresult_monte) then $
         result_monte = allresult_monte[these]
    endelse

    if keyword_set(log) then r23 = 10.0^r231 else r23 = r231
    
; compute the abundances on the lower and upper branches
    result.log12oh_lower = (r23 + 106.4 + 106.8*p - 3.40*p^2) / $
      (17.72 + 6.60*p + 6.95*p^2 - 0.302*r23)
    result.log12oh_upper = (r23 + 726.1 + 842.2*p + 337.5*p^2) / $
      (85.96 + 82.76*p + 43.98*p^2 + 1.793*r23)

    result.log12oh_avg = (result.log12oh_upper+result.log12oh_lower)/2.0
    
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
       pmonte = oiiimonte / (oiimonte + oiiimonte)

       for ii = 0L, nobj-1L do begin
; enforce positivity
          montepos = where((r23monte[*,ii] gt 0.0) and $
            (pmonte[*,ii] gt 0.0),nmontepos)
          if (nmontepos gt 1L) then begin

             r23monte1 = r23monte[montepos,ii]
             pmonte1 = pmonte[montepos,ii]

             logoh_lower_monte = (r23monte1 + 106.4 + 106.8*pmonte1 - 3.40*pmonte1^2) / $
               (17.72 + 6.60*pmonte1 + 6.95*pmonte1^2 - 0.302*r23monte1)
             logoh_upper_monte = (r23monte1 + 726.1 + 842.2*pmonte1 + 337.5*pmonte1^2) / $
               (85.96 + 82.76*pmonte1 + 43.98*pmonte1^2 + 1.793*r23monte1)

             result[ii].log12oh_lower_err = djsig(logoh_lower_monte)
             result[ii].log12oh_upper_err = djsig(logoh_upper_monte)

             good = where((logoh_lower_monte lt logoh_upper_monte),ngood)
             if (ngood ne 0L) then result[ii].log12oh_avg_err = $
               djsig([logoh_lower_monte[good],logoh_upper_monte[good]])

             if arg_present(allresult_monte) then begin
                result_monte[ii].log12oh_upper[montepos] = logoh_upper_monte
                result_monte[ii].log12oh_lower[montepos] = logoh_lower_monte
             endif
          endif
       endfor
    endif
       
    allresult[these] = result
    if (nmonte gt 0L) and arg_present(allresult_monte) then $
      allresult_monte[these] = result_monte
    
return, allresult
end

