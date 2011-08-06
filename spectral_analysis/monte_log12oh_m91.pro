;+
; NAME:
;   MONTE_LOG12OH_M91
;
; PURPOSE:
;   Compute oxygen abundances using the M91 calibration, with
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

function monte_log12oh_m91, allr231, allo321, log=log, nmonte=nmonte, $
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
    result.log12oh_lower = (12.0 - 4.944 + 0.767*x + 0.602*x^2) - y*(0.29 + 0.332*x - 0.331*x^2)
    result.log12oh_upper = (12.0 - 2.939 - 0.2*x - 0.237*x^2 - 0.305*x^3 - 0.0283*x^4) - $
      y*(0.0047 - 0.0221*x - 0.102*x^2 - 0.0817*x^3 - 0.00717*x^4)

; compute the average
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
       o32monte = oiiimonte / oiimonte

       for ii = 0L, nobj-1L do begin
; enforce positivity
          montepos = where((r23monte[*,ii] gt 0.0) and $
            (o32monte[*,ii] gt 0.0),nmontepos)
;         if (nmontepos ne nmonte) then message, 'Problem here'
          if (nmontepos gt 1L) then begin

             xmonte = alog10(r23monte[montepos,ii])
             ymonte = alog10(o32monte[montepos,ii])

             logoh_lower_monte = 12.0 - 4.944 + 0.767*xmonte + 0.602*xmonte^2 - $
               ymonte*(0.29 + 0.332*xmonte - 0.331*xmonte^2)
             logoh_upper_monte = 12.0 - 2.939 - 0.2*xmonte - 0.237*xmonte^2 - $
               0.305*xmonte^3 - 0.0283*xmonte^4 - $
               ymonte*(0.0047 - 0.0221*xmonte - 0.102*xmonte^2 - 0.0817*xmonte^3 - $
               0.00717*xmonte^4)

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

