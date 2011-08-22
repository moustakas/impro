;+
; NAME:
;   MONTE_LOG12OH_ZHK94
;
; PURPOSE:
;   Compute oxygen abundances using the ZHK94 calibration, with
;   proper Monte Carlo of the errors.
;
; INPUTS: 
;   r23 - if you don't know what this is you probably shouldn't be
;     using this routine 
;
; OPTIONAL INPUTS: 
;   nmonte - number of Monte Carlo realizations; note that if
;     NMONTE=0 then the errors will not be computed (default 250) 
;   line - input emission-line structure (see IM_ABUNDANCE for an
;     example); required if NMONTE>0
;
; KEYWORD PARAMETERS: 
;   log - indicates that R23 is in the log, i.e., log(R23)
;
; OUTPUTS: 
;   result - output data structure (mostly self-explanatory)
;
; OPTIONAL OUTPUTS:
;   result_monte - RESULT, for each Monte Carlo realization 
;
; COMMENTS:
;   Only upper-branch metallicity is computed.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Aug 04, UCSD
;
; Copyright (C) 2010, John Moustakas
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

function monte_log12oh_zkh94, allr231, log=log, nmonte=nmonte, $
  line=line, result_monte=allresult_monte

    if (n_elements(nmonte) eq 0L) then nmonte = 250
    nallobj = n_elements(allr231)
    
    allresult = {$
      log12oh:     -999.0, $
      log12oh_err: -999.0}
    allresult = replicate(allresult,nallobj)

    if (nmonte gt 0L) then begin
       allresult_monte = {log12oh: fltarr(nmonte)}
       allresult_monte = replicate(allresult_monte,nallobj)
    endif

; check for good measurements    
    these = where((allr231 gt -900.0),nobj)
    if (nobj eq 0L) then return, allresult else begin
       r231 = allr231[these]
       result = allresult[these]
       if (nmonte gt 0L) then result_monte = allresult_monte[these]
    endelse

    if keyword_set(log) then r23 = 10.0^r231 else r23 = r231

; compute the abundance
    coeff = [9.265D,-0.33D,-0.202D,-0.207D,-0.333D]
    for iobj = 0L, nobj-1L do result[iobj].log12oh = poly(alog10(r23[iobj]),coeff)

; use Monte Carlo to figure out the uncertainties; use the linear
; uncertainties on the individual line-fluxes to construct the Monte
; Carlo distribution of R23, since the two quantities are correlated
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

       for iobj = 0L, nobj-1L do begin
; enforce positivity
          montepos = where((r23monte[*,iobj] gt 0.0),nmontepos)
          if (nmontepos gt 1L) then begin
             for ii = 0L, nmontepos-1L do result_monte[iobj].log12oh[montepos[ii]] = $
               poly(alog10(r23monte[montepos[ii],iobj]),coeff)
             result[iobj].log12oh_err = djsig(result_monte[iobj].log12oh[montepos])
          endif 
       endfor 
    endif 
       
    allresult[these] = result
    if (nmonte gt 0L) then allresult_monte[these] = result_monte
    
return, allresult
end

