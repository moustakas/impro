;+
; NAME:
;       IM_HYPERFIT()
;
; PURPOSE:
;       Fit a hyperbola to data with the option of constraining any of
;       the coefficients. 
;
; CALLING SEQUENCE:
;       coeff = im_hyperfit(x,y,yerr=,coeff_guess=,coeff_fixed,$
;          coeff_limits=,coeff_limited=,coeff_err=,yfit=,chi2=)
;
; INPUTS:
;       x - data along the x axis [NPTS]
;       y - data along the y axis [NPTS]
;
; OPTIONAL INPUTS:
;       yerr          - uncertainty in Y [NPTS]
;       coeff_guess   - coefficients to use as the best-guess starting
;                       values [2-element vector]
;       coeff_fixed   - coefficients to fix (see MPFITFUN) [2-element
;                       byte array]; (0B=free, 1B=fixed)
;       coeff_limits  - coefficients limits (see MPFITFUN)
;                       [2x2-element float array]
;       coeff_limited - see MPFITFUN [2x2-element byte array]
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       coeff - coefficients of the best-fitting hyperbola; the model
;               is  Y = COEFF[0] / (X + COEFF[1]) + COEFF[2]
;
; OPTIONAL OUTPUTS:
;       coeff_err - formal one-sigma uncertainties in COEFF (only
;                   meaningful if YERR is given) [2-element vector] 
;       yfit      - best-fitting model [NPTS]
;       chi2      - reduced chi-squared statistic (only meaningful
;                   if YERR is given)
;
; PROCEDURES USED:
;       MPFITFUN()
;
; COMMENTS:
;       It would be nice to be able to pass uncertainties in X as well
;       as Y. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 July 18, U of A - written
;       jm03aug08uofa - added COEFF_GUESS optional input 
;       jm05mar05uofa - documented and cleaned up
;       jm06feb11uofa - added COEFF_LIMITED and COEFF_LIMITS inputs 
;
; Copyright (C) 2002-2003; 2005-2006, John Moustakas
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

function hyperfunc, x, p
    y = p[0] / (x + p[1]) + p[2]
return, y
end    

function im_hyperfit, x, y, yerr=yerr, coeff_guess=coeff_guess, $
  coeff_fixed=coeff_fixed, coeff_limits=coeff_limits, coeff_limited=coeff_limited, $
  coeff_err=coeff_err, yfit=yfit, chi2=chi2

    nx = n_elements(x)
    ny = n_elements(y)

    if (nx eq 0L) then begin
       print, 'Syntax - coeff = im_hyperfit(x,y,yerr=,coeff_guess=,$'
       print, '   coeff_fixed,coeff_err=,yfit=,chi2=)'
       return, -1L
    endif

    if (nx ne ny) then begin
       print, 'X and Y must have the same number of elements.'
       return, -1L
    endif
    
    nyerr = n_elements(yerr)
    if (nyerr eq 0L) then yerr = y*0.0+1.0 else if (nx ne nyerr) then begin
       print, 'Y and YERR must have the same number of elements.'
       return, -1L
    endif
    
    nparams = 3L

    ncoeff_guess = n_elements(coeff_guess)
    if (ncoeff_guess eq 0L) then coeff_guess = fltarr(nparams)+1.0 else $
      if (ncoeff_guess ne nparams) then begin
       print, 'COEFF_GUESS must be a 3-element array.'
       return, -1L
    endif
    
    ncoeff_fixed = n_elements(coeff_fixed)
    if (ncoeff_fixed eq 0L) then coeff_fixed = bytarr(nparams)+0B else $
      if (ncoeff_fixed ne nparams) then begin
       print, 'COEFF_FIXED must be a 3-element byte array.'
       return, -1L
    endif

; add this error checking later

    ncoeff_limits = size(coeff_limits,/dim)
    if (ncoeff_limits[0] eq 0L) then coeff_limits = fltarr(2,nparams)+0.0 else $
      if (ncoeff_limits[0] ne 2L) and (ncoeff_limits[1] ne nparams) then begin
       print, 'COEFF_LIMITS must be a 3x3-element array.'
       return, -1L
    endif
    
    ncoeff_limited = size(coeff_limited,/dim)
    if (ncoeff_limited[0] eq 0L) then coeff_limited = bytarr(2,nparams)+0B else $
      if (ncoeff_limited[0] ne 2L) and (ncoeff_limited[1] ne nparams) then begin
       print, 'COEFF_LIMITED must be a 3x3-element byte array.'
       return, -1L
    endif
    
; do the fit
    
    parinfo = replicate({value: 0.0, fixed: 0B, limits: [0.0,0.0], limited: [0,0]},nparams)
    parinfo.value = coeff_guess
    parinfo.fixed = coeff_fixed
    parinfo.limits = coeff_limits
    parinfo.limited = coeff_limited

    coeff = mpfitfun('hyperfunc',x,y,yerr,parinfo=parinfo,$
      perror=coeff_err,yfit=yfit,bestnorm=bestnorm,quiet=1,$
      status=status)
    splog, 'Status = ', status

; compute chi-squared
    
    chi2 = bestnorm/(n_elements(x)-nparams)

return, coeff
end    
