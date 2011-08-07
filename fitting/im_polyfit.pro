;+
; NAME:
;       IM_POLYFIT()
;
; PURPOSE:
;       Fit a polynomial to data with the option of constraining any
;       of the coefficients. 
;
; CALLING SEQUENCE:
;       coeff = im_polyfit(x,y,order,yerr=,coeff_guess=,$
;          coeff_fixed,coeff_limited,coeff_limits=,$
;          coeff_err=,yfit=,chi2=)
;
; INPUTS:
;       x     - data along the x axis [NPTS]
;       y     - data along the y axis [NPTS]
;       order - order of the polynomial
;
; OPTIONAL INPUTS:
;       yerr          - uncertainty in Y [NPTS] 
;       coeff_guess   - coefficients to use as the best-guess starting 
;                       values [(ORDER+1)-element vector]
;       coeff_fixed   - coefficients to fix (see MPFITFUN)
;                       [(ORDER+1)-element byte array]; (0B=free, 
;                       1B=fixed) 
;       coeff_limited - coefficients limited (see MPFITFUN)
;       coeff_limits  - coefficient limits (see MPFITFUN)
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       coeff - coefficients of the best-fitting polynomial 
;
; OPTIONAL OUTPUTS:
;       coeff_err - formal one-sigma uncertainties in COEFF (only
;                   meaningful if YERR is given) [(ORDER+1)-element
;                   vector]  
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
;       J. Moustakas, 2004 March 24, U of A, written, based entirely
;          on IM_LINEFIT()
;       jm05mar24uofa - added COEFF_LIMITED and COEFF_LIMITS inputs 
;
; Copyright (C) 2004-2005, John Moustakas 
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

function polyfunc, x, p

    model = p[0]
    order = n_elements(p)-1L
    for i = 1L, order do model = model + p[i]*x^i

return, model
end    

function im_polyfit, x, y, order, yerr=yerr, coeff_guess=coeff_guess, $
  coeff_fixed=coeff_fixed, coeff_limited=coeff_limited, coeff_limits=coeff_limits, $
  coeff_err=coeff_err, yfit=yfit, chi2=chi2

    nx = n_elements(x)
    ny = n_elements(y)
    norder = n_elements(order)

    if (nx eq 0L) or (ny eq 0L) or (norder eq 0L) then begin
       print, 'Syntax - coeff = im_polyfit(x,y,order,yerr=,coeff_guess=,$'
       print, '   coeff_fixed,coeff_limited=,coeff_limits=,coeff_err=,$'
       print, '   yfit=,chi2=)'
       return, -1L
    endif

    if (nx ne ny) then begin
       print, 'X and Y must have the same number of elements.'
       return, -1L
    endif

    if (order lt 1L) then begin
       print, 'ORDER must be greater than or equal to one.'
       return, -1L
    endif
    
    nyerr = n_elements(yerr)
    if (nyerr eq 0L) then yerr = y*0.0+1.0 else if (nx ne nyerr) then begin
       print, 'Y and YERR must have the same number of elements.'
       return, -1L
    endif
    
    nparams = order+1L ; number of coefficients/parameters

    ncoeff_guess = n_elements(coeff_guess)
    if (ncoeff_guess eq 0L) then coeff_guess = fltarr(nparams)+1.0 else $
      if (ncoeff_guess ne nparams) then begin
       print, 'COEFF_GUESS must be an (ORDER+1)-element array.'
       return, -1L
    endif
    
    ncoeff_fixed = n_elements(coeff_fixed)
    if (ncoeff_fixed eq 0L) then coeff_fixed = bytarr(nparams)+0B else $
      if (ncoeff_fixed ne nparams) then begin
       print, 'COEFF_FIXED must be an (ORDER+1)-element byte array.'
       return, -1L
    endif
    
    ncoeff_limited = n_elements(coeff_limited)
    if (ncoeff_limited eq 0L) then coeff_limited = bytarr(2,nparams)+0B else $
      if (ncoeff_limited ne nparams) then begin
       print, 'COEFF_LIMITED must be a 2 by (ORDER+1)-element byte array.' 
       return, -1L
    endif
    
    ncoeff_limits = n_elements(coeff_limits)
    if (ncoeff_limits eq 0L) then coeff_limits = fltarr(2,nparams) else $
      if (ncoeff_limits ne nparams) then begin
       print, 'COEFF_LIMITS must be an 2 by (ORDER+1)-element floating-point array.'
       return, -1L
    endif
    
; do the fit
    
    parinfo = replicate({value: 0.0, fixed: 0B, limited: [0L,0L], limits: [0.0,0.0]},nparams)
    parinfo.value = coeff_guess
    parinfo.fixed = coeff_fixed
    parinfo.limited = coeff_limited
    parinfo.limits = coeff_limits

    coeff = mpfitfun('polyfunc',x,y,yerr,parinfo=parinfo,$
      perror=coeff_err,yfit=yfit,bestnorm=bestnorm,quiet=1,$
      covar=covar)

; compute chi-squared
    
    chi2 = bestnorm/(n_elements(x)-nparams)

return, coeff
end    
