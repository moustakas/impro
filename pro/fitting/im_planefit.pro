;+
; NAME:
;   IM_PLANEFIT()
;
; PURPOSE:
;   Fit a plane to data with the option of constraining any of the
;   coefficients. 
;
; CALLING SEQUENCE:
;   coeff = im_planefit(x,y,z,zerr=,coeff_guess=,coeff_fixed,$
;      coeff_err=,zfit=,chi2=)
;
; INPUTS:
;   x - data along the x axis [NPTS]
;   y - data along the y axis [NPTS]
;   z - data along the z axis [NPTS]
;
; OPTIONAL INPUTS:
;   zerr        - uncertainty in Z [NPTS]
;   coeff_guess - coefficients to use as the best-guess starting
;                 values [3-element vector]
;   coeff_fixed - coefficients to fix (see MPFITFUN) [3-element
;                 byte array]; (0B=free, 1B=fixed)
;
; KEYWORD PARAMETERS:
;   None.
;
; OUTPUTS:
;   coeff - coefficients of the best-fitting plane; the model is 
;           Z = COEFF[0] + COEFF[1]*X + COEFF[2]*Y
;
; OPTIONAL OUTPUTS:
;   coeff_err - formal one-sigma uncertainties in COEFF (only
;               meaningful if ZERR is given) [3-element vector] 
;   zfit      - best-fitting model [NPTS]
;   chi2      - reduced chi-squared statistic (only meaningful
;               if ZERR is given)
;
; PROCEDURES USED:
;   MPFITFUN()
;
; COMMENTS:
;   It would be nice to be able to pass uncertainties in X and Y
;   as well as Z.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Mar 04, U of A - written
;
; Copyright (C) 2005, John Moustakas
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

function planefunc, xy, p
    x = reform(xy[0,*])
    y = reform(xy[1,*])
    z = p[0] + p[1]*x + p[2]*y
return, z
end    

function im_planefit, x, y, z, zerr=zerr, coeff_guess=coeff_guess, $
  coeff_fixed=coeff_fixed, coeff_err=coeff_err, zfit=zfit, chi2=chi2

    nx = n_elements(x)
    ny = n_elements(y)
    nz = n_elements(z)

    if (nx eq 0L) then begin
       doc_library, 'im_planefit'
       return, -1
    endif

    if (nx ne ny) or (nx ne nz) then begin
       print, 'X, Y, and Z must have the same number of elements.'
       return, -1
    endif
    
    nzerr = n_elements(zerr)
    if (nzerr eq 0L) then zerr = z*0.0+1.0 else if (nx ne nzerr) then begin
       print, 'Z and ZERR must have the same number of elements.'
       return, -1
    endif
    
    nparams = 3
    ncoeff_guess = n_elements(coeff_guess)
    if (ncoeff_guess eq 0L) then coeff_guess = fltarr(nparams)+1.0 else $
      if (ncoeff_guess ne nparams) then begin
       print, 'COEFF_GUESS must be a 3-element array.'
       return, -1
    endif
    
    ncoeff_fixed = n_elements(coeff_fixed)
    if (ncoeff_fixed eq 0L) then coeff_fixed = bytarr(nparams)+0B else $
      if (ncoeff_fixed ne nparams) then begin
       print, 'COEFF_FIXED must be a 3-element byte array.'
       return, -1
    endif
    
; do the fit
    nparams = 3
    parinfo = replicate({value: 0.0, fixed: 0B},nparams)
    parinfo.value = coeff_guess
    parinfo.fixed = coeff_fixed

    xy = transpose([[x],[y]])
    coeff = mpfitfun('planefunc',xy,z,zerr,parinfo=parinfo,$
      perror=coeff_err,yfit=zfit,bestnorm=bestnorm,quiet=1)

; compute chi-squared
    chi2 = bestnorm/(n_elements(xy)-nparams)

return, coeff
end    
