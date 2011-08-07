;+
; NAME:
;   IM_FITEXY
;
; PURPOSE:
;   Calls FITEXY (linear least squares with errors in both axes) with
;   iterative outlier rejection. 
;
; INPUTS:
;   x, y - input data [NPTS]
;
; OPTIONAL INPUTS:
;   xerr, yerr - corresponding uncertainties on X,Y [NPTS] 
;   nsig - sigma rejection threshold (default 3)
;   niter - number of iterations (default 10)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   yfit - best-fitting model [NPTS]
;   coeff - best-fit linear coefficients [2]
;   good - indices of points not rejected
;   reject - indices of rejected points
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Jan 19, U of A - based on POLY_ITER 
;
; Copyright (C) 2006, John Moustakas
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

pro im_fitexy_iter, x, y, xerr=xerr, yerr=yerr, nsig=nsig, niter=niter, $
  yfit=yfit, coeff=coeff, good=good, reject=reject, maxrej=maxrej

    good = bytarr(n_elements(x))+1B
    nw = n_elements(x)
    w = lindgen(nw)

    if (n_elements(niter) eq 0) then niter = 10
    if (n_elements(nsig) eq 0) then nsig = 3.0
    
    for i=1, niter do begin

       fitexy, x[w], y[w], a_intercept, b_slope, sigma_a_b, chi_sq, $
         x_sigma=xerr[w], y_sigma=yerr[w]
       coeff = [a_intercept,b_slope]
       yfit = poly(x,coeff)

       res = y[w]-yfit[w]
       sig = stddev(res)

       rej = where(abs(res) gt nsig*sig,nrej)
       if (nrej ge maxrej) then begin
          res = reverse(res[sort(res)])
          res[maxrej:nw-1L] = 0.0
       endif
       
       good[w] = good[w]*(abs(res) lt nsig*sig)
       w = where(good,nw)

    endfor

    fitexy, x[w], y[w], a_intercept, b_slope, sigma_a_b, chi_sq, $
      x_sigma=xerr[w], y_sigma=yerr[w]
    coeff = [a_intercept,b_slope]
    yfit = poly(x,coeff)

    reject = good eq 0B
    
return
end
