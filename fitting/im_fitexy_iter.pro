;+
; NAME:
;       IM_FITEXY
;
; PURPOSE:
;       Calls FITEXY (linear least squares with errors in both axes) 
;       with iterative outlier rejection.
;
; CALLING SEQUENCE:
;       im_fitexy_iter, x, y, xerr=, yerr=, nsig=, niter=, yfit=, coeff=
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Jan 19, U of A - based on POLY_ITER
;-

pro im_fitexy_iter, x, y, xerr=xerr, yerr=yerr, nsig=nsig, niter=niter, $
  yfit=yfit, coeff=coeff, good=good, reject=reject, maxrej=maxrej

    good = bytarr(n_elements(x))+1B
    nw = n_elements(x)
    w = lindgen(nw)

    if (n_elements(niter) eq 0L) then niter = 10L
    if (n_elements(nsig) eq 0L) then nsig = 3.0
    
    for i=1L, niter do begin

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
