;+
; NAME:
;   im_poly_iter
;
; PURPOSE:
;   Calls IDL poly_fit iteratively with outlier rejection
;
; CALLING SEQUENCE:
;   im_poly_iter, x, y, ndeg, nsig, yfit, coeff=coeff
;
; INPUTS:
;   x, y    - indep, dep variables
;   ndeg    - degree of polynomial 
;   yfit    - fit y at given x values
;
; OUTPUTS:
;   coeff   - array of coefficients
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   poly_fit
;
; REVISION HISTORY:
;   moved from hoggpt 10-Jan-2002
;   jm03mar16uofa - added GOOD and REJECT output
;   jm06jan19uofa - added YERR optional input
;   jm07jun19nyu  - if all the points have been rejected then return
;                   without crashing; also the rejection criterion has
;                   been changed to take into account machine
;                   precision (e.g., if you try to fit an array
;                   of 1's) 
;
;-

PRO im_poly_iter, x, y, ndeg, yfit, coeff=coeff, good=good, reject=reject, $
  yerr=yerr, niter=niter, nsig=nsig, e_coeff=e_coeff, sigma=sigma

    good = bytarr(n_elements(x))+1B
    w = lindgen(n_elements(x))

    if (n_elements(niter) eq 0L) then niter = 10L
    if (n_elements(nsig) eq 0L) then nsig = 3.0
    if (n_elements(yerr) eq 0L) then doerror = 0L else doerror = 1L

    for i=1L, niter do begin

       if doerror then $
         coeff = poly_fit(x[w],y[w],ndeg,yfit,measure_errors=yerr[w],/double,sigma=e_coeff) else $
         coeff = poly_fit(x[w],y[w],ndeg,yfit,/double,sigma=e_coeff)

       res = y[w]-yfit
       sig = stddev(res)
       if (sig gt 0.0) then begin ; special case
          good[w] = good[w]*(abs(res) LT nsig*sig)
          wnew = where(good,nw)
          if (nw eq 0L) then splog, 'All points rejected!' else $ ; keep the old "w" array
            w = wnew
       endif else continue

    ENDFOR 

    if doerror then $
      coeff = poly_fit(x[w],y[w],ndeg,yfit,measure_errors=yerr[w],/double,sigma=e_coeff) else $
      coeff = poly_fit(x[w],y[w],ndeg,yfit,/double,sigma=e_coeff)

    yfit = poly(x,coeff)
    sigma = djsig(y-yfit)

    reject = good eq 0B
    
return
end
