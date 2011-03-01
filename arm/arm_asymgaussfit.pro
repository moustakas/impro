;+
; NAME: arm_asymgaussfit
;       
; CATEGORY: math
;
; PURPOSE: fit an asymmetric gaussian+polynomial to data
;
; CALLING SEQUENCE: result = arm_asymgaussfit(x, y, [error, mask=,
;                     lterms=, rterms=, lpars=, rpars=, chisq=, /plot,
;                     /demo, _extra=] 
;
; INPUTS:
;   x - array of independent values
;   y - data points corresponding to X
;       
; OPTIONAL INPUTS:
;   error  - 1-sigma uncertainties in Y (fit will be weighted by
;            1/ERROR^2, if ERROR if not supplied then weights=1)
;   mask   - array of mask values corresponding to X (zeroes indicate
;            values of X to ignore in the fittin process)
;   lterms - number of parameters ([gaussian amplitude, gaussian
;            center, gaussian sigma, 0th order polynomial coefficient,
;            1st order polynomial coefficient, ..., nth order
;            polynomial coefficient]) to fit to left side of curve (in
;            the sense of ascending values of x)
;   rterms - number of parameters to fit to right side of curve
;   _extra - additional parameters/keywords to pass to MPFITFUN()
;
; KEYWORDS:
;   plot - display data, fit and fit components
;   demo - run ARM_ASYMGAUSSFIT in demonstration mode
;
; OUTPUTS: 
;   returns fit values evaluated at X
; 
; OPTIONAL OUTPUTS:
;   lpars - fitted parameters for left side of curve
;   rpars - fitted parameters for right side of curve
;   chisq - Chi squared statistic (BESTNORM parameter returned by
;           MPFITFUN.PRO) 
;
; EXAMPLE:
;
; PROCEDURES USED: MPFITFUN(), MPFITPEAK(), ARM_GETINDEX(),
;                  ARM_ASYMGAUSS(), DJS_OPLOT
;
; COMMENTS:
; 
; BUG REPORT: Please report any bugs to Andrew R. Marble.
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2004 Oct 12
;    display given higher resolution than data, ARM, 2004 Oct 13
;    CHISQ variable added, ARM, 2004 Oct 13
;-

function arm_asymgaussfitfunc, x, p

    n = N_ELEMENTS(p)
    return, ARM_ASYMGAUSS(x, p[0:n/2-1], p[n/2:n-1])

end

pro arm_asymgaussfitdemo

    PRINT
    PRINT, ' [demo] IDL> x = FINDGEN(201)/10'
    PRINT, ' [demo] IDL> y = ARM_ASYMGAUSS(x, [62,13,3,25,1], [90,13,1,10])'
    PRINT, ' [demo] IDL> result = ARM_ASYMGAUSSFIT(x, y, lterms=5, rterms=4, $'
    PRINT, ' [demo]        /quiet, /plot, lpars=lp, rpars=rp)'

    lp = [62, 13, 3, 25, 1]
    rp = [90, 13, 1, 10]
    x = FINDGEN(201)/10
    y = ARM_ASYMGAUSS(x, lp, rp)
    dummy = ARM_ASYMGAUSSFIT(x, y, lterms=5, rterms=4, /quiet, /plot)

    PRINT, ' [demo] IDL> PRINT, lp, rp'
    PRINT, ' [demo]       62.0000      13.0000      3.00000      25.0000      1.00000'
    PRINT, ' [demo]       90.0000      13.0000      1.00000      10.0000'
    PRINT

 end
 
function arm_asymgaussfit, x, y, error, mask=mask, lterms=lterms, rterms=rterms, chisq=chisq, $
                           lpars=lpars, rpars=rpars, plot=plot, _extra=extra, demo=demo

    if KEYWORD_SET(demo) then begin
       ARM_ASYMGAUSSFITDEMO
       return, -1
    endif

; defaults and error-checking

    if N_ELEMENTS(mask) eq 0L then mask = x * 0 + 1

    if N_ELEMENTS(error) eq 0L then weights = x * 0 + 1 else weights = 1d0 / (error^2d)
    
    n = N_ELEMENTS(x)
    if N_ELEMENTS(y) ne n or N_ELEMENTS(weights) ne n or N_ELEMENTS(mask) ne n then $
      MESSAGE, 'X, Y, MASK and WEIGHTS have incompatible dimensions.'

    if N_ELEMENTS(lterms) eq 0L then lterms = 3
    if N_ELEMENTS(rterms) eq 0L then rterms = 3

    lterms = ROUND(lterms) > 3         ; minimal gaussian parameters
    rterms = ROUND(rterms) > 3
    nterms = lterms > rterms

; ignore masked values

    xall = x
    yall = y
    masked = WHERE(mask eq 0, nmasked)
    if nmasked gt 0 then REMOVE, masked, x, y, weights

; fit gaussians to determine initial value guesses

    x0 = (x[ARM_GETINDEX(y, MAX(y))])[0]

    left = WHERE(x lt x0, nleft)
    if nleft eq 0 then lfit = REPLICATE(0,lterms) else $
      lgauss = MPFITPEAK([x[left],2*x0-REVERSE(x[left])], [y[left],REVERSE(y[left])], lguess, nterms=lterms<4)

    right = WHERE(x gt x0, nright)
    if nright eq 0 then rfit = REPLICATE(0,rterms) else $
      rgauss = MPFITPEAK([REVERSE(2*x0-x[right]),x[right]], [REVERSE(y[right]),y[right]], rguess, nterms=rterms<4)

    value = DBLARR(2*nterms)
    value[0:(lterms-1)<3] = lguess
    value[nterms:nterms+((rterms-1)<3)] = rguess

; fix extraneous parameters to be zero

    fixed = INTARR(2*nterms) + 1
    fixed[0:lterms-1] = 0
    fixed[nterms:nterms+rterms-1] = 0

; constrain left and right gaussian centers to be equal

    tied = STRARR(2*nterms)
    tied[nterms+1] = 'P[1]'

; constrain left and right curves be equal where they meet

    tied[nterms] = 'P[0]'
    for i = 3,lterms-1 do tied[nterms] = $
      tied[nterms] + '+P['+STRN(i,f="(i)")+']*P[1]^'+STRN(i-3,f="(i)")
    for i = 3,rterms-1 do tied[nterms] = $
      tied[nterms] + '-P['+STRN(nterms+i,f="(i)")+']*P['+STRN(nterms+1,f="(i)")+']^'+STRN(i-3,f="(i)")

    parinfo = {value: value, fixed: fixed, tied: tied}

    result = MPFITFUN('arm_asymgaussfitfunc', x, y, weights=weights, parinfo=parinfo, $
                      bestnorm=bestnorm, _extra=extra)

    dof = TOTAL(mask) - (lterms + rterms -2) ; degrees of freedom
    chisq = SQRT(bestnorm / dof)

    lpars = result[0:lterms-1]
    rpars = result[nterms:nterms+rterms-1]

    xhighres = DINDGEN(10*n+1)/(10*n)*(MAX(xall)-MIN(xall))+MIN(xall)

    x = xall
    y = yall

    yfit = ARM_ASYMGAUSS(x, lpars, rpars)
    yfitall = ARM_ASYMGAUSS(xhighres, lpars, rpars)

    if KEYWORD_SET(plot) then begin

       PLOT, xall, yall, xst=3, yst=3, /nodata

       left = WHERE(xhighres lt result[1], nleft)
       if nleft gt 0 then begin
          DJS_OPLOT, xhighres[left], GAUSSIAN(xhighres[left], lguess), line=2, color='green'
          DJS_OPLOT, xhighres[left], GAUSSIAN(xhighres[left], lpars[0:2]), color='green'
          if lterms gt 3 then begin
             if lterms gt 4 then DJS_OPLOT, xhighres[left], POLY(xhighres[left], lpars[3:lterms-1]), color='green' $ 
             else DJS_OPLOT, xhighres[left], left*0+lpars[3], color='green'
          endif
       endif

       right = WHERE(xhighres gt result[1], nright)
       if nright gt 0 then begin
          DJS_OPLOT, xhighres[right], GAUSSIAN(xhighres[right], rguess), line=2, color='blue'
          DJS_OPLOT, xhighres[right], GAUSSIAN(xhighres[right], rpars[0:2]), color='blue'
          if rterms gt 3 then begin
             if rterms gt 4 then DJS_OPLOT, xhighres[right], POLY(xhighres[right], rpars[3:rterms-1]), color='blue' $
             else DJS_OPLOT, xhighres[right], right*0+rpars[3], color='blue'
          endif
       endif

       DJS_OPLOT, x, y, psym=2, symsize=0.5, symthick=2
       DJS_OPLOT, xhighres, yfitall, color='red', thick=2

    endif
    
    return, yfit
    
 end
