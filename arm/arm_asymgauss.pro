;+
; NAME: arm_asymgauss
;       
; CATEGORY: math
;
; PURPOSE: evaluate an asymmetric gaussian curve for desired points 
;
; CALLING SEQUENCE: result = ARM_ASYMGAUSS(x, lpar, rpar)
;
; INPUTS:
;   x    - array of values to evaluate
;   lpar - parameters for left side of curve (in the sense of
;          ascending values of x), see COMMENTS
;   rpar - parameter for right side of curve
;
; OPTIONAL INPUTS:
;
; KEYWORDS:
;
; OUTPUTS:
; 
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;
; PROCEDURES USED: GAUSSIAN(), POLY()
;
; COMMENTS: Each side of the curve is a gaussian function plus an
;           underlying polynomial.  LPAR and RPAR should be
;           constructed such that:
;              LPAR[0]: gaussian amplitude
;              LPAR[1]: gaussian center
;              LPAR[2]: gaussian sigma
;              LPAR[3]: 0th order polynomial coefficient (floor)
;              LPAR[4]: 1st order polynomial coefficient (slope)
;              LPAR[n]: (n-3) order polynomial coefficient
;           The LPAR parameters are used for values of x smaller than
;           (LPAR[1]+RPAR[1])/2 and the RPAR parameters are used for
;           larger values.
; 
; BUG REPORT: Please report any bugs to Andrew R. Marble.
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2004 Oct 5
;-

function arm_asymgauss, x, lpar, rpar

    y = x

    x0 = (lpar[1] + rpar[1]) / 2.0

    left  = WHERE(x lt x0, nleft)
    right = WHERE(x ge x0, nright)

    if nleft  gt 0 then y[left]  = GAUSSIAN(x[left], lpar)
    for i = 4, N_ELEMENTS(lpar)-1 do y[left] = y[left] + lpar[i] * x[left] ^ (i-3)
    if nright gt 0 then y[right] = GAUSSIAN(x[right], rpar)
    for i = 4, N_ELEMENTS(rpar)-1 do y[right] = y[right] + rpar[i] * x[right] ^ (i-3)

    return, y

 end
