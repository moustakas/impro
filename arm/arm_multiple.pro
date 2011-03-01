;+
; NAME: ARM_MULTIPLE
;       
; CATEGORY: math
;
; PURPOSE: return 'true' if 2nd number is a multiple of the 1st
;
; CALLING SEQUENCE: result = ARM_MULTIPLE(x1, x2)
;
; INPUTS:
;   x1 - integer
;   x2 - array of one or more "candidate multiples" of x1
;       
; OUTPUTS: returns array with ones indicating multiples
;
; COMMENTS: numbers must be integers
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2001
;-

function arm_multiple, x1, x2

    if x1-FIX(x1) ne 0L or TOTAL(x2-FIX(x2)) ne 0L then $
      MESSAGE, 'numbers must be integers!'
    
    x1 = FIX(x1)
    x2 = FIX(x2)

    n = N_ELEMENTS(x2)

    result = INTARR(n)

    zero = WHERE(x2 eq 0L, zerocount)
    okay = WHERE(x2 ne 0L, okaycount)

    if zerocount gt 0L then result[zero] = 1L

    if okaycount gt 0L then begin
       multiple = WHERE(x2[okay]/x1 ne (x2[okay]-1)/x1, count)
       if count gt 0L then result[okay[multiple]] = 1L
    endif

    return, result

 end
 
