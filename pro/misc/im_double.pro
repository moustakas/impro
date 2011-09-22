;+
; NAME: 
;   IM_DOUBLE
;       
; PURPOSE: 
;   Convert a float to double-precision avoiding roundoff issues. 
;
; CALLING SEQUENCE: 
;   result = IM_DOUBLE(x)
;
; INPUTS: 
;   x - a number
;
; OUTPUTS: 
;   returns x as double-precision number
;
; COMMENTS: 
;   If you convert a non double-precision number to double-precision
;   using DOUBLE(), you will get a numerical error (e.g., DOUBLE(0.01)
;   --> 0.0099999998).  This routine was written in order to "up the
;   precision" of a number without encountering this problem.
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., December 3, 2003
;    modified to accomodate exponential format, ARM, December 5, 2003
;    return scalar rather than array of 1 value, ARM, 2004 June 4
;    jm11apr25ucsd - ported into IMPRO
;-

function im_double, x

    nx = N_ELEMENTS(x)          ; number of elements in passed array

    xd = DBLARR(nx)       ; array of double-precision values to return

    for i=0L,nx-1L do begin

       parts = STRSPLIT(x[i], 'e', /extract) ; check for exponential format

       xnew = parts[0] + 'd0'
       if N_ELEMENTS(parts) gt 1L then xnew = xnew + '*10d0^(' + parts[1] + 'd)' 

       dummy = EXECUTE('xd[i] = ' + xnew)  ; redefine x to be double-precision
       
    endfor

    if nx eq 1L then xd = xd[0] ; return scalar rather an array of 1 value

    return, xd
    
 end
