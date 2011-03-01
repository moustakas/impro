;+
; NAME: ARM_ORDEROFMAG
;       
; CATEGORY: math
;
; PURPOSE: returns the order of magnitude of a number
;
; CALLING SEQUENCE: p = ARM_ORDEROFMAG(x)
;
; INPUTS: x - array of one or more numbers
;
; EXAMPLE: IDL> p = ARM_ORDEROFMAG(123.456)  ; p ->  2
;               p = ARM_ORDEROFMAG(-123.456) ; p ->  2
;               p = ARM_ORDEROFMAG(0.00123)  ; p -> -2
;
; ROUTINES CALLED:
;
; COMMENTS: This routine works for numbers 2.47e-324 < # < 1.80e308.
;           Returns zero when x = 0.
;          
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., September 2001
;    x = 0 bug fixed, ARM, December 3, 2003
;    return scalar value instead of array of 1 value, ARM, Jan 2 2004
;    completely rewritten (thanks to Jeremiah Murphy), ARM, Jul 6 2004
;-

function arm_orderofmag, x

; check for zeroes

   zero = WHERE(x eq 0, count)
   if count gt 0 then x[zero] = 1

   power = ALOG10(ABS(x))

; round negative powers down

   neg = WHERE(power lt 0, count)
   if count gt 0 then power[neg] = power[neg] - 1

   return, FIX(power, type=3)

 end

