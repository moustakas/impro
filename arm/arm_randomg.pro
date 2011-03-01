;+
; NAME: ARM_RANDOMG
;       
; CATEGORY: statistics
;
; PURPOSE: return a sample drawn from a parent population
;
; CALLING SEQUENCE: sample = ARM_RANDOMG(func_name, n, min, max)
;
; INPUTS:
;   func_name - string name of IDL function to evaluate for a given
;               value of x, the returned value should reflect the 
;               relative probability of drawing x from the range of
;               possible values
;   min       - minimum of range of possible x values
;   max       - maximum of range of possible x values
;   n         - number of values in returned sample
;       
; EXAMPLE: create a sample of 1000 values normally distributed about
;          zero (sigma=1), ranging between -2 and 5 
;
;          function GAUSS, x
;            return, exp(-1*x^2 / (2*1^2))
;          end 
;
;          sample = ARM_RANDOMG('gauss', -2, 5, 1d3)
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., November 2003
;-

function ARM_RANDOMG, func_name, min, max, n
    
    vals = DBLARR(n)            ; array of selected values
    cntr = 0L                   ; counter

    while cntr lt n do begin

; randomly select a value withing permitted range

       x = RANDOMU(seed, 1) * (max-min) + min

; evaluate probability of selecting that value

       p = CALL_FUNCTION(func_name, x)

; if a second number randomly selected between 0 and 1 is less than
; that probability then keep the value, otherwise discard

       rnd = RANDOMU(seed, 1)
       if rnd[0] le p then begin
          vals[cntr] = x
          cntr = cntr + 1
       endif

    endwhile

    return, vals

 end
