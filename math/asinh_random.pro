;+
; NAME:
;   ASINH_RANDOM()
;
; PURPOSE:
;   Draw a random variate from an arc sine distribution.
;
; INPUTS: 
;   x1, x2 - minimum, maximum values of the distribution 
;   num - number of variates desired
;
; OPTIONAL INPUTS: 
;   soft - softening parameter; describes the position at which the
;     asinh() distribution transitions from linear to logarithmic
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   out - random variates [NUM]
;
; COMMENTS:
;   Deals with multidimensional arrays, although it's pretty
;   slow.  Could use better error checking.
; 
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Mar 08, UCSD
;
; Copyright (C) 2011, John Moustakas
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


function asinh_random, x1, x2, num, soft=soft

    if (n_params() ne 3) then begin
       doc_library, 'asinh_random'
       return, -1
    endif
    
    if (n_elements(soft) gt 0) then ss = soft else ss = 0.3D

    out = make_array(num,type=size(x1,/type))
    count = 0
    while (count lt cmproduct(num)) do begin
       delvarx, seed1, seed2
       xx = randomu(seed1,1)
       sh = sinh(xx[0]/ss)
       prob = x1 + (x2-x1)*sh/(sinh(1.0/ss))

; if a second number randomly selected between 0 and 1 is less than
; that probability then keep the value, otherwise discard
       rnd = randomu(seed2,1)
       if (rnd[0] le prob) then begin
          out[count] = prob
          count++
       endif
    endwhile

return, out
end
