;+
; NAME:
;   MONOTONIC()
;
; PURPOSE:
;   Determine whether a vector is monotonically increasing or
;   decreasing. 
;
; INPUTS:
;   x - one-dimensional vector of any datatype
;
; OUTPUTS:
;   Returns 1 if the vector is monotonically increasing or
;   decreasing, and 0 otherwise.
;
; COMMENTS:
;   Does not distinguish between increasing and decreasing
;   vectors.  Could be extended to two-dimensional vectors
;   easily. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2001 July 23, U of A, written
;   jm08jan11nyu - check for uniqueness
;   jm13sep17siena - code replaced with FINDEX() algorithm which is
;     much more stable; also optionally returns whether the array is
;     monotonically increasing or decreasing
;
; Copyright (C) 2001, 2008, 2013, John Moustakas
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

function monotonic, x, increasing=increasing

    xx = x-shift(x,+1)
    xx = xx[1:*]
    mx = max(xx,min=mn)
    if mx gt 0 and mn lt 0 then true = 0 else true = 1

    if mx gt 0 then increasing = 1 else increasing = 0
    
;; old code below!    
;    nx = n_elements(x)
;    vec = lindgen(nx)   ; comparison vector
;    srt = sort(x)       ; sort the vector       
;    rsrt = reverse(srt) ; reverse the sorted vector       
;    
;    true = 0
;    if total(abs(srt-vec)) eq float(0) then true = 1 else $ ; monotonically increasing
;      if total(abs(rsrt-vec)) eq float(0) then true = 1     ; monotonically decreasing
;
;; make sure all the values aren't the same
;    if (total(x-x[0]) eq 0.0) then true = 0L
    
return, true
end
