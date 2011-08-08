;+
; NAME:
;   FIND_DIFFS()
;
; PURPOSE:
;   Find all differences of a vector X, or between two vectors X and
;   Y. 
;
; INPUTS: 
;   x - input vector
;
; OPTIONAL INPUTS: 
;   y - optional second input vector
;
; OUTPUTS: 
;   diffs - differences within X or between X and Y
;
; COMMENTS:
;   Slow and dumb; should use HISTOGRAM().
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2001 Jul 23, U of A
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

function find_diffs, x, y=y

    nx = n_elements(x)

    if keyword_set(y) then begin

       ny = n_elements(y)
       ndiffs = nx*ny
       diffs = fltarr(ndiffs)

       count = 0L
       for i = 0L, nx-1L do begin

          for j = 0L, ny-1L do begin
             diffs[count] = abs(x[i]-y[j])
             count = count+1L
          endfor
       
       endfor

    endif else begin
    
       ndiffs = nx*(nx-1L)/2L
       diffs = fltarr(ndiffs)

       count = 0L
       for i = 0L, nx-2L do begin

          for j = i+1L, nx-1L do begin
             diffs[count] = abs(x[i]-x[j])
             count = count+1L
          endfor
       
       endfor

    endelse

return, diffs
end
