;+
; NAME:
;   GET_ELEMENT
;
; PURPOSE:
;   Compute the index of a value, given an input array.
;
; INPUTS: 
;   x - input vector
;   value - input value (can be an array)
;
; OUTPUTS: 
;   position - index or indices of VALUE in X
;
; COMMENTS:
;   This routine is pretty silly.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2001 Jan 28, U of A
;
; Copyright (C) 2001, John Moustakas
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

pro get_element, x, value, position
    position = long(value-value)
    for i = 0L, n_elements(value)-1L do begin
       array_value = min((abs(x-value[i])),temp)
       position[i] = temp
    endfor
return
end
