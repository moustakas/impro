;+
; NAME:
;   GET_MANGLE_AREA()
;
; PURPOSE:
;   Compute the total area of a set of mangle polygons. 
;
; INPUTS: 
;   polyfile - polygon file name
;
; KEYWORD PARAMETERS: 
;   deg - convert AREA to deg^2
;
; OUTPUTS: 
;   area - total area in steradians, unless /DEG 
;
; COMMENTS:
;   Loosely based on K. Wong's PRIMUS_SURVEY_AREA.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Sep 10, UCSD
;
; Copyright (C) 2010, John Moustakas
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

function get_mangle_area, polyfile, deg=deg

    if (n_elements(polyfile) eq 0) then begin
       doc_library, 'get_mangle_area'
       return, 0D
    endif
    if (file_test(polyfile) eq 0) then begin
       splog, 'Polygon file '+polyfile+' not found!'
       return, 0D
    endif

    splog, 'Reading '+polyfile
    read_mangle_polygons, polyfile, poly
    area = total(poly.str,/double)

;   area = 0D & for ii = 0L, n_elements(poly)-1 do area += garea(poly[ii])
    if keyword_set(deg) then area *= (180D/!pi)^2

return, area
end
