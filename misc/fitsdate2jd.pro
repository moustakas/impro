;+
; NAME:
;   FITSDATE2JD()
;
; PURPOSE:
;   Convert FITS format date to Julian date.
;
; INPUTS: 
;   date - FITS-format date (YYYY-MM-DD)
;
; OUTPUTS: 
;   jd - input Julian date
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 Apr 23, U of A
;
; Copyright (C) 2002, John Moustakas
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

function fitsdate2jd, date

    ndate = n_elements(date)
    if ndate eq 0L then begin
       print, 'Syntax - date = fitsdate2jd(date)'
       return, -1
    endif

    jd = dblarr(ndate)
    for i = 0L, ndate-1L do begin
    
       d = strsplit(date[i],'-',/extract)
       jd[i] = julday(d[1],d[2],d[0])

    endfor

    if ndate eq 1L then jd = jd[0]

return, jd
end
