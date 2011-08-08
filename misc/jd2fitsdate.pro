;+
; NAME:
;   JD2FITSDATE()
;
; PURPOSE:
;   Convert a Julian date to a FITS format date.
;
; INPUTS: 
;   jd - input Julian date
;
; KEYWORD PARAMETERS: 
;   noyear - exclude the year from the output date
;
; OUTPUTS: 
;   date - FITS-format date (YYYY-MM-DD)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 Oct 17, U of A
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

function jd2fitsdate, jd, noyear=noyear

    njd = n_elements(jd)
    if njd eq 0L then begin
       print, 'Syntax - date = jd2fitsdate(jd)'
       return, -1
    endif

    daycnv, jd, yr, mn, day, hr

    yr = string(yr,format='(I4)')
    mn = string(mn,format='(I2.2)')
    day = string(day,format='(I2.2)')

    date = strarr(njd)
    if keyword_set(noyear) then $
      for i = 0L, njd-1L do date[i] = strjoin([mn[i],day[i]],'-') else $
      for i = 0L, njd-1L do date[i] = strjoin([yr[i],mn[i],day[i]],'-')

    if njd eq 1L then date = date[0]

return, date
end
