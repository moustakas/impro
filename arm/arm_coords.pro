;+
; NAME: ARM_COORDS
;       
; CATEGORY: astronomy
;
; PURPOSE: convert RA/DEC coordinates to/from decimal format
;
; CALLING SEQUENCE: result = ARM_COORDS(coord, [delimiter=,
;                                       /to_hours, /from_hours)  
;
; INPUTS:
;  coord  - coordinate array of right ascensions or declinations
;           (format can be decimal or colon/space delimited strings)
;       
; OPTIONAL INPUTS:
;   delimiter - character string delimiter of coordinate components
;               (colon by defult)
;   
; KEYWORDS:
;   to_hours   - input coordinates are in hours rather than degrees
;   from_hours - output coordinates in hours rather than degrees
;
; OUTPUTS: converted coordinate returned
; 
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;
; PROCEDURES USED: ARM_DOUBLE()
;
; COMMENTS:  This was written to be more general than the standard
;            RADEC routine.  The input format is automatically
;            determined (either a character delimited string or
;            decimal) and the alternate format is returned.
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2004 June 8
;    ARM_DOUBLE() called to prevent numerical errors; ARM, 2005 Jul 26
;    additional numerical error checking; ARM, 2005 Aug 10
;
; Copyright (C) 2004, 2005, Andrew R. Marble
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

function arm_coords, coord, delimiter=delimiter, to_hours=to_hours, $
            from_hours=from_hours

    if N_ELEMENTS(delimiter) eq 0 then delimiter = ':'

    if KEYWORD_SET(to_hours) and KEYWORD_SET(from_hours) then begin
       to_hours = 0
       from_hours = 0
    endif

    n = N_ELEMENTS(coord)

    datatype = DATATYPE(coord)

    if datatype ne 'STR' then begin

       if KEYWORD_SET(from_hours) then factor = 15d0 else factor = 1d0
       if KEYWORD_SET(to_hours) then factor = factor / 15d0
       
       copy = ARM_DOUBLE(coord)
       coord2 = STRARR(n)
       
       sign = STRARR(n)
       negative = WHERE(copy lt 0, count)
       if count gt 0 then sign[negative] = '-'

       d  = FLOOR(ABS(copy*factor))
       m  = FLOOR((ABS(copy*factor) - d) * 6d1)
       s1 = FLOOR(((ABS(copy*factor) - d) * 6d1 - m) * 6d1)
       s2 = ROUND(1d2 * ((((ABS(copy*factor) - d) * 6d1 - m) * 6d1) - s1))

 ; check for precision errors

       s2bad = WHERE(s2 eq 100L, nbad)
       if nbad gt 0L then begin
          s2[s2bad] = 0L
          s1[s2bad] = s1[s2bad] + 1L
          s1bad = WHERE(s1[s2bad] eq 60L, nbad)
          if nbad gt 0L then begin
             (s1[s2bad])[s1bad] = 0L
             (m[s2bad])[s1bad] = (m[s2bad])[s1bad] + 1L
             mbad = WHERE((m[s2bad])[s1bad] eq 60L, nbad)
             if nbad gt 0L then begin
                ((m[s2bad])[s1bad])[mbad] = 0L
                ((d[s2bad])[s1bad])[mbad] = ((d[s2bad])[s1bad])[mbad] + 1L
             endif
          endif
       endif

       coord2 = sign + STRING(d,f='(i2.2)') + delimiter + STRING(m,f='(i2.2)') + $
         delimiter + STRING(s1,f='(i2.2)')+'.'+STRING(s2,f='(i2.2)')

    endif else begin
       
       coord2 = DBLARR(n)

       for i = 0L, n-1L do begin

          parts = STRSPLIT(coord[i], delimiter, /extract)
          if STRMID(STRTRIM(parts[0], 2), 0, 1) eq '-' then sign = -1d0 else sign = 1d0
          coord2[i] = ABS(parts[0]) + parts[1] / 6d1
          if N_ELEMENTS(parts) gt 2 then coord2[i] = coord2[i] + parts[2] / 36d2
          coord2[i] = sign * coord2[i]

       endfor

       if KEYWORD_SET(from_hours) then coord2 = coord2 * 15d0
       if KEYWORD_SET(to_hours) then coord2 = coord2 / 15d0

    endelse

    return, coord2

end
