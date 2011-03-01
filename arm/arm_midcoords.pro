;+
; NAME: ARM_MIDCOORDS
;       
; CATEGORY: astronomy
;
; PURPOSE: compute center for paired coordinates and vice-versa
;
; CALLING SEQUENCE: result = ARM_MIDCOORDS(ra, dec, [ra2, dec2, pa=,
;                            sep=, /arcmin, /arcsec, /silent, /other]) 
;
; INPUTS:
;  ra  - right ascension coordinate array for center or first of
;        paired positions (can be in degrees, 'hh:mm:ss.ss', 'hh:mm',  
;        'hh mm ss.ss' or 'hh mm' format)
;  dec - declination coordinate array for center or first of paired
;        posittions (can be in degrees, 'dd:mm:ss.ss', 'dd:mm', 
;        'dd mm ss.ss' or 'dd mm' format) 
;       
; OPTIONAL INPUTS:
;   ra2  - right ascension coordinate array for second of paired
;          positions (same format requirements as RA)
;   dec2 - declination coordinate array for second of paired positions
;          (same format requirements as DEC)
;   sep  - separation in degrees (unless ARCMIN or ARCSEC keyword is
;          set); required if only one set of coordinates is passed
;   pa   - position angle (rotation angle left of North); required if
;          only onse set of coordinates is passed
;
; KEYWORDS:
;   arcmin - separations given/returned in arcminutes
;            (default=degrees) 
;   arcsec - separations given/returned in arcseconds
;            (default=degrees) 
;   silent - suppress output messages
;   other  - compute other set of coordinates for pair (RA and DEC
;            interpreted as coordinates for one pair member rather
;            than as central coordinates)
;
; OUTPUTS:
; 
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;
; PROCEDURES USED: ARM_DOUBLE(), ARM_COORDS()
;
; COMMENTS:  If one set of coordinates is passed with separations and
;            position angles, they are interpreted as center
;            coordinates and the corresponding paired coordinates are
;            returned (structure with RA1, DEC1, RA2 and DEC2
;            extensions).  If two sets of coordinates are passed, the
;            center coordinates, separations, and position angles are
;            returned (structure with RA, DEC, PA and SEP
;            extensions). 
;
;            NOTE: This has NOT yet been thoroughly beta-tested!
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2004 May 27
;    modified to call ARM_COORDS(), ARM, 2004 June 8
;-

function arm_midcoords, ra, dec, ra2, dec2, pa=pa, sep=sep, $
            arcmin=arcmin, arcsec=arcsec, silent=silent, other=other

; determine desired operation

    if N_ELEMENTS(ra2) eq 0L then begin
       if not KEYWORD_SET(silent) then MESSAGE, $
         'Central coordinates defined, computing coordinate pairs...', /continue
       getmid = 0L
    endif else begin
       if not KEYWORD_SET(silent) then MESSAGE, $
         'Two coordinate sets defined, computing central coordinates...', /continue
       getmid = 1L
    endelse

; error checking

    ON_ERROR, 1                 ; return to main level if error occurs

    if not getmid and (N_ELEMENTS(pa) eq 0L or N_ELEMENTS(sep) eq 0L) then $
      MESSAGE, 'Parameters PA and SEP must be defined.'

    n = N_ELEMENTS(ra)
    if n eq 0L then MESSAGE, 'No coordinates given.'

    if KEYWORD_SET(arcsec) and KEYWORD_SET(arcmin) then $
      MESSAGE, 'Keywords ARCMIN and ARCSEC cannot both be set.'

    if not getmid then begin
       if KEYWORD_SET(arcmin) then sep = sep / 6d1
       if KEYWORD_SET(arcsec) then sep = sep / 36d2 
    endif

; convert coordinates from string format to decimal degrees

    if STRMATCH(STRTRIM(ra[0], 2), '*:*') eq 1L then delimiter = ':'
    if STRMATCH(STRTRIM(ra[0], 2), '* *') eq 1L then delimiter = ' '
    if N_ELEMENTS(delimiter) ne 0L then begin
       rr1 = ARM_COORDS(ra,  delimiter=delimiter, /from_hours)
       dd1 = ARM_COORDS(dec, delimiter=delimiter)
       if getmid then begin
          rr2 = ARM_COORDS(ra2,  delimiter=delimiter, /from_hours)
          dd2 = ARM_COORDS(dec2, delimiter=delimiter)
       endif
    endif else begin
       rr1 = ARM_DOUBLE(ra)
       dd1 = ARM_DOUBLE(dec)
       if getmid then begin
          rr2  = ARM_DOUBLE(ra2)
          dd2 = ARM_DOUBLE(dec2)
       endif
    endelse

; perform desired operation

    if getmid then begin
       
       r1 = rr1
       r2 = rr2
       d1 = dd1
       d2 = dd2

; prevent "wrap-around error"

       wh = where(ABS(r1 - r2) gt 18d1, count)
       if count gt 0L then for i = 0L, count-1L do begin
          j = wh[i]
          if MIN([r1[j], r2[j]]) lt 360 - MAX([r1[j], r2[j]]) then $
            if r1[j] lt r2[j] then r1[j] = r1[j] + 360 else r2[j] = r2[j] + 360 $
          else if r1[j] gt r2[j] then r1[j] = r1[j] - 360 else r2[j] = r2[j] - 360
       endfor

       r = (r1 + r2) / 2.
       d = (d1 + d2) / 2.
       dr = (r1 - r2) * COS(d*!pi/18d1)
       dd = (d1 - d2)
       p = ATAN(dr / dd) * 18d1 / !pi
       s = SQRT(dr^2 + dd^2)

       if N_ELEMENTS(delimiter) ne 0L then begin
          r = ARM_COORDS(r, delimiter=delimiter, /to_hours)
          d = ARM_COORDS(d, delimiter=delimiter)
       endif

       if KEYWORD_SET(arcmin) then s = s * 6d1 else $
         if KEYWORD_SET(arcsec) then s = s * 36d2

      return, {ra: r, dec: d, pa: p, sep: s}

    endif else begin

       dd = sep / SQRT((TAN(pa*!pi/18d1))^2 + 1d0)
       dr = SQRT(sep^2 - dd^2) / COS(dd*!pi/18d1)

       if KEYWORD_SET(other) then number = 1d0 else number = 2d0

       r1 = rr1 + dr / number
       r2 = rr1 - dr / number
       d1 = dd1 - dd / number
       d2 = dd1 + dd / number

       if N_ELEMENTS(delimiter) ne 0L then begin
          r1 = ARM_COORDS(r1, delimiter=delimiter, /to_hours)
          d1 = ARM_COORDS(d1, delimiter=delimiter)
          if not KEYWORD_SET(other) then begin
             r2 = ARM_COORDS(r2, delimiter=delimiter, /to_hours)
             d2 = ARM_COORDS(d2, delimiter=delimiter)
          endif
       endif

       if KEYWORD_SET(other) then st = {ra1: ra, dec1: dec, ra2: r1, dec2: d1} $
         else st = {ra1: r1, dec1: d1, ra2: r2, dec2: d2}

       return, st

    endelse

 end
