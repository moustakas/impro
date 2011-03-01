;+
; NAME:
;       SPECLINEFIT_LOCATE()
;
; PURPOSE:
;       Locate a particular object (galaxy) in the parsed data
;       structure as written by PARSE_ISPECLINEFIT() and read by
;       IRDSPECLINEFIT(). 
;
; CALLING SEQUENCE:
;       indx - speclinefit_locate(speclinefit,object,objtagname=,/exact)
;
; INPUTS:
;       speclinefit - data structure read by IRDSPECLINEFIT()
;       object      - locate this object in the tag OBJTAGNAME
;                     (string array or scalar) 
;
; OPTIONAL INPUTS:
;       objtagname - SPECLINEFIT structure tag name in which to search
;                    for OBJECT (default 'GALAXY')
;
; KEYWORD PARAMETERS:
;       exact - match exactly
;
; OUTPUTS:
;       indx - index number in SPECLINEFIT corresponding to OBJECT 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;       IDL> linefit = irdspeclinefit()
;       IDL> print, speclinefit_locate(linefit,'NGC4736')
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Mar 14, U of A - written, based on an
;         earlier code
;       jm06oct19nyu - added EXACT optional input
;       jm08feb10nyu - bug fix: failed to pass EXACT when OBJECT was
;                      an array
; 
; Copyright (C) 2004, 2006, 2008, John Moustakas
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

function speclinefit_locate, speclinefit, object, objtagname=objtagname, exact=exact

    nspeclinefit = n_elements(speclinefit)
    nobject = n_elements(object)

    if (n_elements(objtagname) eq 0L) then objtagname = 'GALAXY'

    if (nspeclinefit eq 0L) then begin
       print, 'Syntax - indx = speclinefit_locate(speclinefit,object,objtagname=)'
       return, -1L
    endif

    if (nobject eq 0L) then return, lindgen(nspeclinefit)

    tagindx = where(objtagname eq tag_names(speclinefit),ntagindx)
    if (ntagindx ne 1L) then begin
       message, 'Multiple or no tags found matching '+OBJTAGNAME+'.', /info
       return, -1L
    endif
    
    if (nobject gt 1L) then begin
       indx = lonarr(nobject)
       for k = 0L, nobject-1L do indx[k] = speclinefit_locate(speclinefit,object[k],exact=exact)
    endif else begin
       if keyword_set(exact) then indx = where(strmatch(strtrim(speclinefit.(tagindx),2),strtrim(object,2),/fold) eq 1B) else $
         indx = where(strmatch(strtrim(speclinefit.(tagindx),2),'*'+strtrim(object,2)+'*',/fold) eq 1B)
    endelse

return, indx
end
