;+
; NAME:
;   IM_SPHEREMATCH()
;
; PURPOSE:
;   Simple iterative wrapper on spherematch. 
;
; INPUTS: 
;   cat1 - first catalog
;   cat2 - second catalog
;
; OPTIONAL INPUTS: 
;   ratagname1 - structure tag name containing the right ascension (in
;     decimal degrees) in the first catalog (default 'RA')
;   dectagname1 - structure tag name containing the declination (in
;     decimal degrees) in the first catalog (default 'DEC') 
;   ratagname2 - structure tag name containing the right ascension (in
;     decimal degrees) in the second catalog (default 'RA')
;   dectagname2 - structure tag name containing the declination (in
;     decimal degrees) in the second catalog (default 'DEC')
;   radius - spherematching radius in arcseconds (default 1 arcsec) 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   match1 - index array of matching objects in CAT1
;   match2 - index array of matching objects in CAT2
;
; OPTIONAL OUTPUTS:
;   raoff - median offset in right ascension (decimal degrees) 
;   decoff - median offset in declination (decimal degrees) 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   2010-Jan-16, UCSD - J. Moustakas, UCSD
;
; Copyright (C) 2009, John Moustakas
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

function im_spherematch, cat1, cat2, match2=match2, $
  ratagname1=ratagname1, dectagname1=dectagname1, $
  ratagname2=ratagname2, dectagname2=dectagname2, $
  radius=radius, raoff=raoff, decoff=decoff

    if (n_elements(cat1) eq 0) or (n_elements(cat2) eq 0) then begin
       doc_library, 'im_spherematch'
       return, -1
    endif
    
    if (n_elements(ratagname1) eq 0) then ratagname1 = 'RA'
    if (n_elements(dectagname1) eq 0) then dectagname1 = 'DEC'
    if (n_elements(ratagname2) eq 0) then ratagname2 = 'RA'
    if (n_elements(dectagname2) eq 0) then dectagname2 = 'DEC'
    if (n_elements(radius) eq 0) then radius = 1.0 ; [arcsec]

    ratag1 = tag_indx(cat1[0],ratagname1)
    dectag1 = tag_indx(cat1[0],dectagname1)
    if (ratag1[0] eq -1) or (dectag1[0] eq -1) then $
      message, 'RA and/or DEC tags in CAT1 found'

    ratag2 = tag_indx(cat2[0],ratagname2)
    dectag2 = tag_indx(cat2[0],dectagname2)
    if (ratag2[0] eq -2) or (dectag2[0] eq -2) then $
      message, 'RA and/or DEC tags in CAT2 found'

; match twice, computing and applying the median astrometric offset
; after the first iteration
    spherematch, cat1.(ratag1), cat1.(dectag1), cat2.(ratag2), $
      cat2.(dectag2), radius/3600.0, match1, match2
;   help, match1, match2
    if (match1[0] eq -1) then begin
       raoff = 0.0
       decoff = 0.0
    endif else begin
       raoff = djs_median(cat2[match2].(ratag2)-cat1[match1].(ratag1))
       decoff = djs_median(cat2[match2].(dectag2)-cat1[match1].(dectag1))
       spherematch, cat1.(ratag1)+raoff, cat1.(dectag1)+decoff, $
         cat2.(ratag2), cat2.(dectag2), radius/3600.0, match1, match2
;      help, match2, match1
    endelse
    
return, match1
end
