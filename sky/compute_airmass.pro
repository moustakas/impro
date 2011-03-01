;+
; NAME:
;	COMPUTE_AIRMASS()
;
; PURPOSE:
;	Compute the airmass of an observation.
;
; INPUTS:
;	latitude - observatory latitude [decimal degrees]
;	dec      - declination [decimal degrees]
;	ha       - hour angle [decimal degrees]
;
; OUTPUTS:
;	airmass  - airmass of the observation
;
; OPTIONAL OUTPUTS:
;	zd       - zenith distance of the observation [degrees]
;
; COMMENTS:
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 Feb 13, U of A
;
; Copyright (C) 2002, 2007, John Moustakas
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

function compute_airmass, latitude, dec, ha, zd=zd

    if (n_params() ne 3L) then begin
       doc_library, 'compute_airmass'
       return, 0.0
    endif
    
; atmospheric scale height [Allen 1973, p. 125, 133]

    scale = 750.0D 

; compute the zenith distance

    cos_zd = sin(latitude*!dtor)*sin(dec*!dtor) + $
      cos(latitude*!dtor)*cos(dec*!dtor)*cos(ha*!dtor)
    zd = acos(cos_zd)*!radeg
       
    x = scale * cos_zd

    airmass = sqrt(x^2.0D + 2.0*scale + 1.0) - x

return, airmass
end
