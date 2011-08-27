;+
; NAME:
;   MAKE_LNWAVE()
;
; PURPOSE:
;   Calculate a natural-log wavelength vector of a one-dimensional
;   spectrum using the header parameters.
;
; INPUTS:
;   header - FITS header
;
; OUTPUTS:
;   wave - output alog(wavelength) vector
;
; OPTIONAL OUTPUTS:
;   cd1_1  - dispersion [ln-Angstrom/pixel]
;   crval1 - starting wavelength [ln-Angstrom]
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009-Nov-10, UCSD - based on MAKE_WAVE()
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

function make_lnwave, header, cd1_1=cd1_1, crval1=crval1

    nhead = n_elements(header)
    if (nhead eq 0) then begin
       doc_library, 'make_lnwave'
       return, -1
    endif
    
    npix = sxpar(header,'NAXIS1',count=count1)
    crval1 = sxpar(header,'CRVAL1',count=count2)
    refpix = sxpar(header,'CRPIX1',count=count3)-1 ; IDL is zero-indexed
    cd1_1 = sxpar(header,'CD1_1',count=count4)	
    cdelt1 = sxpar(header,'CDELT1',count=count5)

    if (count4 eq 0) and (count5 eq 1) then begin
       cd1_1 = cdelt1
       count4 = count5
    endif

    if (float(count1+count2+count3+count4))[0] ne 4.0 then begin
       print, 'Insufficient wavelength information'
       return, -1
    endif

    wave = crval1 + (dindgen(npix)-refpix)*cd1_1

return, wave
end
