;+
; NAME:
;   MAKE_WAVE()
;
; PURPOSE:
;   Calculate a wavelength vector of a one-dimensional spectrum
;   using the header parameters.
;
; INPUTS:
;   header - FITS header
;
; OUTPUTS:
;   wave - output wavelength vector
;
; OPTIONAL OUTPUTS:
;   cd1_1  - dispersion [Angstrom/pixel]
;   crval1 - starting wavelength [Angstrom]
;
; PROCEDURES USED:
;   SXPAR()
;
; MODIFICATION HISTORY:
;   J. Moustakas, 1999 August 10, UCB, written
;   jm01jul11uofa - documented and made a function
;   jm01oct19uofa - simplified
;   jm02mar12uofa - added CD1_1 and CRVAL1 as optional outputs 
;   jm03apr12uofa - added error checking
;   jm05apr18uofa - check if CDELT1 is present instead of CD1_1 
;   jm10apr28ucsd - double precision, by default
;
; Copyright (C) 1999, John Moustakas
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

function make_wave, header, cd1_1=cd1_1, crval1=crval1

    nhead = n_elements(header)
    if nhead eq 0L then begin
       doc_library, 'make_wave'
       return, -1L
    endif
    
    npix = sxpar(header,'NAXIS1',count=count1)
    crval1 = sxpar(header,'CRVAL1',count=count2)
    refpix = sxpar(header,'CRPIX1',count=count3)-1L ; IDL is zero-indexed
    cd1_1 = sxpar(header,'CD1_1',count=count4)	
    cdelt1 = sxpar(header,'CDELT1',count=count5)

    if (count4 eq 0L) and (count5 eq 1L) then begin
       cd1_1 = cdelt1
       count4 = count5
    endif

    if (float(count1+count2+count3+count4))[0] ne 4.0 then begin
       print, 'Insufficient wavelength information'
       return, -1L
    endif
    
    indx = dindgen(npix)-refpix
    wave = crval1 + indx*cd1_1

return, wave
end
