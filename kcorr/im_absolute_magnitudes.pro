;+
; NAME:
;       IM_ABSOLUTE_MAGNITUDES()
;
; PURPOSE:
;       Given an apparent magnitude and luminosity distance, compute
;       the absolute magnitude and luminosity.
;
; INPUTS:
;       band     - scalar string specifying the bandpass (e.g., 'U')  
;       mag      - apparent magnitude in BAND [NOBJECT, mag] 
;       distance - luminosity distance [NOBJECT, Mpc]
;
; OPTIONAL INPUTS:
;       mag_err      - error in MAG [NOBJECT, mag]
;       distance_err - error in DISTANCE [NOBJECT, Mpc]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       result - 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine could be generalized to take a vector BAND.  
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jun 02, U of A
;
; Copyright (C) 2005, John Moustakas
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

function im_absolute_magnitudes, band, mag, distance, mag_err=mag_err, $
  distance_err=distance_err

    if (n_elements(band) eq 0L) or (n_elements(mag) eq 0L) or (n_elements(distance) eq 0L) then begin
       print, 'Syntax - result = im_absolute_magnitudes(band,mag,distance,$'
       print, '   mag_err=, distance_err=)'
       return, -1L
    endif

    nobject = n_elements(mag)
    if (n_elements(distance) ne nobject) then begin
       print, 'MAG and DISTANCE must have the same number of elements.'
       return, -1L
    endif

    if (n_elements(mag_err) eq 0L) then mag_err = mag*0.0 else begin
       if (n_elements(mag_err) ne nobject) then begin
          print, 'MAG and MAG_ERR must have the same number of elements.'
          return, -1L
       endif
    endelse
    
    if (n_elements(distance_err) eq 0L) then distance_err = distance*0.0 else begin
       if (n_elements(distance_err) ne nobject) then begin
          print, 'DISTANCE and DISTANCE_ERR must have the same number of elements.'
          return, -1L
       endif
    endelse

; initialize the output data structure

    result = {$
      mag:          0.0, $
      mag_err:      0.0, $
      distance:     0.0, $
      distance_err: 0.0, $
      absmag:       0.0, $
      absmag_err:   0.0, $
      lum:          0.0, $ ; log [L_sun]
      lum_err:      0.0}
    result = replicate(result,nobject)

; select the right band

    case strupcase(strcompress(band,/remove)) of
       'U'          : filter = 'bessell_U.par'
       'B'          : filter = 'bessell_B.par'
       'V'          : filter = 'bessell_V.par'
       'R'          : filter = 'bessell_R.par'
       'I'          : filter = 'bessell_I.par'
       'TWOMASS_J'  : filter = 'twomass_J.par'
       'TWOMASS_H'  : filter = 'twomass_H.par'
       'TWOMASS_KS' : filter = 'twomass_Ks.par'
       'SDSS_U'     : filter = 'sdss_u0.par'
       'SDSS_G'     : filter = 'sdss_g0.par'
       'SDSS_R'     : filter = 'sdss_r0.par'
       'SDSS_I'     : filter = 'sdss_i0.par'
       'SDSS_Z'     : filter = 'sdss_z0.par'
       else: begin
          print, 'Unrecognized bandpass.'
          return, -1L
       end
    endcase

    filtinfo = im_filterspecs(filterlist=filter)
    absmsun = filtinfo.solarmags
    
; do the calculation    
    
    result.mag     = mag
    result.mag_err = mag_err
    
    result.distance     = distance
    result.distance_err = distance_err
    
    result.absmag     = mag - 5.0*alog10(distance) - 25.0
    result.absmag_err = sqrt( mag_err^2.0 + ((5.0/alog(10))*distance_err/distance)^2.0 )

    lum = 10.0D^(-0.4*(result.absmag-absmsun))
    lum_err = alog(10.0) * 0.4D * lum * result.absmag_err 
    
    result.lum_err = lum_err/lum/alog(10.0)
    result.lum = alog10(lum)

return, result
end
    
