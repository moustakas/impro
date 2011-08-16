;+
; NAME:
;   COMPUTE_INCLINATION()
;
; PURPOSE:
;   Derive the inclination angle from the major-to-minor axis ratio of
;   a galaxy.  
;
; INPUTS: 
;   diam - input structure as returned by NED_WEBGET_DIAMETERS 
;
; KEYWORD PARAMETERS: 
;   twomass - use the 2MASS diameters (default is to use RC3) 
;
; OUTPUTS: 
;   incl - inclination angle, measured positive from North to East
;     (degrees) 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2006 Feb 16, U of A
;
; Copyright (C) 2006, John Moustakas
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

function compute_inclination, diam, twomass=twomass

    ndiam = n_elements(diam)
    if (ndiam eq 0L) then begin
       print, 'Syntax - incl = compute_inclination(diam,/twomass)'
       return, -1L
    endif

    incl = fltarr(ndiam)-999.0
    
    if keyword_set(twomass) then begin

       good = where((diam.twomass_k20_major_axis gt -900.0) and (diam.twomass_k20_minor_axis gt -900.0),ngood)
       if (ngood ne 0L) then begin

          d25_maj = float(diam[good].twomass_k20_major_axis)
          d25_min = float(diam[good].twomass_k20_minor_axis)

          ratio = d25_min/d25_maj ; "b" divided by "a"
          quantity = sqrt((1.0-ratio^2)/0.96)

          incl[good] = asin(quantity<1)*!radeg

       endif

    endif else begin

       good = where((diam.rc3_major_axis gt -900.0) and (diam.rc3_minor_axis gt -900.0),ngood)
       if (ngood ne 0L) then begin

          d25_maj = float(diam[good].rc3_major_axis)
          d25_min = float(diam[good].rc3_minor_axis)

          ratio = d25_min/d25_maj ; "b" divided by "a"
          quantity = sqrt((1.0-ratio^2)/0.96)

          incl[good] = asin(quantity<1)*!radeg
       
       endif

    endelse
       
return, incl
end
