;+
; NAME:
;   ISEDFIT_LINECONTINUUM()
;
; PURPOSE:
;   Get the stellar continuum around the position of prominent nebular
;   emission lines.
;
; INPUTS: 
;   wave - *rest-frame* wavelength array [NPIX]
;   flux - corresponding *rest-frame* flux vector [NPIX]
;
; OPTIONAL INPUTS: 
;   linewave - *rest-frame* wavelengths of the emission lines of
;     interest (defaults to 3727.420, 4861.325, 5006.843, 6562.800,
;     corresponding to [OII], H-beta, [OIII], and H-alpha); the units
;     should be the same as WAVE (Angstrom by default) [NLINE]
;   vsigma - maximum emission-line velocity width (default 300 km/s) 
;
; KEYWORD PARAMETERS: 
;   debug - make a simple plot and wait for a keystroke 
;
; OUTPUTS: 
;   linecontinuum - stellar continuum at the position of each line in
;     LINEWAVE in the same units as FLUX [NLINE]
;
; OPTIONAL OUTPUTS:
;
; MODIFICATION HISTORY:
;    J. Moustakas, 2014 Jun 06, Siena
;
; Copyright (C) 2014, John Moustakas
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

function isedfit_linecontinuum, wave, flux, linewave=linewave, vsigma=vsigma, debug=debug

    if n_elements(linewave) eq 0 then linewave = $
      [3727.420,4861.325,5006.843,6562.800]
    linecontinuum = linewave*0

    if n_elements(vsigma) eq 0 then vsigma = 300.0 ; [km/s]

; build the emission-line mask    
    mask = fix(wave*0+1)
    for ll = 0, n_elements(linewave)-1 do begin
       factor = vsigma/im_light(/kms)*linewave[ll]
       mask = mask and (wave lt (linewave[ll]-4*factor) or (wave gt (linewave[ll]+4*factor)))
    endfor

; now loop through and get the continuum level    
    for ll = 0, n_elements(linewave)-1 do begin
       factor = vsigma/im_light(/kms)*linewave[ll]
       these = where(mask ne 0 and (((wave gt (linewave[ll]-10*factor)) and $
         (wave lt (linewave[ll]-4*factor))) or ((wave gt (linewave[ll]+4*factor)) and $
         (wave lt (linewave[ll]+10*factor)))),npix)
       linecontinuum[ll] = djs_median(medsmooth(flux[these],20))
       
       if keyword_set(debug) then begin
          ww = where(mask eq 0)
          djs_plot, wave, flux, xsty=3, ysty=3, xrange=linewave[ll]+20*[-1,1]*factor
          djs_oplot, wave[ww], flux[ww], psym=8, color='blue'
          djs_oplot, wave[these], flux[these], psym=8, color='orange'

          djs_oplot, !x.crange, linecontinuum[ll]*[1,1], color='red', line=0
          djs_oplot, linewave[ll]*[1,1], !y.crange, color='red', line=0
          cc = get_kbrd(1)
       endif
    endfor

return, linecontinuum
end
