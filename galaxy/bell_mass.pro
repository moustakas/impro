;+
; NAME:
;   BELL_MASS()
;
; PURPOSE:
;   Compute stellar mass-to-light ratios and stellar masses using
;   the Bell et al. color-based method. 
;
; INPUTS:
;   color  - photometric color (see METHOD)
;   absmag - absolute magnitude (see METHOD)
;
; OPTIONAL INPUTS:
;   method - technique to use, which also specifies COLOR and ABSMAG 
;     1: B-V color, V-band luminosity
;     2: B-R color, R-band luminosity
;     3: ....
;
;   color_err - uncertainty in COLOR
;   absmag_err - uncertainty in ABSMAG
;
; KEYWORD PARAMETERS:
;   log - return the logarithmic mass and error
;   salpeter - convert the output to a Salpeter (1955) IMF (default is
;     to use the diet Salpeter IMF)
;
; OUTPUTS:
;   mass - stellar mass (M_sun)
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Jan 27, NYU - based on a earlier, crappy piece
;     of code 
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

function mass_bell, ml, absmag, solar, absmag_err=absmag_err, $
  ml_err=ml_err, mass_err=mass_err
; compute the stellar mass
    
    if (n_elements(ml_err) eq 0L) then ml_err = ml*0.0
    if (n_elements(absmag_err) eq 0L) then absmag_err = absmag*0.0

    lum = 10.0^(-0.4*(absmag-solar))
    lum_err = alog(10.0)*0.4*absmag_err*lum

    mass = ml*lum
    mass_err = sqrt(ml_err^2.0+lum_err^2.0)
          
return, mass
end

function ml_bell, color, coeff, color_err=color_err, ml_err=ml_err
; compute the M/L ratio

    if (n_elements(color_err) eq 0L) then color_err = color*0.0

    ml = 10.0^(coeff[0] + coeff[1]*color)
    ml_err = ml*coeff[1]*color_err*alog(10.0)

return, ml
end

function bell_mass, color, absmag, method=method, color_err=color_err, $
  absmag_err=absmag_err, ml=ml, err_ml=err_ml, mass_err=mass_err, log=log, $
  salpeter=salpeter

    if (n_elements(color) eq 0L) or (n_elements(absmag) eq 0L) then begin
       doc_library, 'bell_mass'
       return, -1L
    endif

    if (n_elements(color_err) eq 0L) then color_err = color*0.0
    if (n_elements(lum_err) eq 0L) then lum_err = color*0.0

    bell01 = read_01bell() ; Bell & de Jong 2001

    if (n_elements(method) eq 0L) then method = 1

    case method of
       1: begin ; B-V, V
          coeff = [bell01[0].av,bell01[0].bv]
          filter = 'bessell_V.par'
       end
       2: begin ; B-R, R
          coeff = [bell01[1].ar,bell01[1].br]
          filter = 'bessell_R.par'
       end
       else: begin
          splog, 'Method not yet supported'
          return, color*0.0
       end
    endcase

; compute M/L and mass    
    solar = k_solar_magnitudes(filterlist=filter,/silent)
    ml = ml_bell(color,coeff,color_err=color_err,ml_err=err_ml)
    mass = mass_bell(ml,absmag,solar[0],absmag_err=absmag_err,$
      ml_err=err_ml,mass_err=mass_err)

    if keyword_set(salpeter) then begin
       mass = mass*10.0^(0.15)
       ml = ml*10.0^(0.15)
    endif
    
    if keyword_set(log) then begin
       mass_err = mass_err/alog(10.0)/(mass+(mass eq 0.0))*(mass ne 0.0)
       mass = alog10(mass+(mass eq 0.0))
    endif

return, mass
end    
