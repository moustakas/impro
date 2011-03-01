;+
; NAME:
;       CONTINUUM_PROPERTIES()
;
; PURPOSE:
;       Compute properties of the stellar continuum based on
;       template-fitting with the Bruzual & Charlot (2003) models.  
;
; CALLING SEQUENCE:
;       properties = continuum_properties(linefit,ancillary=,disttag=)
;
; INPUTS:
;       linefit - ISPECLINEFIT() data linefiture with a DISTTAG field 
;
; OPTIONAL INPUTS:
;       ancillary - ancillary data (specifically, the distance)
;       disttag   - name of the distance tag name (default DISTANCE)  
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       properties - 
;
; OPTIONAL OUTPUTS:
;
;
; PROCEDURES USED:
;       SPLOG, TSUM()
;
; COMMENTS:
;
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Mar 09, U of A, written, based on an
;          earlier code 
;       jm04jul29uofa - additional error checking
;       jm05jul25uofa - output changed to floating-point precision 
;       jm05nov02uofa - added ANCILLARY optional input
;
; Copyright (C) 2004-2005, John Moustakas
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

function continuum_properties, linefit, ancillary=ancillary, disttag=disttag
    
    nspec = n_elements(linefit)

    if (nspec eq 0L) then begin
       print, 'Syntax - properties = continuum_properties(linefit,ancillary=,disttag=)'
       return, -1L
    endif

    if (n_elements(disttag) eq 0L) then disttag = 'DISTANCE'
    
    if (tag_exist(ancillary,disttag) eq 0) then begin
       splog, 'No DISTTAG field in ANCILLARY!'
       return, -1L
    endif else begin
       match = where(disttag eq tag_names(ancillary),nmatch)
       if (nmatch eq 1L) then distance = ancillary.(match) else begin
          splog, 'Problem finding DISTTAG field in ANCILLARY.'
          return, -1L
       endelse
    endelse 
    
    Mpc2cm = 3.086D24 ; [cm/Mpc]

    if (tag_exist(linefit,'NTEMPLATE') eq 0) then begin
       splog, 'No NTEMPLATE field in LINEFIT.'
       return, -1L
    endif 
    
    ntemplate = linefit[0].ntemplate
    init = fltarr(ntemplate) ; dblarr(ntemplate)
    properties = {$
      continuum_age:            0.0, $ ; mass-weighted age [Myr]
      continuum_age_U:          0.0, $ ; luminosity-weighted age [Myr]
      continuum_age_B:          0.0, $
      continuum_age_V:          0.0, $
      continuum_age_R:          0.0, $
;     continuum_age_I:          0.0, $
;     continuum_age_J:          0.0, $
;     continuum_age_H:          0.0, $
;     continuum_age_Ks:         0.0, $
      continuum_fraction_U:    init, $ ; light fraction
      continuum_fraction_B:    init, $
      continuum_fraction_V:    init, $
      continuum_fraction_R:    init, $
;     continuum_fraction_I:    init, $
;     continuum_fraction_J:    init, $
;     continuum_fraction_H:    init, $
;     continuum_fraction_Ks:   init, $
      continuum_template_mass: init, $ ; total mass in each template [M_sun]
      continuum_total_mass:     0.0, $ ; total galaxy mass [M_sun]
      continuum_mass_fraction: init, $ ; mass fraction
      continuum_sfr:            0.0, $ ; star formation rate [M_sun/yr]
;     continuum_sfh:           init, $ ; star formation history [M_sun/yr]
      continuum_birthrate:      0.0  $ ; Scalo birthrate parameter
      }
    properties = replicate(properties,nspec)

    age = linefit[0].template_age ; [yr]
    maxage = max(age,maxageindx)
    nage = n_elements(age)
;   logage = alog10(age)
;
;   dlogage = fltarr(nage-1L)
;   for i = 0L, nage-2L do dlogage[i] = (logage[i+1]-logage[i])/2.0
;
;   deltat = fltarr(nage)
;   for i = 0L, nage-1L do begin
;
;      case i of
;         0L: deltat[0L] = 10^dlogage[0L]
;         nage-1L: deltat[nage-1L] = 10^dlogage[nage-2L]
;         else: deltat[i] = 10^dlogage[i] + 10^dlogage[i-1L]
;      endcase
;         
;   endfor
;
;   stop
    
    for k = 0L, nspec-1L do begin

       if (distance[k] gt -900.0) then begin

          dist = distance[k]*Mpc2cm ; [cm]
          area = 4.0*!dpi*dist*dist ; [cm2]

          properties[k].continuum_template_mass = linefit[k].continuum_coeff*area           ; [M_sun]
          properties[k].continuum_total_mass = total(properties[k].continuum_template_mass) ; [M_sun]

;         properties[k].continuum_sfh = properties[k].continuum_template_mass/age       ; [M_sun/yr]

          properties[k].continuum_sfr = properties[k].continuum_template_mass[0]/1D7    ; [M_sun/yr] <-- NOTE!!
;         properties[k].continuum_sfr = properties[k].continuum_template_mass[0]/age[0] ; [M_sun/yr]

          properties[k].continuum_birthrate = properties[k].continuum_sfr * age[maxageindx] / $
            properties[k].continuum_total_mass

          mfraction = properties[k].continuum_template_mass/properties[k].continuum_total_mass
          properties[k].continuum_mass_fraction  = mfraction
          properties[k].continuum_age  = tsum(age,age*mfraction)/tsum(age,mfraction)/1D6
          
; luminosity-weighted ages

          lum = mfraction/linefit[k].template_ml_u
          properties[k].continuum_age_u = tsum(age,age*lum)/tsum(age,lum)/1D6
          
          lum = mfraction/linefit[k].template_ml_b
          properties[k].continuum_age_b = tsum(age,age*lum)/tsum(age,lum)/1D6

          lum = mfraction/linefit[k].template_ml_v
          properties[k].continuum_age_v = tsum(age,age*lum)/tsum(age,lum)/1D6

          lum = mfraction/linefit[k].template_ml_r
          properties[k].continuum_age_r = tsum(age,age*lum)/tsum(age,lum)/1D6
          
;         lum = mfraction/linefit[k].template_ml_i
;         properties[k].continuum_age_i = tsum(age,age*lum)/tsum(age,lum)/1D6
;
;         lum = mfraction/linefit[k].template_ml_j
;         properties[k].continuum_age_j = tsum(age,age*lum)/tsum(age,lum)/1D6
;
;         lum = mfraction/linefit[k].template_ml_h
;         properties[k].continuum_age_h = tsum(age,age*lum)/tsum(age,lum)/1D6
;
;         lum = mfraction/linefit[k].template_ml_ks
;         properties[k].continuum_age_ks = tsum(age,age*lum)/tsum(age,lum)/1D6

; light fractions in various bands

          lfraction = properties[k].continuum_template_mass/linefit[k].template_ml_u
          properties[k].continuum_fraction_u = lfraction/total(lfraction)

          lfraction = properties[k].continuum_template_mass/linefit[k].template_ml_b
          properties[k].continuum_fraction_b = lfraction/total(lfraction)

          lfraction = properties[k].continuum_template_mass/linefit[k].template_ml_v
          properties[k].continuum_fraction_v = lfraction/total(lfraction)

          lfraction = properties[k].continuum_template_mass/linefit[k].template_ml_r
          properties[k].continuum_fraction_r = lfraction/total(lfraction)

;         lfraction = properties[k].continuum_template_mass/linefit[k].template_ml_i
;         properties[k].continuum_fraction_i = lfraction/total(lfraction)
;
;         lfraction = properties[k].continuum_template_mass/linefit[k].template_ml_j
;         properties[k].continuum_fraction_j = lfraction/total(lfraction)
;
;         lfraction = properties[k].continuum_template_mass/linefit[k].template_ml_h
;         properties[k].continuum_fraction_h = lfraction/total(lfraction)
;
;         lfraction = properties[k].continuum_template_mass/linefit[k].template_ml_ks
;         properties[k].continuum_fraction_ks = lfraction/total(lfraction)
          
       endif

    endfor

return, properties
end    
