;+
; NAME:
;   IM_READ_MARASTON05()
;
; PURPOSE:
;   Read a Maraston+05 SSP into a structure.
;
; INPUTS:
;   None required.
;
; OPTIONAL INPUTS:
;   metallicity - stellar metallicity (default 'z002')
;     z007  => [Z/H] = +0.67 (corresponding to Z=3.5 Zsun)
;     z004  => [Z/H] = +0.35 (corresponding to Z=2.0 Zsun)
;     z002  => [Z/H] = +0.00 (corresponding to Z=1.0 Zsun)
;     z001  => [Z/H] = -0.33 (corresponding to Z=0.5 Zsun)
;     z0001 => [Z/H] = -1.35 (corresponding to Z=1/50 Zsun)
;     z10m4 => [Z/H] = -2.25 (corresponding to Z=1/200 Zsun)
; 
;   hbmorph - horizontal branch morphology, either blue ('bhb') or red
;     ('rhb', default)
;
; KEYWORD PARAMETERS:
;   kroupa - read the Kroupa IMF files (default is to read the
;     Salpeter ones)
;   abmag - convert the output spectra to AB mag at 10 pc
;   flambda - convert the output spectra to F_lambda units (erg/s/cm^2/A) at 10 pc
;   fnu - convert the output spectra to F_nu units (erg/s/cm^2/Hz) at 10 pc 
;
; OUTPUTS:
;   mara - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   See http://www.icg.port.ac.uk/~maraston/Claudia's_Stellar_Population_Model.html
;   for additional relevant details.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 23, UCSD
;
; Copyright (C) 2011, John Moustakas
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

function im_read_maraston05, metallicity=metallicity, hbmorph=hbmorph, $
  kroupa=kroupa, abmag=abmag, flambda=flambda, fnu=fnu

    ssppath = getenv('IM_RESEARCH_DIR')+'/synthesis/maraston05/'

    if (n_elements(metallicity) eq 0) then metallicity = 'z002'
    case metallicity of
       'z007': begin
          zz = 0.07
          bracket_zh = +0.67
       end
       'z004': begin
          zz = 0.04
          bracket_zh = +0.35
       end
       'z002': begin
          zz = 0.02 ; solar
          bracket_zh = 0.0
       end
       'z001': begin
          zz = 0.01
          bracket_zh = -0.33
       end
       'z0001': begin
          zz = 0.001
          bracket_zh = -1.35
       end
       'z10m4': begin
          zz = 1E-4
          bracket_zh = -2.25
       end
       else: begin
          splog, 'Metallicity '+metallicity+' not recognized!'
          return, -1
       end
    endcase
    
    if (n_elements(hbmorph) eq 0) then hbmorph = 'rhb'
    if keyword_set(kroupa) then begin
       imfshort = 'kr'
       imf = 'Kroupa'
       massfile = ssppath+'stellarmass.kroupa'
    endif else begin
       imfshort = 'ss'
       imf = 'Salpeter'
       massfile = ssppath+'stellarmass.salpeter'
    endelse

; read the mass file
    if (file_test(massfile) eq 0) then begin
       splog, 'Mass file '+massfile+' not found!'
       return, -1
    endif
    readfast, massfile, allmass, skip=2
    
; read the SSP, parse, and pack into a data structure
    sspfile = ssppath+'sed.'+imfshort+metallicity+'.'+hbmorph
    if (file_test(sspfile) eq 0) then begin
       splog, 'SSP '+sspfile+' not found!'
       return, -1
    endif

    splog, 'Reading '+sspfile
    readfast, sspfile, data, skip=0, /double
    allage = reform(data[0,*])
    age = allage[uniq(allage,sort(allage))]
    nage = n_elements(age)

    for ii = 0L, nage-1 do begin
       these = where(age[ii] eq allage,npix)
       if (ii eq 0) then begin
          mara = {$
            imf:                imf,$
            Z:                   zz,$
            age:            age*1D9,$ ; [yr]
            mstar:     fltarr(nage),$ ; [Msun]
            wave: float(reform(data[2,these])),$
            flux:      fltarr(npix,nage)} ; [erg/s/AA/Msun]
       endif
       mara.flux[*,ii] = reform(data[3,these])
; get the total stellar mass (living stars plus remnants)  
       mthese = where(allmass[0,*] eq bracket_zh)
       mara.mstar = interpol(allmass[2,mthese],allmass[1,mthese]*1D9,mara.age)
    endfor

; convert the units of the spectra as desired
    pc10 = 10.0*3.085678D18     ; =10 pc
    light = 2.99792458D18       ; [Angstrom/s]
    if keyword_set(flambda) then mara.flux = mara.flux/(4.0*!dpi*pc10^2)
    if keyword_set(fnu) then for ii = 0, nage-1 do $
      mara.flux[*,ii] = mara.flux[*,ii]/(4.0*!dpi*pc10^2)*mara.wave^2/light
       
    if keyword_set(abmag) then begin
       for ii = 0, nage-1 do begin
          mara.flux[*,ii] = mara.flux[*,ii]/(4.0*!dpi*pc10^2)*mara.wave^2/light
          gd = where(mara.flux[*,ii] gt 0.0,ngd)
          if (ngd ne 0) then mara.flux[gd,ii] = -2.5*alog10(mara.flux[gd,ii])-48.6
       endfor
    endif
    
return, mara
end
