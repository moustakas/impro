;+
; NAME:
;   IM_MEASURE_PEG()
;
; PURPOSE:
;   Given a Pegase data structure (read by IM_READ_PEG), measure some
;   broadband and spectroscopic properties of interest. 
;
; INPUTS:
;   peg - input structure from IM_READ_PEG
;
; OPTIONAL INPUTS:
;   age  - desired output ages [Myr]
;   mass - scale to this total mass
;
; KEYWORD PARAMETERS:
;   silent - suppress messages to STDOUT 
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   A disaster beyond your imagination will occur if MASS is
;   passed to IM_READ_PEG() *and* to this routine!
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Oct 18, U of A
;   jm06apr27uofa - added MASS optional input
;   jm08mar28nyu - consolidated the photometry into arrays; improved
;     the speed a bit, and cleaned up the documentation 
;   jm10mar13ucsd - some of the spectroscopic measurements (e.g.,
;     H-delta_A, and other EWs have been commented out)
;
; Copyright (C) 2005-2006, 2008, 2010, John Moustakas
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

function peg_qpint1d_func, x, wave=wave, flux=flux
return, interpol(flux,wave,x)
end

function im_measure_peg, peg, age=age, mass=mass, silent=silent

    if (n_elements(peg) eq 0) then begin
       doc_library, 'im_measure_peg'
       return, -1
    endif

; initialize the output data structure; optionally restrict the sample
; to a desired set of ages
    if (n_elements(age) eq 0) then begin
       nage = peg[0].nage
       ageindx = lindgen(nage)
    endif else begin
       get_element, peg.age, age, ageindx
       uindx = uniq(ageindx)
       ageindx = ageindx[uindx]
       age = age[uindx]
       nage = n_elements(age)
    endelse
    peg1 = peg[ageindx]

; scale by the stellar mass    
    if (n_elements(mass) eq 0) then mass = 1.0 ; [M_sun]
    
    peg1.flux        = peg1.flux*mass
    peg1.lineflux    = peg1.lineflux*mass
    peg1.mgalaxy     = peg1.mgalaxy*mass
    peg1.mstar       = peg1.mstar*mass
    peg1.mwd         = peg1.mwd*mass
    peg1.mnsbh       = peg1.mnsbh*mass
    peg1.msubstellar = peg1.msubstellar*mass
    peg1.mallstars   = peg1.mallstars*mass
    peg1.mgas        = peg1.mgas*mass
    peg1.lbol        = peg1.lbol*mass
    peg1.sfr         = peg1.sfr*mass
    peg1.nlyc        = peg1.nlyc*mass
    peg1.sniirate    = peg1.sniirate*mass
    peg1.sniarate    = peg1.sniarate*mass
    outpeg = peg1

; compute additional properties of interest; note that the magnitudes
; computed below are apparent magnitudes, absolute magnitudes, *and*,
; if the normalization is equal to one solar mass, the M/L ratio in
; each band; also note that the metallicities, etc., may not be
; defined for all star formation scenarios
    result = {$
      ubvrijhk: fltarr(8),$ ; Vega, z=0.0
      ugriz:    fltarr(5),$ ; AB, z=0.1
      galex:    fltarr(2),$ ; AB, z=0.0
      uvlum:    fltarr(2),$ ; continuum luminosity at 1500,2800 A [erg/s/A]
      uvfactor: fltarr(2),$ ; conversion factor from L(1500,2800 to total L(UV) (see Bell+05) 

      log12oh:       -99.0,$ ; oxygen metallicity, by number
      zoxygen:       -99.0,$ ; oxygen metallicity, by mass
      fgas:          -99.0,$ ; gas fraction
      yeff_tot:      -99.0,$ ; effective yield (all metals)
      yeff:          -99.0,$ ; effective yield (oxygen)

      rfrac:           0.0,$ ; return fraction

      lha:             0.0,$ ; [erg/s]
      lhb:             0.0}  ; [erg/s]
    result = struct_addtags(outpeg,replicate(result,nage))

; compute some gas and metallicity properties
    nonzero = where((result.zgas gt 0.0),nnzero)
    if (nnzero ne 0) then begin
       result[nonzero].log12oh = alog10(result[nonzero].zgas/0.0134) + 8.69 ; Z_sun,(O/H)_sun from Asplund+09
       result[nonzero].zoxygen = 12.0*10.0^(result[nonzero].log12oh-12.0) ; 

       result[nonzero].fgas = result[nonzero].mgas/result[nonzero].mgalaxy
;      result[nonzero].fgas = result.mgas/(result[nonzero].mgas+result[nonzero].mstar) ; which one??
             
       result[nonzero].yeff_tot = - result[nonzero].zgas/alog(result[nonzero].fgas)
       result[nonzero].yeff = -result[nonzero].zoxygen/alog(result[nonzero].fgas)
    endif

    result.rfrac = result.mgas/result.mgalaxy

; emission lines
    result.lha  = peg1.lineflux[1]*3.826D33 ; [erg/s]
    result.lhb  = peg1.lineflux[0]*3.826D33 ; [erg/s]
    
;add the SFR conversions the UV slope the total UV luminosity, etc.!    

; compute synthesized magnitudes and some other neat quantities
    ubvrijhk_filters = [bessell_filterlist(),twomass_filterlist()]
    ugriz_filters = sdss_filterlist()
    galex_filters = galex_filterlist()

    for iage = 0, nage-1 do begin
       if (keyword_set(silent) eq 0) then print, $
         format='("IM_MEASURE_PEG: Age ",I0,"/",I0,A0,$)', $
         iage+1, nage, string(13b)

       wave_edges = k_lambda_to_edges(peg1[iage].wave)
       newflux = peg1[iage].flux/(4*!dpi*(10.0*3.085678D18)^2) ; [erg/s/cm2/A] at 10 pc

       result[iage].ubvrijhk = -2.5*alog10(k_project_filters(wave_edges,$ ; Vega, z=0.0
         newflux,filterlist=ubvrijhk_filters,band_shift=0.0,/silent)) - $
         k_vega2ab(filterlist=ubvrijhk_filters,/kurucz,/silent)
       result[iage].ugriz = -2.5*alog10(k_project_filters(wave_edges,$    ; AB, z=0.1
         newflux,filterlist=ugriz_filters,band_shift=0.1,/silent))
       result[iage].galex = -2.5*alog10(k_project_filters(wave_edges,$    ; AB, z=0.0
         newflux,filterlist=galex_filters,band_shift=0.0,/silent))
       
; compute the UV luminosity at 1500 and 2800 A and the ratio of the
; luminosity at these wavelengths to the total UV luminosity
       luv = qpint1d('peg_qpint1d_func',1216D,3000D,functargs=$
         {wave:peg1[iage].wave,flux:peg1[iage].flux})
;      luv = qpint1d('peg_qpint1d_func',0.0D,3000.0D,functargs=$
;        {wave:peg1[iage].wave,flux:peg1[iage].flux})
       result[iage].uvlum = [1500.0,2800.0]*interpol(peg1[iage].flux,$
         peg1[iage].wave,[1500.0,2800.0])
       result[iage].uvfactor = luv/result[iage].uvlum
    endfor
    
return, result
end
