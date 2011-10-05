;+
; NAME:
;   IM_READ_PEGASE()
;
; PURPOSE:
;   Read the output of CREATE_PEGASE_SSPs into a structure. 
;
; INPUTS:
;   None required.
;
; OPTIONAL INPUTS:
;   metallicity - stellar metallicity to read; to see the full list 
;     just call this routine with
;
;     IDL> ssp = im_read_pegase(metallicity='Z')
;
; KEYWORD PARAMETERS:
;   kroupa - read the Kroupa+01 IMF files (default is to read the 
;     Salpeter ones)
;   abmag - convert the output spectra to AB mag at 10 pc
;   flambda - convert the output spectra to F_lambda units (erg/s/cm^2/A) at 10 pc
;   fnu - convert the output spectra to F_nu units (erg/s/cm^2/Hz) at 10 pc 
;
; OUTPUTS:
;   pegase - output data structure
;     Z - metallicity
;     age - age vector [NAGE]
;     mstar - mass in stars [NAGE]
;     wave - wavelength vector [NPIX] (Angstrom)
;     flux - flux vector [NPIX,NAGE] (erg/s/A)
;     ... - and many other goodies
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Oct 05, UCSD
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

function im_read_pegase, metallicity=metallicity, kroupa=kroupa, $
  abmag=abmag, flambda=flambda, fnu=fnu

    ssppath = getenv('PEGASE_HR_DIR')+'/SSPs/'

; defaults
    imf = 'salp'
    if keyword_set(kroupa) then imf = 'kroupa01'

; metallicity
    if (n_elements(metallicity) eq 0) then metallicity = 'Z0.02'
    allz = 'Z'+['0.0001','0.0004','0.004','0.008','0.02','0.05','0.1']
    match, allZ, metallicity, m1, m2
    if (m1[0] eq -1) then begin
       splog, 'Supported values of METALLICITY:'
       for ii = 0, n_elements(allZ)-1 do print, '  '+allZ[ii]
       return, -1
    endif
    zz = float(strmid(metallicity,1))

    sspfile = ssppath+'SSP_'+imf+'_'+metallicity+'.fits'
    if (file_test(sspfile) eq 0) then begin
       splog, 'SSP '+sspfile+' not found!'
       return, -1
    endif

; read the SSP using IM_READ_PEG(), but convert to a better format 
    splog, 'Reading '+sspfile
    peg = im_read_peg(sspfile)
    info = struct_trimtags(peg,except=['file','nage','ncont','nlines',$
      'wave','flux','linewave','lineflux'])
    info.age *= 1D6 ; [yr]
    nage = n_elements(info)

    pegase = {Z: zz}
    tags = tag_names(info)
    for ii = 0, n_elements(tags)-1 do pegase = create_struct(pegase,$
      tags[ii],make_array(nage,type=size(info[0].(ii),/type)))
    for ii = 0, n_elements(tags)-1 do pegase.(ii+1) = info.(ii)
    pegase = struct_addtags(pegase,{wave: peg[0].wave, flux: peg.flux, $
      lineflux: peg.lineflux})

; convert the units of the spectra as desired
    pc10 = 10.0*3.085678D18     ; =10 pc
    light = 2.99792458D18       ; [Angstrom/s]
    if keyword_set(flambda) then pegase.flux = pegase.flux/(4.0*!dpi*pc10^2)
    if keyword_set(fnu) then for ii = 0, nage-1 do $
      pegase.flux[*,ii] = pegase.flux[*,ii]/(4.0*!dpi*pc10^2)*pegase.wave^2/light
       
    if keyword_set(abmag) then begin
       for ii = 0, nage-1 do begin
          pegase.flux[*,ii] = pegase.flux[*,ii]/(4.0*!dpi*pc10^2)*pegase.wave^2/light
          gd = where(pegase.flux[*,ii] gt 0.0,ngd)
          if (ngd ne 0) then pegase.flux[gd,ii] = -2.5*alog10(pegase.flux[gd,ii])-48.6
       endfor
    endif
    
return, pegase
end
