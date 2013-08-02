;+
; NAME:
;   IM_READ_FSPS()
;
; PURPOSE:
;   Read an FSPS SSP into a structure.
;
; INPUTS:
;   None required.
;
; OPTIONAL INPUTS:
;   metallicity - stellar metallicity to read; the choices depends on
;     the stellar library and isochrones desired; to see the full list
;     just call this routine with
;
;     IDL> ssp = im_read_fsps(metallicity='Z')
;
; KEYWORD PARAMETERS:
;   basti - read the BaSTI isochrones (default is to use Padova) (not
;     supported as of v2.3 of FSPS!)
;   miles - read the MILES stellar library (available just in the
;     optical) (default is to use just the low-resolution BaSeL
;     stellar library, which has much broader wavelength coverage) 
;   kroupa - read the Kroupa+01 IMF files (default is to read the 
;     Salpeter ones)
;   chabrier - read the Chabrier+03 IMF files
;   abmag - convert the output spectra to AB mag at 10 pc
;   flambda - convert the output spectra to F_lambda units (erg/s/cm^2/A) at 10 pc
;   fnu - convert the output spectra to F_nu units (erg/s/cm^2/Hz) at 10 pc 
;
; OUTPUTS:
;   fsps - output data structure
;     Z - metallicity
;     age - age vector [NAGE]
;     mstar - mass in stars [NAGE]
;     lbol - bolometric luminosity [NAGE]
;     wave - wavelength vector [NPIX] (Angstrom)
;     flux - flux vector [NPIX,NAGE] (erg/s/A)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   See https://www.cfa.harvard.edu/~cconroy/FSPS.html for additional
;   relevant details. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 30, UCSD
;   jm11mar14ucsd - updated to the latest SSPs
;   jm11mar28ucsd - read the BaSeL library by default 
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

function im_read_fsps, metallicity=metallicity, basti=basti, $
  miles=miles, kroupa=kroupa, chabrier=chabrier, abmag=abmag, $
  flambda=flambda, fnu=fnu

    fspspath = getenv('IM_RESEARCH_DIR')+'/synthesis/fsps/'
    ssppath = fspspath+'SSP/'

; defaults
    if keyword_set(miles) then lib = 'MILES' else lib = 'BaSeL' ; stellar library
    if keyword_set(basti) then begin
       iso = 'BaSTI' 
       splog, 'The BaSTI isochrones are not supported as of FSPS v2.3!!'
       return, -1
    endif
    iso = 'Padova' ; isochrones

    imf = 'Salpeter'
    if keyword_set(kroupa) then imf = 'Kroupa'
    if keyword_set(chabrier) then imf = 'Chabrier'

; metallicity
    case strlowcase(iso) of
       'padova': begin
          if (n_elements(metallicity) eq 0) then metallicity = 'Z0.0190'
          if keyword_set(miles) then begin ; MILES
             allZ = ['Z0.0008','Z0.0031','Z0.0096','Z0.0190','Z0.0300']
          endif else begin ; BaSeL
             allZ = ['Z0.0002','Z0.0003','Z0.0004','Z0.0005','Z0.0006',$
               'Z0.0008','Z0.0010','Z0.0012','Z0.0016','Z0.0020','Z0.0025',$
               'Z0.0031','Z0.0039','Z0.0049','Z0.0061','Z0.0077','Z0.0096',$
               'Z0.0120','Z0.0150','Z0.0190','Z0.0240','Z0.0300']
          endelse
       end
       'basti': begin
          if (n_elements(metallicity) eq 0) then metallicity = 'Z0.0200'
          if keyword_set(hires) then begin ; MILES
             allZ = ['Z0.0006','Z0.0040','Z0.0100','Z0.0200','Z0.0300']
          endif else begin ; BaSeL
             allZ = ['Z0.0003','Z0.0006','Z0.0010','Z0.0020','Z0.0040',$
               'Z0.0080','Z0.0100','Z0.0200','Z0.0300','Z0.0400']
          endelse 
       end
    endcase
    match, allZ, metallicity, m1, m2
    if (m1[0] eq -1) then begin
       splog, 'Supported values of METALLICITY for the '+$
         lib+' stellar library and the '+iso+' isochrones:'
       for ii = 0, n_elements(allZ)-1 do print, '  '+allZ[ii]
       return, -1
    endif
    zz = float(strmid(metallicity,1))

    sspfile = ssppath+'SSP_'+iso+'_'+lib+'_'+imf+'_'+metallicity+'.out.spec'
    if (file_test(sspfile) eq 0) then begin
       splog, 'SSP '+sspfile+' not found!'
       return, -1
    endif

;; read the wavelength array - only for v2.2 and earlier!
;    if strlowcase(lib) eq 'miles' then $
;      wavepath = fspspath+'SPECTRA/MILES/' else $
;        wavepath = fspspath+'SPECTRA/BaSeL3.1/'
;    wavefile = wavepath+strlowcase(lib)+'.lambda'
;    if (file_test(wavefile) eq 0) then begin
;       splog, 'Wavelength file '+wavefile+' not found!'
;       return, -1
;    endif
;    splog, 'Reading '+wavefile
;    readcol, wavefile, wave, format='F', /silent
;    npix = n_elements(wave)

; read the SSP -- see Conroy's code READ_SPEC
    splog, 'Reading '+sspfile
    openr, lun, sspfile, /get_lun
    char = '#' ; burn the header
    while (strmid(char,0,1) eq '#') do readf, lun, char
    char = strsplit(char,' ',/extr)

    nage = long(char[0])
    npix = long(char[1])

    wave = fltarr(npix)
    readf, lun, wave

; convert to air!
    airwave = wave
    vactoair, airwave
    
    fsps = {Z: zz, age: dblarr(nage), mstar: fltarr(nage), $
      lbol: fltarr(nage), wave: airwave, $
      flux: fltarr(npix,nage)}

    tspec = fltarr(npix)
    t = 0.0D & m = 0.0 & l = 0.0 & s = 0.0

    for ii = 0, nage-1 do begin
       readf, lun, t, m, l, s
       readf, lun, tspec
       fsps.age[ii]  = 10.0^t  ; [yr]
       fsps.mstar[ii] = 10.0^m ; [Msun]
       fsps.lbol[ii] = l
       fsps.flux[*,ii] = 3.826D33*tspec*im_light(/ang)/airwave^2 ; [Lsun/Hz]-->[erg/s/A]
    endfor
    free_lun,lun

; convert the units of the spectra as desired
    pc10 = 10.0*3.085678D18     ; =10 pc
    light = 2.99792458D18       ; [Angstrom/s]
    if keyword_set(flambda) then fsps.flux = fsps.flux/(4.0*!dpi*pc10^2)
    if keyword_set(fnu) then for ii = 0, nage-1 do $
      fsps.flux[*,ii] = fsps.flux[*,ii]/(4.0*!dpi*pc10^2)*fsps.wave^2/light
       
    if keyword_set(abmag) then begin
       for ii = 0, nage-1 do begin
          fsps.flux[*,ii] = fsps.flux[*,ii]/(4.0*!dpi*pc10^2)*fsps.wave^2/light
          gd = where(fsps.flux[*,ii] gt 0.0,ngd)
          if (ngd ne 0) then fsps.flux[gd,ii] = -2.5*alog10(fsps.flux[gd,ii])-48.6
       endfor
    endif
    
return, fsps
end
