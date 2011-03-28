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
;   basti - read the BaSTI isochrones (default is to use Padova) 
;   hires - read the MILES+BaSeL stellar library, which is
;     high-resolution in the optical (default is to use just the 
;     low-resolution BaSeL stellar library), but note that the
;     MILES+BaSeL library is incomplete at low metallicity and at
;     young ages (see Vazdekis+10).
;   kroupa - read the Kroupa+01 IMF files (default is to read the 
;     Salpeter ones)
;   chabrier - read the Chabrier+03 IMF files
;
; OUTPUTS:
;   fsps - output data structure
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
  hires=hires, kroupa=kroupa, chabrier=chabrier 

    ssppath = getenv('IM_DATA_DIR')+'/synthesis/fsps/SSP/'

; defaults
    if keyword_set(hires) then lib = 'MILES+BaSeL' else lib = 'BaSeL' ; stellar library
    if keyword_set(basti) then iso = 'BaSTI' else iso = 'Padova' ; isochrones

    imf = 'Salpeter'
    if keyword_set(kroupa) then imf = 'Kroupa'
    if keyword_set(chabrier) then imf = 'Chabrier'

; metallicity
    case strlowcase(iso) of
       'padova': begin
          if (n_elements(metallicity) eq 0) then metallicity = 'Z0.0190'
          if keyword_set(lowres) then begin ; BaSeL+MILES
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
          if keyword_set(hires) then begin ; BaSeL+MILES
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

; read the wavelength array
    wavefile = ssppath+strlowcase(lib)+'.lambda'
    if (file_test(wavefile) eq 0) then begin
       splog, 'Wavelength file '+wavfile+' not found!'
       return, -1
    endif
    splog, 'Reading '+wavefile
    readcol, wavefile, wave, format='F', /silent
    npix = n_elements(wave)
    
; read the SSP using Conroy's code READ_SPEC
    splog, 'Reading '+sspfile
    openr, lun, sspfile, /get_lun
    char = '#' ; burn the header
    while (strmid(char,0,1) eq '#') do readf, lun, char

    nage = long(char)
    fsps = {Z: zz, age: dblarr(nage), mstar: fltarr(nage), $
      lbol: fltarr(nage), wave: wave, $
      flux: fltarr(npix,nage)}

    tspec = fltarr(npix)
    t = 0.0D & m = 0.0 & l = 0.0 & s = 0.0

    for ii = 0, nage-1 do begin
       readf, lun, t, m, l, s
       readf, lun, tspec
       fsps.age[ii]  = 10.0^t  ; [yr]
       fsps.mstar[ii] = 10.0^m ; [Msun]
       fsps.lbol[ii] = l
       fsps.flux[*,ii] = tspec*im_light(/ang)/wave^2 ; [Lsun/Hz]-->[erg/s/A]
    endfor
    free_lun,lun

return, fsps
end
