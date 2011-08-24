;+
; NAME:
;   BUILD_FSPS_SSP
;
; PURPOSE:
;   Generate a grid of FSPS simple stellar populations (SSPs)
;   to be used by BUILD_ISEDFIT_SFHGRID.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   kroupa - use the Kroupa+01 IMF (default is the Salpeter IMF) 
;   chabrier - use the Chabrier+03 IMF
;
; OUTPUTS:
;   An information structure is written to getenv('ISEDFIT_SSP_DIR')
;   and the models themselves are written to the 'fsps' subdirectory.
;
; COMMENTS:
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

pro build_fsps_ssp, kroupa=kroupa, chabrier=chabrier

    splog, 'Building the FSPS SSPs'

    ssppath = getenv('ISEDFIT_SSP_DIR')+'/'
    outpath = ssppath+'fsps/'

    imfstr = 'salp'
    if keyword_set(kroupa) then imfstr = 'kroupa01'
    if keyword_set(chabrier) then imfstr = 'chab'

    dist = 10.0*3.085678D18 ; fiducial distance [10 pc in cm]
    
; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out
;   Zstr = ['Z0.0008','Z0.0031','Z0.0096','Z0.0190','Z0.0300'] ; MILES+BaSeL
    Zstr = ['Z0.0002','Z0.0003','Z0.0004','Z0.0005','Z0.0006',$
      'Z0.0008','Z0.0010','Z0.0012','Z0.0016','Z0.0020','Z0.0025',$
      'Z0.0031','Z0.0039','Z0.0049','Z0.0061','Z0.0077','Z0.0096',$
      'Z0.0120','Z0.0150','Z0.0190','Z0.0240','Z0.0300']
    nZ = n_elements(Zstr)

    Z = fltarr(nZ)
    sspfile = strarr(nZ)
    for iZ = 0, nZ-1 do begin 
       fsps = im_read_fsps(metallicity=Zstr[iZ],$ ; Padova/BaSeL
         kroupa=kroupa,chabrier=chabrier)
       nage = n_elements(fsps.age)
       npix = n_elements(fsps.wave)
       ssp = init_isedfit_ssp(nage=nage,npix=npix)

       ssp.age = fsps.age
       ssp.wave = fsps.wave
       ssp.flux = fsps.flux
       ssp.mstar = fsps.mstar
       ssp.Z = fsps.Z

; normalize the spectrum by the stellar mass, put it at a fiducial
; distance of 10 pc, and convert to erg/s
;      ssp.flux = ssp.flux/rebin(reform(ssp.mstar,1,nage),npix,nage)
       ssp.flux = ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/Msun]

       sspfile1 = 'fsps_'+imfstr+'_'+Z2string(ssp.Z)+'.fits'
       im_mwrfits, ssp, outpath+sspfile1, /clobber

       Z[iZ] = ssp.Z
       sspfile[iZ] = sspfile1
    endfor

; write out an information structure
    info = {$
      imf:               imfstr,$
      Z:                      Z,$
      sspfile:          sspfile+'.gz'}

    infofile = ssppath+'info_fsps_'+imfstr+'.fits'
    im_mwrfits, info, infofile, /clobber
    
return
end
