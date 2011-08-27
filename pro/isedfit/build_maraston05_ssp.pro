;+
; NAME:
;   BUILD_MARASTON05_SSP
;
; PURPOSE:
;   Generate a grid of Maraston+05 simple stellar populations (SSPs)
;   to be used by BUILD_ISEDFIT_SFHGRID.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   kroupa - use the Kroupa+01 IMF (default is the Salpeter+05 IMF) 
;
; OUTPUTS:
;   An information structure is written to getenv('ISEDFIT_SSP_DIR')
;   and the models themselves are written to the 'maraston05'
;   subdirectory.
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

pro build_maraston05_ssp, kroupa=kroupa

    splog, 'Building the MARASTON05 SSPs'

    ssppath = getenv('ISEDFIT_SSP_DIR')+'/'
    outpath = ssppath+'maraston05/'

    hbmorph = 'rhb' ; red horizontal branch morphology
    if keyword_set(kroupa) then imfstr = 'kroupa01' else $
      imfstr = 'salp'

    dist = 10.0*3.085678D18 ; fiducial distance [10 pc in cm]
    
; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out; do not use the two most extreme metallicities
; because they lack full age coverage
    Zstr = reverse(['z004','z002','z001','z0001'])
;   Zstr = reverse(['z007','z004','z002','z001','z0001','z10m4'])
    nZ = n_elements(Zstr)

    Z = fltarr(nZ)
    sspfile = strarr(nZ)
    for iZ = 0, nZ-1 do begin 
       mara = im_read_maraston05(metallicity=Zstr[iZ],$
         kroupa=kroupa,hbmorph=hbmorph)
       nage = n_elements(mara.age)
       npix = n_elements(mara.wave)
       ssp = init_isedfit_ssp(nage=nage,npix=npix)

       ssp.age = mara.age
       ssp.wave = mara.wave
       ssp.flux = mara.flux
       ssp.mstar = mara.mstar
       ssp.Z = mara.Z

; normalize the spectrum by the stellar mass, put it at a fiducial
; distance of 10 pc, and convert to erg/s
;      ssp.flux = ssp.flux/rebin(reform(ssp.mstar,1,nage),npix,nage)
       ssp.flux = ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/Msun]

       sspfile1 = 'maraston05_'+imfstr+'_'+Z2string(ssp.Z)+'.fits'
       im_mwrfits, ssp, outpath+sspfile1, /clobber

       Z[iZ] = ssp.Z
       sspfile[iZ] = sspfile1
    endfor

; write out an information structure
    info = {$
      imf:               imfstr,$
      Z:                      Z,$
      sspfile:          sspfile+'.gz'}

    infofile = ssppath+'info_maraston05_'+imfstr+'.fits'
    im_mwrfits, info, infofile, /clobber
    
return
end
