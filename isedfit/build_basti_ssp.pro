;+
; NAME:
;   BUILD_BASTI_SSP
;
; PURPOSE:
;   Generate a grid of BaSTI simple stellar populations (SSPs) to be
;   used by BUILD_ISEDFIT_SFHGRID. 
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   enhanced - write out the alpha-enhanced models
;
; OUTPUTS:
;   An information structure is written to
;   getenv('ISEDFIT_SFHGRID_DIR')+'/ssp/' and the models themselves
;   are written to the 'basti' subdirectory.
;
; COMMENTS:
;   See http://albione.oa-teramo.inaf.it for additional relevant
;   details.  The models downloaded are as of 2010-Nov.
; 
;   Only the Kroupa+01 IMF is available. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 24, UCSD
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

pro build_basti_ssp, enhanced=enhanced

    splog, 'Building the BASTI SSPs'

    imfstr = 'kroupa01'
    if keyword_set(enhanced) then enh = 'ae' else enh = 'ss'

    ssppath = getenv('ISEDFIT_SFHGRID_DIR')+'/ssp/'
    outpath = ssppath+'basti_'+enh+'/'

    dist = 10.0*3.085678D18 ; fiducial distance [10 pc in cm]
    
; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out; do not use the full range of metallicities since
; the spectra are somewhat incomplete 
    Zstr = reverse(['z103','z203','z403','z803','z102','zsun','z302','z402'])
    nZ = n_elements(Zstr)

    Z = fltarr(nZ)
    sspfile = strarr(nZ)
    for iZ = 0, nZ-1 do begin 
       basti = im_read_basti(metallicity=Zstr[iZ],enhanced=enhanced)
       nage = n_elements(basti.age)
       npix = n_elements(basti.wave)
       ssp = init_isedfit_ssp(nage=nage,npix=npix)

       ssp.age = basti.age
       ssp.wave = basti.wave
       ssp.flux = basti.flux
       ssp.mstar = basti.mstar
       ssp.Z = basti.Z

; normalize the spectrum by the stellar mass, put it at a fiducial
; distance of 10 pc, and convert to erg/s
;      ssp.flux = ssp.flux/rebin(reform(ssp.mstar,1,nage),npix,nage)
       ssp.flux = ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/Msun]

       sspfile1 = 'basti_'+enh+'_'+imfstr+'_'+Z2string(ssp.Z)+'.fits'
       im_mwrfits, ssp, outpath+sspfile1, /clobber

       Z[iZ] = ssp.Z
       sspfile[iZ] = sspfile1
    endfor

; write out an information structure
    info = {$
      imf:               imfstr,$
      Z:                      Z,$
      sspfile:          sspfile+'.gz'}

    infofile = ssppath+'info_basti_'+enh+'_'+imfstr+'.fits'
    im_mwrfits, info, infofile, /clobber
    
return
end
