;+
; NAME:
;   BUILD_BC03_SSP
;
; PURPOSE:
;   Generate a grid of BC03 simple stellar populations (SSPs) to be
;   used by BUILD_ISEDFIT_SFHGRID.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   chabrier - use the Chabrier+03 IMF (default is the Salpeter+55 IMF) 
;   lowres - use the low-resolution BC03 models (default is the
;     high-resolution models)
;
; OUTPUTS:
;   An information structure is written to getenv('ISEDFIT_SSP_DIR')
;   and the models themselves are written to the 'bc03' subdirectory.
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan  22, UCSD
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

pro build_bc03_ssp, chabrier=chabrier, lowres=lowres

    splog, 'Building the BC03 SSPs'

    ssppath = getenv('ISEDFIT_SSP_DIR')+'/'
    outpath = ssppath+'bc03/'
    bc03path = getenv('bc03_dir')+'/models/Padova1994/'

    if keyword_set(lowres) then begin
       resstr = 'lr'
       lsuffix = '_lowres'
       outpath = ssppath+'bc03_lowres/'
    endif else begin
       resstr = 'hr'
       lsuffix = ''
    endelse

    if keyword_set(chabrier) then begin
       bc03path = bc03path+'chabrier/' 
       imfstr = 'chab'
       salpeter = 0
    endif else begin
       bc03path = bc03path+'salpeter/'
       imfstr = 'salp'
       salpeter = 1
    endelse

    lsun = 3.826D33         ; [erg/s]
    dist = 10.0*3.085678D18 ; fiducial distance [10 pc in cm]
    
; metallicity grid
    Z = [0.0004,0.004,0.008,0.02,0.05]
    Zstr = 'm'+['32','42','52','62','72']
    nZ = n_elements(Z)

; input/output file names
    bc03file = 'bc2003_'+resstr+'_'+Zstr+'_'+imfstr+'_ssp.ised'
    sspfile = 'bc03'+lsuffix+'_'+imfstr+'_'+Z2string(Z)+'.fits'

; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out
    for iZ = 0, nZ-1 do begin 
       bc03 = im_read_bc03(isedpath=bc03path,isedfile=bc03file[iZ],$
         bc03_extras=ext,/array_extras,salpeter=salpeter,lr=keyword_set(lowres))

       nage = n_elements(bc03.age)
       npix = n_elements(bc03.wave)
       ssp = init_isedfit_ssp(nage=nage,npix=npix)

       ssp.age = bc03.age
       ssp.wave = bc03.wave
       ssp.flux = bc03.flux
       ssp.mstar = ext.m_ ; stellar mass [Msun]
       ssp.Z = Z[iZ]

; put the model at a fiducial distance of 10 pc, and convert to erg/s 
       ssp.flux = lsun*ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/Msun]
;      ssp.flux = ssp.flux/rebin(reform(ssp.mstar,1,nage),npix,nage)

       im_mwrfits, ssp, outpath+sspfile[iZ], /clobber
    endfor

; write out an information structure
    info = {$
      imf:               imfstr,$
      Z:                      Z,$
      sspfile:          sspfile+'.gz'}

    infofile = ssppath+'info_bc03'+lsuffix+'_'+imfstr+'.fits'
    im_mwrfits, info, infofile, /clobber
    
return
end
