;+
; NAME:
;   BUILD_BC03_SSP
;
; PURPOSE:
;   Generate a grid of BC03 simple stellar populations (SSPs) to be
;   to be used by ISEDFIT_MONTEGRIDS.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   chabrier - use the Chabrier+03 IMF (default is the Salpeter+55 IMF) 
;   basel - use the low-resolution BC03 models based on the BaSeL
;     stellar library (default is to use the high-resolution STELIB
;     models) 
;   doitall - generate all the combinations of IMFs and SSPs
;
; OUTPUTS:
;   An information structure is written to getenv('ISEDFIT_SSP_DIR')
;   and the models themselves are written to the appropriate 'bc03'
;   subdirectory. 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan  22, UCSD
;   jm13aug01siena - compute the number of Lyman continuum photons;
;     added DOITALL keyword; additional minor tweaks
;
; Copyright (C) 2011, 2013, John Moustakas
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

pro build_bc03_ssp, chabrier=chabrier, basel=basel, doitall=doitall

    if keyword_set(doitall) then begin
; Salpeter high and low resolution
       build_bc03_ssp, chabrier=0, basel=0
       build_bc03_ssp, chabrier=0, basel=1
; Chabrier high and low resolution
       build_bc03_ssp, chabrier=1, basel=0
       build_bc03_ssp, chabrier=1, basel=1
       return
    endif
    
    splog, 'Building the BC03 SSPs'
    outpath = getenv('ISEDFIT_SSP_DIR')+'/'

    bc03path = getenv('bc03_dir')+'/models/Padova1994/'
    if keyword_set(chabrier) then begin
       bc03path = bc03path+'chabrier/' 
       imfstr = 'chab'
       salpeter = 0
    endif else begin
       bc03path = bc03path+'salpeter/'
       imfstr = 'salp'
       salpeter = 1
    endelse

    if keyword_set(basel) then begin
       resstr = 'lr'
       lsuffix = 'basel'
    endif else begin
       resstr = 'hr'
       lsuffix = 'stelib'
    endelse

    ssppath = outpath+'bc03_'+lsuffix+'/'
    if file_test(ssppath,/dir) eq 0 then file_mkdir, ssppath

    lsun = 3.826D33                     ; [erg/s]
    const = 1D/(6.626D-27*2.9979246D18) ; [erg*Angstrom]
    dist = 10.0*3.085678D18             ; fiducial distance [10 pc in cm]

; estimate the instrumental velocity dispersion; for the STELIB+BaSeL
; SSPs assume 3 A in the optical (ignore the change in resolution at
; longer and shorter wavelengths); for BaSeL assume 20 A
;   light = 2.99792458D5        ; speed of light [km/s]
;   fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
    if keyword_set(basel) then begin 
;      inst_vsigma = 20.0/5500.0/fwhm2sig*light 
       inst_vsigma = 450.0       ; [km/s]
    endif else begin
;      inst_vsigma = 3.0/5500.0/fwhm2sig*light 
       inst_vsigma = 70.0        ; [km/s]
    endelse

; metallicity grid
    Z = [0.0004,0.004,0.008,0.02,0.05]
    Zstr = 'm'+['32','42','52','62','72']
    nZ = n_elements(Z)

; input/output file names
    bc03file = 'bc2003_'+resstr+'_'+Zstr+'_'+imfstr+'_ssp.ised'
    sspfile = 'bc03_'+lsuffix+'_'+imfstr+'_Z'+string(Z,format='(G0.0)')+'.fits'

; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out
    for iZ = 0, nZ-1 do begin 
       bc03 = im_read_bc03(isedpath=bc03path,isedfile=bc03file[iZ],$
         bc03_extras=ext,/array_extras,salpeter=salpeter,lr=keyword_set(basel))
       nage = n_elements(bc03.age)
       npix = n_elements(bc03.wave)
       ssp = init_isedfit_ssp(nage=nage,npix=npix)

       ssp.age = bc03.age
       ssp.wave = bc03.wave
       ssp.flux = bc03.flux
       ssp.mstar = ext.m_ ; stellar mass [Msun]
       ssp.Zmetal = Z[iZ]

; compute the number of hydrogen-ionizing photons
       for jj = 0, nage-1 do ssp.nlyc[jj] = alog10(const*im_integral(ssp.wave,$
         lsun*ssp.wave*ssp.flux[*,jj],0D,912D))

; put the model at a fiducial distance of 10 pc, and convert to erg/s 
       ssp.flux = lsun*ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/Msun]
;      ssp.flux = ssp.flux/rebin(reform(ssp.mstar,1,nage),npix,nage)

;; get the r-band mass-to-light ratio
;       rmaggies = k_project_filters(k_lambda_to_edges(ssp.wave),$
;         ssp.flux,filterlist='sdss_r0.par')

; write out       
       im_mwrfits, ssp, ssppath+sspfile[iZ], /clobber
    endfor

; write out an information structure
    info = {$
      imf:               imfstr,$
      Zmetal:                 Z,$
      inst_vsigma:  inst_vsigma,$ ; [km/s]
      sspfile:          sspfile+'.gz'}

    infofile = outpath+'info_bc03_'+lsuffix+'_'+imfstr+'.fits'
    im_mwrfits, info, infofile, /clobber
    
return
end
