;+
; NAME:
;   BUILD_MARASTON05_SSP
;
; PURPOSE:
;   Generate a grid of Maraston+05 simple stellar populations (SSPs)
;   to be used by ISEDFIT_MONTEGRIDS.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   kroupa - use the Kroupa+01 IMF (default is the Salpeter+05 IMF) 
;   doitall - generate all the combinations of IMFs and SSPs
;
; OUTPUTS:
;   An information structure is written to getenv('ISEDFIT_SSP_DIR')
;   and the models themselves are written to the appropriate
;   'maraston05' subdirectory.
;
; COMMENTS:
;   Adopt the 'red horizontal branch' (rhb) morphology models. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 23, UCSD
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

pro build_maraston05_ssp, kroupa=kroupa, doitall=doitall

    if keyword_set(doitall) then begin
       build_maraston05_ssp
       build_maraston05_ssp, /kroupa
       return
    endif
    
    splog, 'Building the MARASTON05 SSPs'
    outpath = getenv('ISEDFIT_SSP_DIR')+'/'

    hbmorph = 'rhb' ; red horizontal branch morphology
    if keyword_set(kroupa) then imfstr = 'kroupa01' else $
      imfstr = 'salp'

    ssppath = outpath+'maraston05/'
    if file_test(ssppath,/dir) eq 0 then file_mkdir, ssppath

    const = 1D/(6.626D-27*2.9979246D18) ; [erg*Angstrom]
    dist = 10.0*3.085678D18             ; fiducial distance [10 pc in cm]
    
; estimate the instrumental velocity dispersion; for BaSeL assume 20 A 
;   light = 2.99792458D5        ; speed of light [km/s]
;   fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
;   inst_vsigma = 20.0/5500.0/fwhm2sig*light 
    inst_vsigma = 450.0         ; [km/s]

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
       ssp.Zmetal = mara.Z

; compute the number of hydrogen-ionizing photons
       for jj = 0, nage-1 do ssp.nlyc[jj] = alog10(const*im_integral(ssp.wave,$
         1D*ssp.wave*ssp.flux[*,jj],0D,912D))

; normalize the spectrum by the stellar mass, put it at a fiducial
; distance of 10 pc, and convert to erg/s
;      ssp.flux = ssp.flux/rebin(reform(ssp.mstar,1,nage),npix,nage)
       ssp.flux = ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/Msun]

;; get the r-band mass-to-light ratio
;      rmaggies = k_project_filters(k_lambda_to_edges(ssp.wave),$
;        ssp.flux,filterlist='sdss_r0.par')

       sspfile1 = 'maraston05_'+imfstr+'_Z'+string(ssp.Zmetal,format='(G0.0)')+'.fits'
       im_mwrfits, ssp, ssppath+sspfile1, /clobber

       Z[iZ] = ssp.Zmetal
       sspfile[iZ] = sspfile1
    endfor

; write out an information structure
    info = {$
      imf:               imfstr,$
      Zmetal:                 Z,$
      inst_vsigma:  inst_vsigma,$ ; [km/s]
      sspfile:          sspfile+'.gz'}

    infofile = outpath+'info_maraston05_'+imfstr+'.fits'
    im_mwrfits, info, infofile, /clobber
    
return
end
