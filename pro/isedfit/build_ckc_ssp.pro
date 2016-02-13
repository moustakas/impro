;+
; NAME:
;   BUILD_CKC_SSP
;
; PURPOSE:
;   Generate a grid of CKC simple stellar populations (SSPs)
;   to be used by ISEDFIT_MONTEGRIDS.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   An information structure is written to getenv('ISEDFIT_SSP_DIR')
;   and the models themselves are written to the appropriate 'ckc'
;   subdirectory. 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2016 Feb 2, Siena College
;
; Copyright (C) 2016 John Moustakas
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

pro build_ckc_ssp

    splog, 'Building the CKC14z SSPs'
    outpath = getenv('ISEDFIT_SSP_DIR')+'/'

    ckc_ver = '14z' ; hard-coded version number!
    imfstr = 'kroupa01'

    ssppath = outpath+'ckc'+ckc_ver+'/'
    if file_test(ssppath,/dir) eq 0 then file_mkdir, ssppath
    
    const = 1D/(6.626D-27*2.9979246D18) ; [erg*Angstrom]
    dist = 10.0*3.085678D18             ; fiducial distance [10 pc in cm]
    light = 2.99792458D5                ; [km/s]

; median pixel size of the models is 15 km/s, corresponding to roughly 0.7 A
; FWHM resolution at 5500 A
;  ii = im_read_fsps(/ckc)    
;  print, djs_median(alog10(ii.wave)-shift(alog10(ii.wave),1))*alog(10)*3E5
;  print, 0.7D/5500/2.35*3E5
    ckc_velpixsize = 15D        ; [km/s]

; resample the models to constant log10-lambda rest-frame wavelengths (in
; *air*); this also specifies the instrumental velocity dispersion
    inst_vsigma = 50D                    ; [km/s]
    pixsize = inst_vsigma/light/alog(10) ; [pixel size in log-10 A]
    minwave = alog10(950D)               ; minimum wavelength [log10-A]
    maxwave = alog10(55000D)             ; maximum wavelength [log10-A]
    npix = round((maxwave-minwave)/pixsize+1L)
    restwave = minwave+dindgen(npix)*pixsize ; constant log-10 spacing

; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out
    Zstr = ['Z0.0003','Z0.0006','Z0.0012','Z0.0025',$
      'Z0.0049','Z0.0096','Z0.0190','Z0.0300']
    nZ = n_elements(Zstr)

    Z = fltarr(nZ)
    sspfile = strarr(nZ)
    for iZ = 0, nZ-1 do begin 
       fsps = im_read_fsps(metallicity=Zstr[iZ],/ckc) ; note! no /VACUUM
       nage = n_elements(fsps.age)
;      npix = n_elements(fsps.wave)
       ssp = init_isedfit_ssp(nage=nage,npix=npix)

       ssp.age = fsps.age
       ssp.mstar = fsps.mstar
       ssp.Zmetal = fsps.Z

; not sure why, but we need to sort in wavelength
;      srt = sort(fsps.wave)
;      wave = fsps.wave[srt]
;      flux = fsps.flux[srt,*]

; compute the number of hydrogen-ionizing photons
       for jj = 0, nage-1 do ssp.nlyc[jj] = alog10(const*im_integral(wave,$
         1D*wave*flux[*,jj],0D,912D))

; interpolate (this should be smarter!)
       ssp.wave = 10D^restwave
;      for jj = 0, nage-1 do ssp.flux[*,jj] = rebin_spectrum(flux[*,jj],$
;        alog10(wave),restwave)
       ssp.flux = interpolate(flux,findex(alog10(wave),restwave),lindgen(nage),/grid)

; put each spectrum at a fiducial distance of 10 pc, and convert to erg/s; do
; not normalize by the stellar mass here
;      ssp.flux = ssp.flux/rebin(reform(ssp.mstar,1,nage),npix,nage)
       ssp.flux = ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/Msun]

;; get the r-band mass-to-light ratio
;       rmaggies = k_project_filters(k_lambda_to_edges(ssp.wave),$
;         ssp.flux,filterlist='sdss_r0.par')

; write out       
       sspfile1 = 'ckc'+ckc_ver+'_'+imfstr+'_Z'+$
         string(ssp.Zmetal,format='(G0.0)')+'.fits'
       im_mwrfits, ssp, ssppath+sspfile1, /clobber

       Z[iZ] = ssp.Zmetal
       sspfile[iZ] = sspfile1
    endfor

; write out an information structure
    info = {$
      imf:               imfstr,$
      Zmetal:                 Z,$
      inst_vsigma:  inst_vsigma,$ ; [km/s]
      sspfile:      sspfile+'.gz'}

    infofile = outpath+'info_ckc'+ckc_ver+'_'+imfstr+'.fits'
;   infofile = outpath+'info_fsps_'+imfstr+'.fits'
    im_mwrfits, info, infofile, /clobber
    
return
end
