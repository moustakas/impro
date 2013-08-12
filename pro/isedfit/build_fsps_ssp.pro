;+
; NAME:
;   BUILD_FSPS_SSP
;
; PURPOSE:
;   Generate a grid of FSPS simple stellar populations (SSPs)
;   to be used by ISEDFIT_MONTEGRIDS.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   kroupa - use the Kroupa+01 IMF (default is the Salpeter IMF) 
;   chabrier - use the Chabrier+03 IMF
;   miles - read the MILES+BaSeL SSPs in lieu of the pure BaSeL SSPs
;   doitall - generate all the combinations of IMFs and SSPs
;
; OUTPUTS:
;   An information structure is written to getenv('ISEDFIT_SSP_DIR')
;   and the models themselves are written to the appropriate 'fsps'
;   subdirectory. 
;
; COMMENTS:
;   If there are any updates to the FSPS models all you have to do is
;   run this code with the /DOITALL keyword.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 23, UCSD
;   jm13aug01siena - added the FSPS version number and the SSP names
;     to the output files since FSPS is under development; compute the
;     number of Lyman continuum photons; added DOITALL keyword
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

pro build_fsps_ssp, kroupa=kroupa, chabrier=chabrier, miles=miles, doitall=doitall

    if keyword_set(doitall) then begin
; Salpeter with and without MILES
       build_fsps_ssp, kroupa=0, chabrier=0, miles=0
       build_fsps_ssp, kroupa=0, chabrier=0, miles=1
; Kroupa with and without MILES
       build_fsps_ssp, kroupa=1, chabrier=0, miles=0
       build_fsps_ssp, kroupa=1, chabrier=0, miles=1
; Chabrier with and without MILES
       build_fsps_ssp, kroupa=0, chabrier=1, miles=0
       build_fsps_ssp, kroupa=0, chabrier=1, miles=1
       return
    endif
    
    splog, 'Building the FSPS SSPs'
    outpath = getenv('ISEDFIT_SSP_DIR')+'/'

    fsps_ver = 'v2.4'           ; hard-coded version number!
    imfstr = 'salp'
    if keyword_set(kroupa) then imfstr = 'kroupa01'
    if keyword_set(chabrier) then imfstr = 'chab'
    if keyword_set(miles) then sspstr = 'miles' else sspstr = 'basel'

    ssppath = outpath+'fsps_'+fsps_ver+'_'+sspstr+'/'
    if file_test(ssppath,/dir) eq 0 then file_mkdir, ssppath
    
    const = 1D/(6.626D-27*2.9979246D18) ; [erg*Angstrom]
    dist = 10.0*3.085678D18             ; fiducial distance [10 pc in cm]

; estimate the instrumental velocity dispersion; for the MILES+BaSeL
; SSPs assume 3 A in the optical (ignore the change in resolution at
; longer and shorter wavelengths); for BaSeL assume 20 A
;   light = 2.99792458D5        ; speed of light [km/s]
;   fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
    if keyword_set(miles) then begin 
;      inst_vsigma = 3.0/5500.0/fwhm2sig*light 
       inst_vsigma = 70.0        ; [km/s]
    endif else begin
;      inst_vsigma = 20.0/5500.0/fwhm2sig*light 
       inst_vsigma = 450.0       ; [km/s]
    endelse
    
; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out
    if keyword_set(miles) then begin
       Zstr = ['Z0.0008','Z0.0031','Z0.0096','Z0.0190','Z0.0300'] ; MILES+BaSeL
    endif else begin
       Zstr = ['Z0.0002','Z0.0003','Z0.0004','Z0.0005','Z0.0006',$ ; BaSeL
         'Z0.0008','Z0.0010','Z0.0012','Z0.0016','Z0.0020','Z0.0025',$
         'Z0.0031','Z0.0039','Z0.0049','Z0.0061','Z0.0077','Z0.0096',$
         'Z0.0120','Z0.0150','Z0.0190','Z0.0240','Z0.0300']
    endelse
    nZ = n_elements(Zstr)

    Z = fltarr(nZ)
    sspfile = strarr(nZ)
    for iZ = 0, nZ-1 do begin 
       fsps = im_read_fsps(metallicity=Zstr[iZ],$ ; Padova/BaSeL unless /MILES
         kroupa=kroupa,chabrier=chabrier,miles=miles)
       nage = n_elements(fsps.age)
       npix = n_elements(fsps.wave)
       ssp = init_isedfit_ssp(nage=nage,npix=npix)

       ssp.age = fsps.age
       ssp.flux = fsps.flux
       ssp.wave = fsps.wave
       ssp.mstar = fsps.mstar
       ssp.Zmetal = fsps.Z

; compute the number of hydrogen-ionizing photons
       for jj = 0, nage-1 do ssp.nlyc[jj] = alog10(const*im_integral(ssp.wave,$
         1D*ssp.wave*ssp.flux[*,jj],0D,912D))

; normalize the spectrum by the stellar mass, put it at a fiducial
; distance of 10 pc, and convert to erg/s
;      ssp.flux = ssp.flux/rebin(reform(ssp.mstar,1,nage),npix,nage)
       ssp.flux = ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/Msun]

;; get the r-band mass-to-light ratio
;       rmaggies = k_project_filters(k_lambda_to_edges(ssp.wave),$
;         ssp.flux,filterlist='sdss_r0.par')

; write out       
       sspfile1 = 'fsps_'+fsps_ver+'_'+sspstr+'_'+imfstr+'_Z'+$
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

    infofile = outpath+'info_fsps_'+fsps_ver+'_'+sspstr+'_'+imfstr+'.fits'
;   infofile = outpath+'info_fsps_'+imfstr+'.fits'
    im_mwrfits, info, infofile, /clobber
    
return
end
