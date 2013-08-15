;+
; NAME:
;   BUILD_BASTI_SSP
;
; PURPOSE:
;   Generate a grid of BaSTI simple stellar populations (SSPs) to be
;   to be used by ISEDFIT_MONTEGRIDS.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   enhanced - write out the alpha-enhanced models (default is to
;     write out the solar-scaled models)
;   doitall - generate all the combinations of IMFs and SSPs
;
; OUTPUTS:
;   An information structure is written to getenv('ISEDFIT_SSP_DIR')
;   and the models themselves are written to the appropriate 'basti'
;   subdirectory. 
;
; COMMENTS:
;   See http://albione.oa-teramo.inaf.it for additional relevant
;   details.  The models downloaded are as of 2010-Nov.
; 
;   Only the Kroupa+01 IMF is available. 
;   
;   Uses the recommended default models with mass-loss parameter
;   eta=0.4.  Also uses the low-resolution theoretical models, which
;   have greater wavelength coverage, although in principle the
;   high-resolution theoretical models could be incorporated as well.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 24, UCSD
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

pro build_basti_ssp, enhanced=enhanced, doitall=doitall

    if keyword_set(doitall) then begin
; Kroupa+01 default (solar-scaled) and alpha-enhanced 
       build_basti_ssp
       build_basti_ssp, /enhanced
       return
    endif

    splog, 'Building the BaSTI SSPs'
    outpath = getenv('ISEDFIT_SSP_DIR')+'/'

    imfstr = 'kroupa01'
    if keyword_set(enhanced) then enh = 'ae' else enh = 'ss'

    ssppath = outpath+'basti_'+enh+'/'
    if file_test(ssppath,/dir) eq 0 then file_mkdir, ssppath

    const = 1D/(6.626D-27*2.9979246D18) ; [erg*Angstrom]
    dist = 3.085678D19 ; fiducial distance [10 pc in cm]
    
; estimate the instrumental velocity dispersion; assume 20 A FWHM
; resolution for the low-resolution models (need to double-check
; this!); the relevant reference is Percival+09
;   light = 2.99792458D5        ; speed of light [km/s]
;   fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
;   inst_vsigma = 20.0/5500.0/fwhm2sig*light 
    inst_vsigma = 450.0         ; [km/s]
    
; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out; do not use the full range of metallicities since
; the spectra are somewhat incomplete; also, cut the models at 14.5
; Gyr because some 15 Gyr models are missing (e.g., z803,
; sss_agb.t615000) 
    Zstr = ['z103','z203','z403','z803','z102','zsun','z302','z402']
    nZ = n_elements(Zstr)

    Z = fltarr(nZ)
    sspfile = strarr(nZ)
    for iZ = 0, nZ-1 do begin 
       basti = im_read_basti(metallicity=Zstr[iZ],enhanced=enhanced,/silent)
       keep = where(basti.age lt 15D9,nage)
       npix = n_elements(basti.wave)
       ssp = init_isedfit_ssp(nage=nage,npix=npix)

       ssp.age = basti.age[keep]
       ssp.mstar = basti.mstar[keep]
       ssp.wave = basti.wave
       ssp.flux = basti.flux[*,keep]
       ssp.Zmetal = basti.Z

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
       sspfile1 = 'basti_'+enh+'_'+imfstr+'_Z'+string(ssp.Zmetal,format='(G0.0)')+'.fits'
       im_mwrfits, ssp, ssppath+sspfile1, /clobber

       Z[iZ] = ssp.Zmetal
       sspfile[iZ] = sspfile1
    endfor

; write out an information structure
    info = {$
      imf:               imfstr,$
      Zmetal:                 Z,$
      inst_vsigma:  inst_vsigma,$ ; [km/s]
      sspfile:     sspfile+'.gz'}

    infofile = outpath+'info_basti_'+enh+'_'+imfstr+'.fits'
    im_mwrfits, info, infofile, /clobber
    
return
end
