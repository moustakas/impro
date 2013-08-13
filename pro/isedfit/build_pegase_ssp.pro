;+
; NAME:
;   BUILD_PEGASE_SSP
;
; PURPOSE:
;   Generate a grid of PEGASE simple stellar populations (SSPs)
;   to be used by ISEDFIT_MONTEGRIDS.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   kroupa - use the Kroupa+01 IMF (default is the Salpeter IMF) 
;   cosmic_imf - generate SSPs for my PRIMUS cosmic-IMF project 
;   doitall - generate all the combinations of IMFs and SSPs
;
; OUTPUTS:
;   An information structure is written to getenv('ISEDFIT_SSP_DIR')
;   and the models themselves are written to the appropriate 'pegase'
;   subdirectory. 
;
; COMMENTS:
;   The code CREATE_PEGASE_SSPS needs to have been run.  The SSPs have
;   been built using the BaSeL stellar library (without emission
;   lines) with the additional free parameters in Pegase set to their 
;   fiducial recommended values.  The highest stellar metallicity
;   value (Z=0.01) is ignored.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Mar 29, UCSD
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

pro build_pegase_ssp, kroupa=kroupa, cosmic_imf=cosmic_imf, doitall=doitall

    if keyword_set(doitall) then begin
       build_pegase_ssp
       build_pegase_ssp, /kroupa
       return
    endif

    splog, 'Building the PEGASE SSPs'
    outpath = getenv('ISEDFIT_SSP_DIR')+'/'
    pegpath = getenv('PEGASE_HR_DIR')+'/SSPs/'
    
    imfstr = 'salp'
    if keyword_set(kroupa) then imfstr = 'kroupa01'
    if keyword_set(cosmic_imf) then begin
       alpha2 = primus_cosmicimf_slope() ; high-mass slope
       imfstr = 'cosmicimf_'+string(alpha2,format='(F4.2)')
    endif

    ssppath = outpath+'pegase/'
    if file_test(ssppath,/dir) eq 0 then file_mkdir, ssppath

    const = 1D/(6.626D-27*2.9979246D18) ; [erg*Angstrom]
    dist = 10.0*3.085678D18             ; fiducial distance [10 pc in cm]
    
; estimate the instrumental velocity dispersion; for BaSeL assume 20 A 
;   light = 2.99792458D5        ; speed of light [km/s]
;   fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
;   inst_vsigma = 20.0/5500.0/fwhm2sig*light 
    inst_vsigma = 450.0         ; [km/s]

; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out
    Zstr = ['0.0001','0.0004','0.004','0.008','0.02','0.05']
;   Zstr = ['0.0001','0.0004','0.004','0.008','0.02','0.05','0.1']
    nZ = n_elements(Zstr)

    Z = fltarr(nZ)
    sspfile = strarr(nZ)
    for ii = 0, n_elements(imfstr)-1 do begin ; need this for /COSMIC_IMF
       for iZ = 0, nZ-1 do begin
          infile = pegpath+'SSP_'+imfstr[ii]+'_Z'+Zstr[iZ]+'.fits'
          pegase = im_read_peg(infile)
          nage = n_elements(pegase.age)
          npix = n_elements(pegase[0].wave)
          ssp = init_isedfit_ssp(nage=nage,npix=npix)

          ssp.age = pegase.age*1D6 ; [yr]
          ssp.wave = pegase[0].wave
          ssp.flux = pegase.flux
          ssp.mstar = pegase.mstar
          ssp.Zmetal = Zstr[iZ]

; compute the number of hydrogen-ionizing photons
          for jj = 0, nage-1 do ssp.nlyc[jj] = alog10(const*im_integral(ssp.wave,$
            1D*ssp.wave*ssp.flux[*,jj],0D,912D))

; normalize the spectrum by the stellar mass, put it at a fiducial
; distance of 10 pc, and convert to erg/s
;         ssp.flux = ssp.flux/rebin(reform(ssp.mstar,1,nage),npix,nage)
          ssp.flux = ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/Msun]

;; get the r-band mass-to-light ratio
;         rmaggies = k_project_filters(k_lambda_to_edges(ssp.wave),$
;           ssp.flux,filterlist='sdss_r0.par')

          sspfile1 = 'pegase_'+imfstr[ii]+'_Z'+string(ssp.Zmetal,format='(G0.0)')+'.fits'
          im_mwrfits, ssp, ssppath+sspfile1, /clobber

          Z[iZ] = ssp.Zmetal
          sspfile[iZ] = sspfile1
       endfor

; write out an information structure
       info = {$
         imf:               imfstr[ii],$
         Zmetal:                     Z,$
         inst_vsigma:  inst_vsigma,$ ; [km/s]
         sspfile:          sspfile+'.gz'}

       infofile = outpath+'info_pegase_'+imfstr[ii]+'.fits'
       im_mwrfits, info, infofile, /clobber
    endfor
    
return
end
