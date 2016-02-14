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

function ckc_blur, oldflux, old_velscale, new_velscale
; convolve each model to a new velocity resolution
    ntemp = (size(oldflux,/dim))[1]
    vdisp = sqrt(new_velscale^2 - old_velscale^2)
    smoothing = vdisp/old_velscale  ; [pixel]
;   kernel = psf_gaussian(npix=10.0*smoothing,$
;   fwhm=2.35*smoothing,/norm,ndim=1)
    nkpix = long(4.0*ceil(smoothing))*2L+3
    klam = findgen(nkpix)-float(nkpix-1.0)/2.0
    kernel = exp(-0.5*(klam/smoothing)^2)/sqrt(2.0*!dpi)/smoothing
    kernel = kernel/total(kernel)
    newflux = oldflux
    for ii = 0, ntemp-1 do newflux[*,ii] = $
      convol(oldflux[*,ii],kernel,/edge_truncate)
return, newflux
end

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
;   ckc = im_read_fsps(/ckc)    
;   ckc_velscale = djs_median(alog10(ckc.wave)-shift(alog10(ckc.wave),1))*alog(10)*light
    ckc_velscale = 15D          ; [km/s]
    ckc_pixsize = ckc_velscale/light/alog(10)

; resample to be constant in log10-lambda (in *air*) at a lower pixel scale
    velscale = 50D                    ; [km/s]
    pixsize = velscale/light/alog(10) ; [pixel size in log-10 A]
    minwave = 500D                    ; minimum wavelength [log10-A]
    maxwave = 60000D                  ; maximum wavelength [log10-A]

;   minwave = alog10(950D)               ; minimum wavelength [log10-A]
;   maxwave = alog10(55000D)             ; maximum wavelength [log10-A]
;   npix = round((maxwave-minwave)/pixsize+1L)
;   restwave = minwave+dindgen(npix)*pixsize ; constant log-10 spacing

; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out
    Zstr = ['Z0.0003','Z0.0006','Z0.0012','Z0.0025',$
      'Z0.0049','Z0.0096','Z0.0190','Z0.0300']
    nZ = n_elements(Zstr)

    Z = fltarr(nZ)
    sspfile = strarr(nZ)
    for iZ = 0, nZ-1 do begin 
       ckc = im_read_fsps(metallicity=Zstr[iZ],/ckc) ; note! no /VACUUM
       nage = n_elements(ckc.age)
;      npix = n_elements(ckc.wave)

; convolve to a lower spectral resolution with a uniform wavelength array; we do
; this here at the top so we can get the number of pixels
       these = where(ckc.wave gt minwave*0.98 and ckc.wave lt maxwave*1.02)
       ckcwave = ckc.wave[these]
       ckcflux = ckc.flux[these,*]*1D
;      ckcwave = ckc.wave
;      ckcflux = ckc.flux
       
       cflux = ckc_blur(ckcflux,ckc_velscale,velscale)
       for jj = 0, nage-1 do begin
          flux1 = im_log_rebin(ckcwave,reform(cflux[*,jj]),$
            vsc=velscale,outwave=lnwave,minwave=minwave,maxwave=maxwave)
          if jj eq 0 then flux = fltarr(n_elements(flux1),nage)
          flux[*,jj] = flux1
       endfor
;      flux = interpolate(flux,findex(alog10(ckc.wave),wave),lindgen(nage),/grid)
       npix = n_elements(lnwave)

;      djs_plot, ckc.wave, ckc.flux[*,120], xr=[1500,5D4], /xlog, /ylog
;      djs_oplot, exp(lnwave), flux[*,120], color='green'              

       ssp = init_isedfit_ssp(nage=nage,npix=npix)

       ssp.age = ckc.age
       ssp.mstar = ckc.mstar
       ssp.Zmetal = ckc.Z
       ssp.wave = exp(lnwave)
       ssp.flux = flux
       
; compute the number of hydrogen-ionizing photons
       for jj = 0, nage-1 do ssp.nlyc[jj] = alog10(const*im_integral(ckc.wave,$
         1D*ckc.wave*ckc.flux[*,jj],0D,912D))

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
      inst_vsigma:     velscale,$ ; [km/s]
      sspfile:      sspfile+'.gz'}

    infofile = outpath+'info_ckc'+ckc_ver+'_'+imfstr+'.fits'
    im_mwrfits, info, infofile, /clobber
    
return
end
