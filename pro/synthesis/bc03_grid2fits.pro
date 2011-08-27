;+
; NAME:
;       BC03_GRID2FITS
;
; PURPOSE:
;       Convert BC03 data structures to individual FITS files with the
;       option of writing error spectra.
;
; CALLING SEQUENCE:
;       bc03_grid2fits, grid, fitsgrid, mgalaxy=, redshift=, $
;          fwhmres=, starvdisp=, minwave=, maxwave=, specsnr=, snrwave=, $
;          hbeta_ew=, ebv_gas=, dustratio=, outfile=, outpath=, /espectrum, $
;          /wfits, _extra=extra
;
; INPUTS:
;       grid - input data structure from IM_READ_BC03()
;
; OPTIONAL INPUTS:
;       mgalaxy   - galaxy mass (default 1 M_sun) [M_sun]
;       redshift  - galaxy redshift (default 0.0)
;       fwhmres   - *observed* FWHM spectral resolution; must be
;                   greater than or equal to 3*(1+redshift) [Angstrom] 
;       starvdisp - stellar velocity dispersion (default 120) [km/s] 
;       minwave   - *rest-frame* minimum wavelength [Angstrom]
;       maxwave   - *rest-frame* maximum wavelength [Angstrom]
;       specsnr   - *observed* signal-to-noise per resolution element; 
;                   if provided then construct an error spectrum; in
;                   addition, SNRWAVE must be specified 
;       snrwave   - *observed* wavelength corresponding to SPECSNR
;                   [Angstrom]
;       hbeta_ew  - if ESPECTRUM=1 then scale the emission-line
;                   spectrum to this EW(H-beta) (default 10.0); NB:
;                   H-beta must be within the interval
;                   [MINWAVE,MAXWAVE] +/- 50 Angstroms
;       ebv_gas   - redden the emission-line spectrum (if ESPECTRUM=1)
;                   by this amount [mag]
;       dustratio - E(B-V)_stars/E(B-V)_gas; must be in the interval
;                   [0,1], inclusive (default 1.0)
;       outfile   - output file name (default
;                   'bc03_'+AGESUFFIX+'.fits') 
;       outpath   - output path name (default CWD)
;       extra     - inputs for IM_CONSTRUCT_ESPECTRUM() and K_LAMBDA()
;
; KEYWORD PARAMETERS:
;       espectrum - superpose a model emission-line spectrum (keywords
;                   must be passed via _EXTRA)
;
; OUTPUTS:
;       fitsgrid - output informational data structure
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       DLUMINOSITY(), IM_GAUSS_BROADEN(), GCONVOLVE(), MKHDR,
;       SXDELPAR, SXADDPAR, IM_CONSTRUCT_ESPECTRUM(), MWRFITS, SPLOG,
;       SPECTRAL_INDICES(), STRUCT_ADDTAGS(), K_LAMBDA() 
;
; COMMENTS:
;    The wavelength vector is re-gridded to a uniform 1 Angstrom per
;    pixel.  
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Mar 01 - written
;       jm04oct26uofa - totally re-written and expanded 
;       jm05apr13uofa - various bug fixes; added DUSTRATIO optional
;                       input 
;       jm05may05uofa - also return the emission-line and continuum
;                       spectra 
;
; Copyright (C) 2004-2005, John Moustakas
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

pro bc03_grid2fits, grid, fitsgrid, mgalaxy=mgalaxy, redshift=redshift, $
  fwhmres=fwhmres, starvdisp=starvdisp, minwave=minwave, maxwave=maxwave, $
  specsnr=specsnr, snrwave=snrwave, hbeta_ew=hbeta_ew, ebv_gas=ebv_gas, $
  dustratio=dustratio, outfile=outfile, outpath=outpath, espectrum=espectrum, $
  wfits=wfits, _extra=extra

    lsun = 3.826D33             ; [erg/s]
    light = 2.99792458D5        ; speed of light [km/s]
    sfr2energy = 1.0/7.9D-42    ; (erg/s) / (M_sun/yr) [Kennicutt 1998]

    hahb = 2.86
;   hahb = return_tbalmer(/HaHb)
    fwhm2sig = 2.0*sqrt(2.0*alog(2.0))

    ngrid = n_elements(grid)
    if (ngrid eq 0L) then begin
       print, 'Syntax - bc03_grid2fits, grid, fitsgrid, mgalaxy=, redshift=, $'
       print, '   fwhmres=, starvdisp=, minwave=, maxwave=, specsnr=, snrwave=, $'
       print, '   hbeta_ew=, ebv_gas=, dustratio=, outfile=, outpath=, /espectrum, $'
       print, '   /wfits, _extra=extra'
       return
    endif

    if (n_elements(mgalaxy) eq 0L) then mgalaxy = 1.0   ; [M_sun]

    if (n_elements(redshift) eq 0L) then begin
       redshift = 0.0D 
       fluxarea = 1.0
    endif else begin
       dlum = dluminosity(redshift,/cm)
       fluxarea = 4.0*!dpi*dlum*dlum
    endelse

    minfwhmres = 3.0            ; minimum FWHM spectral resolution [Angstrom]

    if (n_elements(fwhmres) eq 0L) then fwhmres = minfwhmres*(1+redshift)
    if (fwhmres lt minfwhmres*(1+redshift)) then begin
       splog, 'FWHMRES (observed) must be larger than 3*(1+REDSHIFT).'
       return
    endif

    if (n_elements(starvdisp) eq 0L) then starvdisp = 120.0 ; [km/s]

    restfwhmres = fwhmres/(1.0+redshift) ; rest-frame FWHM spectral resolution [Angstrom]
    restdwave = 1.0                      ; rest-frame dispersion [Angstrom/pixel]

    coeff1 = 2D-4                        ; desired pixel size, ~60 km/s [log-Angstrom]
    sigma_vdisp = starvdisp/coeff1/light ; stellar velocity dispersion [pixel]

    wave = grid.wave
    
    if (n_elements(minwave) eq 0L) then minwave = min(wave) else minwave = minwave>min(wave)
    if (n_elements(maxwave) eq 0L) then maxwave = max(wave) else maxwave = maxwave<max(wave)

    reddwave = restdwave*(1+redshift) ; observed dispersion [Angstrom/pixel]

    oldwave = wave
    wave = findgen((maxwave-minwave)/restdwave+1)*restdwave+minwave
    redwave = wave*(1+redshift)
    
    npix = n_elements(wave)
    
    nage = n_elements(grid.age)
    nspecsnr = n_elements(specsnr)
    nsnrwave = n_elements(snrwave)

    if (nspecsnr eq 1L) then begin
       errormake = 1L
       if (nsnrwave ne 1L) then begin
          splog, 'SPECSNR and SNRWAVE must be provided.'
          return
       endif
       if (snrwave le minwave*(1+redshift)) or (snrwave ge maxwave*(1+redshift)) then begin
          splog, 'SNRWAVE must be in the interval [MINWAVE,MAXWAVE].'
          return
       endif
    endif else errormake = 0L

    if (n_elements(hbeta_ew) eq 0L) then hbeta_ew = 0.0

    if keyword_set(espectrum) then if (4861.0-50.0 lt minwave) or $
      (4861.0+50.0 gt maxwave) then begin
       splog, 'H-beta not in the interval [MINWAVE,MAXWAVE]: Setting ESPECTRUM=0.'
       espectrum = 0L
    endif
    
    if (n_elements(ebv_gas) eq 0L) then ebv_gas = 0.0
    if (n_elements(dustratio) eq 0L) then dustratio = 1.0
    if (dustratio lt 0.0) or (dustratio gt 1.0) then begin
       splog, 'DUSTRATIO must be between [0,1], inclusive.'
       return
    endif

    ebv_stars = dustratio*ebv_gas
    
    if (n_elements(outpath) eq 0L) then outpath = cwd()

    fitsgrid = {$
      file:            '', $
      agesuffix:       '', $
      age:             0.0, $
      mgalaxy:         mgalaxy, $
      redshift:        redshift, $
      fluxarea:        fluxarea, $
      fwhmres:         fwhmres,  $ ; observed
      starvdisp:       starvdisp,$
      minwave:         minwave,  $ ; rest frame
      maxwave:         maxwave,  $ ; rest frame
      hbeta_ew:        hbeta_ew, $ ; rest frame
      hbeta_continuum: 0.0,      $ ; rest frame
      ebv_gas:         ebv_gas,  $
      ebv_stars:       ebv_stars,$
      sfr:             0.0,      $
      babs_halpha_ew:  0.0,      $ ; rest frame
      babs_hbeta_ew:   0.0,      $ ; rest frame
      babs_hgamma_ew:  0.0,      $ ; rest frame
      babs_hdelta_ew:  0.0,      $ ; rest frame
      wave:            fltarr(npix), $ ; observed
      flux:            fltarr(npix), $ ; observed
      ferr:            fltarr(npix), $ ; observed
      espectrum:       fltarr(npix), $ ; observed
      continuum:       fltarr(npix), $ ; observed
      header:          strarr(20)}
    if errormake then fitsgrid = struct_addtags(fitsgrid,{$
      specsnr: specsnr, $
      snrwave: snrwave})        ; observed frame

    fitsgrid = replicate(fitsgrid,nage)

    for iage = 0L, nage-1L do begin

       age = grid.age[iage]/1E9 ; [Gyr]
       fitsgrid.age = age
       
       oldflux = reform(grid.flux[*,iage])
       flux = interpol(oldflux,oldwave,wave) ; do not use this vector

; redshift the flux vector       

       constant = lsun*mgalaxy/fluxarea/(1+redshift)
       redflux = constant*interpol(oldflux,oldwave,wave) ; [erg/s/cm2/A]
       
; broaden to the *observed* spectral resolution, where the 3 Angstroms
; limiting spectral resolution has already been reduced by the
; redshifting 

       broadflux = im_gauss_broaden(redwave,redflux,minfwhmres*(1+redshift),fwhmres)

; broaden to the stellar velocity dispersion; resample the spectrum
; such that each pixel has a constant width in km/s equal to ~60 km/s;
; finally re-sample back onto the constant dispersion wavelength
; vector; it shouldn't matter whether we do the stellar broadening
; before or after redshifting, since we re-bin to constant pixels in
; km/s 

       logredwave = alog10(redwave)
       coeff0 = logredwave[0]
       lognpix = ceil(1D0 + (max(logredwave)-coeff0)/coeff1)
       logvconstwave = coeff0 + coeff1*findgen(lognpix)
       vconstwave = 10^logvconstwave

       fluxvconst = interpol(broadflux,redwave,vconstwave)
       vbroadflux = gconvolve(fluxvconst,sigma_vdisp,/edge_truncate)
       finalflux = interpol(vbroadflux,vconstwave,redwave)

; redden the continuum
       
       finalflux = finalflux*10^(-0.4*ebv_stars*k_lambda(wave,/charlot))
       fitsgrid[iage].continuum = finalflux

; superpose an emission-line spectrum with the specified EW(H-beta);
; compute the continuum at H-beta to derive the SFR

       if keyword_set(espectrum) then begin

; in this first call to IBALMERABS() all we care about is the
; continuum around H-beta, which is insensitive to BALMERSIGMA 
          
          babs = ibalmerabs(redwave,finalflux,z=redshift,balmerwaves=4861.0,$
            ferr=finalflux*0.1,balmersigma=balmersigma,/silent)
          hbeta_continuum = babs[2].babs_continuum ; rest-frame continuum flux

          sfr = hahb*fluxarea*hbeta_ew*hbeta_continuum/sfr2energy ; [M_sun/yr]
          fitsgrid[iage].sfr = sfr
          fitsgrid[iage].hbeta_continuum = hbeta_continuum

          espec = im_construct_espectrum(redwave,redshift=redshift,linefit=linefit,$
            lineres=fwhmres,sfr=sfr,linewidth=linewidth,ebv=ebv_gas,_extra=extra);,/debug)
          fitsgrid[iage].espectrum = espec

;         niceprint, linefit.linename, linefit.linesigma, $
;           fwhm2sig*linefit.linewave*linefit.linesigma_instr/light, $
;           fwhm2sig*linefit.linewave*linefit.linesigma_total/light
;         niceprint, linefit.linename, linefit.linearea, linefit.linebox
;         print, sfr, linefit[7].linearea*fluxarea/sfr2energy

          trueflux = finalflux
          finalflux = finalflux + espec

; now that the emission-line spectrum has been specified (ie, the
; LINEWIDTH), measure the absorption-line EW's of all the Balmer
; absorption lines; define the absorption window in the *observed*
; frame so that it becomes a redshift-independent quantity

          balmersigma = sqrt((babs.babs_wave*linewidth/light)^2 + (fwhmres/fwhm2sig)^2)
          
          babs = ibalmerabs(redwave,trueflux,z=redshift,ferr=finalflux*0.1,$
            balmerwaves=balmerwaves,balmersigma=balmersigma,/silent);,/debug)

;         struct_print, babs
;         niceprint, babs.babs_line, babs.babs_ew*babs.babs_continuum

          fitsgrid[iage].babs_hdelta_ew = babs[0].babs_ew
          fitsgrid[iage].babs_hgamma_ew = babs[1].babs_ew
          fitsgrid[iage].babs_hbeta_ew = babs[2].babs_ew
          fitsgrid[iage].babs_halpha_ew = babs[3].babs_ew

; measure spectral indices

          indices1 = spectral_indices(wave,trueflux*(1+redshift),$
            ferr=trueflux*0.1,debug=debug,/silent)
          if (iage eq 0L) then indices = indices1 else indices = struct_append(indices,indices1)
          
       endif
       
;      djs_plot, redwave, finalflux, ps=10, xsty=3, ysty=3, thick=1.0, $
;        xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, $
;        xtitle='Wavelength', ytitle='Flux'
;      djs_oplot, redwave, broadflux, ps=10, color='green'
;      djs_oplot, vconstwave, fluxvconst, ps=10, color='cyan'
;      djs_oplot, vconstwave, vbroadflux, ps=10, color='yellow'
;      djs_oplot, redwave, finalflux, ps=10, color='white', thick=2

; construct the error spectrum; assume each resolution element is
; equal to two pixels centered on SNRWAVE
          
       if errormake then begin

          get_element, redwave, snrwave+reddwave*[0,+1], normindx
          finalferr = finalflux/sqrt(finalflux/djs_mean(finalflux[normindx]))/specsnr

          sigmaflux = randomn(seed,npix)*finalferr
          finalflux = finalflux + sigmaflux

       endif else finalferr = finalflux*0.0

       fitsgrid[iage].wave = redwave
       fitsgrid[iage].flux = finalflux
       fitsgrid[iage].ferr = finalferr
       
       mkhdr, header, finalflux, extend=errormake
       sxdelpar, header, 'COMMENT'
       sxaddpar, header, 'OBJECT', 'Object'
       sxaddpar, header, 'CRVAL1', float(redwave[0]), ' central wavelength of first pixel'
       sxaddpar, header, 'CD1_1', float(reddwave), ' dispersion [Angstrom/pixel]'
       sxaddpar, header, 'CRPIX1', 1, ' starting pixel (1-indexed)'
       sxaddpar, header, 'CTYPE1', 'LINEAR'
       sxaddpar, header, 'DC-FLAG', 0, ' log-linear flag'
       sxaddpar, header, 'Z', float(redshift), ' redshift'
       sxaddpar, header, 'VDISP', float(starvdisp), ' stellar velocity dispersion [km/s]'
       sxaddpar, header, 'FWHMRES', float(fwhmres), ' FWHM spectral resolution [Angstrom]'
       sxaddpar, header, 'EBVGAS', float(ebv_gas), ' nebular E(B-V) [mag]'
       sxaddpar, header, 'EBVSTARS', float(ebv_stars), ' stellar E(B-V) [mag]'
       sxaddpar, header, 'EWHBETA', float(hbeta_ew), ' H-beta EW [Angstrom]'

       agesuffix = strtrim(string(age,format='(F12.3)'),2)+'Gyr'
       if (n_elements(outfile) eq 0L) then outfile = 'bc03_'+agesuffix+'.fits'
       
       fitsgrid[iage].file = outfile+'.gz'
       fitsgrid[iage].agesuffix = agesuffix
       fitsgrid[iage].header = header
       
       if keyword_set(wfits) then begin

          splog, 'Writing '+outpath+outfile+'.'
          mwrfits, float(finalflux), outpath+outfile, header, /create
          if errormake then mwrfits, float(finalferr), outpath+outfile
          spawn, ['gzip -f '+outpath+outfile], /sh

       endif

    endfor

; append the Lick index measurements    
    
    fitsgrid = struct_addtags(fitsgrid,indices)
    
return
end
    
