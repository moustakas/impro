;+
; NAME:
;       IM_SPLOTEW()
;
; PURPOSE:
;       Measure the equivalent width of emission or absorption lines
;       in the spirit of IRAF's SPLOT task (fit the continuum locally,
;       but estimate the error in the continuum using Monte Carlo).  
;
; CALLING SEQUENCE:
;       result = im_splotew(wave,flux,ivar,linewave,linename=,$
;          nmonte=, boxwidth=, /absorption, /doplot, /silent)
;
; INPUTS:
;       wave     - spectrum wavelength in Angstroms [NPIX]
;       flux     - data spectrum in erg/s/cm2/A [NPIX]
;       ivar     - inverse variance spectrum in erg/s/cm2/A [NPIX]
;       linewave - emission- or absorption-line wavelength [NLINE] 
;
; OPTIONAL INPUTS:
;       linename - optional name for each line in LINEWAVE [NLINE]  
;       nmonte   - number of Monte Carlo iterations (default 50)
;       boxwidth - derive the EW within LINEWAVE+/-BOXWIDTH Angstroms
;                  (default 35)
;
; KEYWORD PARAMETERS:
;       absorption - set if fitting absorption lines
;       doplot     - generate a diagnostic plot of the fit and wait
;                    for a keystroke
;       silent     - suppress messages to STDOUT 
;
; OUTPUTS:
;       result - output data structure (see LINEEW_INIT for the
;                structure fields) 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       LINEEW_INIT(), MPFITPEAK(), DJS_PLOT, DJS_OPLOT
;
; COMMENTS:
;       This routine assumes that the dispersion in Angstroms per
;       pixel is constant. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 Jul 09, U of A - written
;       jm05jan06uofa - documented, updated, and streamlined 
;
; Copyright (C) 2002, 2005, John Moustakas
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

function lineew_init, linewave, linename
; function for IM_SPLOTEW()
    lineew = {$
      linename:        linename, $
      linewave:        linewave, $
      linecenter:      0.0,$
      linecenter_flux: 0.0,$
      linewidth:       0.0,$
      continuum:       [0.0,-1.0], $
      gaussflux:       [0.0,-1.0], $
      ew:              [0.0,-1.0]}
return, lineew
end

function im_splotew, wave, flux, ivar, linewave, $
  linename=linename, nmonte=nmonte, boxwidth=boxwidth, $
  absorption=absorption, doplot=doplot, silent=silent

    npix = n_elements(wave)
    nline = n_elements(linewave)

    if (nline eq 0L) or (npix eq 0L) then begin
       doc_library, 'im_splotew'
       return, -1L
    endif

    if (n_elements(flux) ne npix) then begin
       splog, 'WAVE and FLUX are not the same size'
       return, -1L
    endif
    
    if (n_elements(ivar) ne npix) then begin
       splog, 'WAVE and IVAR are not the same size'
       return, -1L
    endif

    if (n_elements(linename) eq 0L) then linename = string(linewave,format='(G0)')

; call this routine recursively    
    
    if (nline gt 1L) then begin

       if (n_elements(boxwidth) eq 1L) then boxwidth = replicate(boxwidth,nline)
       if (n_elements(boxwidth) eq 0L) then boxwidth = replicate(35.0,nline)
       
       for i = 0L, nline-1L do begin
          lineew1 = im_splotew(wave,flux,ivar,linewave[i],$
            linename=linename[i],nmonte=nmonte,$
            boxwidth=boxwidth[i],doplot=doplot,silent=silent)
          if (i eq 0L) then lineew = lineew1 else $
            lineew = [lineew,lineew1]
       endfor
       return, lineew
    endif

    if (not keyword_set(silent)) then splog, 'Fitting '+linename
    
; default settings
    if (n_elements(nmonte) eq 0L) then nmonte = 50L      ; number of Monte-Carlo realizations
    if (n_elements(boxwidth) eq 0L) then boxwidth = 35.0 ; Angstrom (should be about 20 Angstrom for H-gamma)

    nterms = 5             ; Gaussian + line
    cd1_1 = wave[1]-wave[0] ; Angstrom/pixel (assumed constant)

    lineew = lineew_init(linewave,linename) ; initialize the output structure
    
; check for an out-of-range line
    if (linewave gt max(wave)) or (linewave lt min(wave)) then begin
       splog, 'Line '+linename+' is beyond the specified wavelength range.'
       return, lineew
    endif
       
; zoom in on the line
    zoomline = where((wave gt linewave-boxwidth) and $
      (wave lt linewave+boxwidth),nwave)
    if (nwave eq 0L) then begin
       splog, 'A problem was encountered zooming into the line.'
       return, lineew
    endif

    lwave = wave[zoomline]
    lflux = flux[zoomline]
    livar = ivar[zoomline]

; constrain the central wavelength of the Gaussian
    parinfo = {value: 1.0D, limited: [0,0], limits: [0.0D,0.0D]}
    parinfo = replicate(parinfo,nterms)
    parinfo[1].value = linewave
    parinfo[1].limited = 1
    parinfo[1].limits = linewave + 5.0*[-1,1] ; [-10.0,+10.0] ; +/- 5 Angstrom

    indx = where((lwave gt parinfo[1].limits[0]) and $
      (lwave lt parinfo[1].limits[1]))
    if keyword_set(absorption) then $
      parinfo[0].value = min(lflux[indx]) else $ ; initial guess
      parinfo[0].value = max(lflux[indx])

; first iteration: fit the line to get a rough estimate of the
; line-centroid and sigma-width
    yfit = mpfitpeak(lwave,lflux,fitcoeff,nterms=nterms,$
      weight=livar,/gaussian,negative=absorption,$
      perror=perror,quiet=quiet,bestnorm=bestnorm,$
      parinfo=parinfo,status=status)

; second iteration: zoom to within +/-5*sigma and then refit
    nsigma = 8.0
    indx = where((wave gt fitcoeff[1]-nsigma*fitcoeff[2]) and $
      (wave lt fitcoeff[1]+nsigma*fitcoeff[2]),nindx)
    if (nindx eq 0) then message, 'Fix me'
    
    lwave = wave[indx]
    lflux = flux[indx]
    livar = ivar[indx]

    parinfo.value = fitcoeff
    parinfo[1].limited = 0 ; lift the restriction on the line-center
    yfit = mpfitpeak(lwave,lflux,fitcoeff,nterms=nterms,$
      weight=livar,/gaussian,negative=absorption,$
      perror=perror,quiet=quiet,bestnorm=bestnorm,$
      parinfo=parinfo,status=status)

    gaussflux = sqrt(2.0*!dpi)*fitcoeff[2]*fitcoeff[0] ; gaussian flux
    linecenter = fitcoeff[1]
    linewidth = fitcoeff[2]
    linecenter_flux = total(lwave*lflux)/total(lflux) ; flux-weighted center
    
; integrate the continuum model to get the EW    
;   lowave = linecenter-3.0*linewidth
;   hiwave = linecenter+3.0*linewidth

    clineflux = poly(lwave,fitcoeff[3:4])
    continuum = interpol(clineflux,lwave,linecenter)
    ew = gaussflux/continuum
    
;   lineflux = lflux/clineflux
;   ivar_lineflux = livar*clineflux^2.0
;   ew = im_integral(lwave,ivar_lineflux*lineflux,lowave,hiwave)/$
;     im_integral(lwave,ivar_lineflux,lowave,hiwave)

; pack it in
    lineew.linecenter = linecenter
    lineew.linecenter_flux = linecenter_flux
    lineew.linewidth = linewidth
    lineew.continuum[0] = continuum
    lineew.gaussflux[0] = gaussflux
    lineew.ew[0] = ew

; debugging plot    
    if keyword_set(doplot) then begin
       power = ceil(abs(alog10(median(flux))))
       scale = 10.0^power
       ytitle1 = 'Flux (10^{-'+string(power,format='(I0)')+'} '+flam_units()+')'
       xrange = linewave+[-1,1]*boxwidth>(linewave-min(lwave))>(max(lwave)-linewave)
       xrange = [(linewave-boxwidth)<min(lwave),(linewave+boxwidth)>max(lwave)]
       djs_plot, wave, scale*flux, xsty=3, ysty=3, ps=10, charsize=1.8, $
         xtitle='Wavelength (\AA)', ytitle=ytitle1, $
         xrange=xrange
       djs_oplot, linewave*[1,1], !y.crange, thick=3, line=5
       djs_oplot, lwave, scale*yfit, ps=10, thick=4.0, color='red'
       djs_oplot, lwave, scale*clineflux, thick=2.0, color='blue'
;      djs_oplot, [lwave[indx[0]],lwave[indx[0]]], !y.crange, line=2, thick=2.0
;      djs_oplot, [lwave[indx[nindx-1L]],lwave[indx[nindx-1L]]], !y.crange, line=2, thick=2.0
    endif

; Monte Carlo to get the errors
    lerr = livar*0.0
    notzero = where(livar gt 0,nnotzero)
    if (nnotzero eq 0) then message, 'Fix me'
    lerr[notzero] = 1.0/sqrt(livar[notzero])

    if (nmonte gt 0) then begin
       monte_gaussflux = fltarr(nmonte)
       monte_continuum = fltarr(nmonte)
       
       for jj = 0L, nmonte-1L do begin
          newflux = lflux+randomn(seed,nindx)*lerr
          yfit = mpfitpeak(lwave,newflux,fitcoeff,nterms=nterms,$
            weight=livar,/gaussian,negative=absorption,$
            perror=perror,quiet=quiet,bestnorm=bestnorm,$
            parinfo=parinfo,status=status)
          
          monte_gaussflux[jj] = sqrt(2.0*!dpi)*fitcoeff[2]*fitcoeff[0] ; gaussian flux

          clineflux = poly(lwave,fitcoeff[3:4])
          monte_continuum[jj] = interpol(clineflux,lwave,linecenter)

;         lineflux = lflux/clineflux
;         ivar_lineflux = livar*clineflux^2.0
;         monte_ew[jj] = im_integral(lwave,ivar_lineflux*lineflux,lowave,hiwave)/$
;           im_integral(lwave,ivar_lineflux,lowave,hiwave)
       endfor

       monte_ew = monte_gaussflux/monte_continuum
       
; pack it in    
       lineew.gaussflux[1] = djsig(monte_gaussflux)
       lineew.continuum[1] = djsig(monte_continuum)
       lineew.ew[1] = djsig(monte_ew)
    endif
       
return, lineew
end
