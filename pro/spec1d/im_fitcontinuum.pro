;+
; NAME:
;   IM_FITCONTINUUM()
;
; PURPOSE:
;   Fit a spectral continuum using a variety of methods.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;   wave   - wavelength vector [NPIX]
;   flux   - corresponding spectrum [NPIX]
;   ferr   - corresponding error spectrum [NPIX] (not incorporated)
;
; OPTIONAL INPUTS:
;   nfind   - if MASK is set then mask NFIND peaks (default 20) 
;   nsmooth - if MASK is set then smooth the mask by NSMOOTH
;             pixels (default 20) 
;   
;   method - 1: b-spline (default)
;            2: median box-car smoothing
;            3: linear combination of polynomials
;            4: Legendre fit (least reliable)
;
;   method 1:
;      boxstat - MIN, MAX, MEAN, or MEDIAN (default), see
;                DKBOXSTATS(); this parameter should be set to MIN
;                for an emission-line spectrum and MAX for an
;                absorption-line spectrum; a compromise is MEDIAN
;      xwidth  - replace each pixel with the DXBOXSTAT pixel value
;                in an array of width XWIDTH (default 30)
;      bsorder - order of the b-spline fit (default 15)
;
;   method 2:
;      medwidth    - width of the median filter, see DJS_MEDIAN()
;                    (default 150 pixels)
;      smoothwidth - width of the smoothing window, see SMOOTH()
;                    (default 49 pixels)
;
;   method 3:
;      npoly   - number of polynomials of increasing order to fit
;                (default 4)
;
;   method 4:
;      ncoeff  - number of coefficients in the Legendre fit
;                (default 4)
;   plotthick  - plot thickness (default 2.0)
;
;
; KEYWORD PARAMETERS:
;   mask     - mask emission or absorption lines and iteratively
;              and  reject outlying pixel values before fitting 
;   tellmask - only mask telluric features
;   doplot   - generate a plot of the fit and the masked pixels 
;   debug    - wait for a keystroke after making the plot 
;   silent   - suppress messages to STDOUT
;
; OUTPUTS:
;   continuum - the fit to the continuum [NPIX]
;
; OPTIONAL OUTPUTS:
;   nocflux - FLUX divided by the continuum [NPIX]
;   sset    - if METHOD=1 then return the b-spline structure
;
; COMMENTS:
;   if MASK is set then use an iterative procedure.  Fit the
;   continuum without masking to do a coarse normalization, then
;   identify peaks in the spectrum, mask those pixels, and re-fit
;   the continuum.
;
; PROCEDURES USED:
;   FIND_NPEAKS(), EMISSION_MASK(), FIND_NPEAKS(), DJS_ITERSTAT,
;   GET_ELEMENT, DKBOXSTATS(), BSPLINE_ITERFIT(), DJS_MEDIAN(),
;   POLY_ARRAY(), MPFITFUN(), POLY_ARRAY(), MPFITFUN(),
;   FUNC_FIT(), POLYLEG(), DJS_PLOT, DJS_OPLOT
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 November 4, U of A
;   jm03feb20uofa - rewritten and documented
;   jm04jul20uofa - added PLOTTHICK input
;   jm04sep12uofa - added TELLMASK keyword
;   jm05jun30uofa - bug fix when TELLMASK=1 and no telluric
;                   features are within the wavelength range 
;
; Copyright (C) 2002-2005, John Moustakas
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

function continuum_model, x, p ; for METHOD 3
return, x # p
end

function im_fitcontinuum, wave, flux, ferr, method=method, nfind=nfind, $
  nsmooth=nsmooth, boxstat=boxstat, xwidth=xwidth, bsorder=bsorder, $
  medwidth=medwidth, smoothwidth=smoothwidth, npoly=npoly, ncoeff=ncoeff, $
  plotthick=plotthick, nocflux=nocflux, sset=sset, ccoeff=ccoeff, mask=mask, $
  tellmask=tellmask, doplot=doplot, debug=debug, silent=silent, $
  _extra=extra

    npix = n_elements(wave)

    if n_elements(method) eq 0L then method = 1L ; default to the b-spline fit
    if n_elements(plotthick) eq 0L then plotthick = 2.0

; mask out peaks (emission and absorption lines).  in the mask 0B is
; good, 1B is bad 
    
    if keyword_set(mask) then begin 

       cfit = im_fitcontinuum(wave,flux,method=method,boxstat=boxstat,$
         xwidth=xwidth,bsorder=bsorder,medwidth=medwidth,$
         smoothwidth=smoothwidth,npoly=npoly,ncoeff=ncoeff,$
         plotthick=plotthick,mask=0,doplot=0,/silent)

       normflux = flux/cfit
       
       if n_elements(nfind) eq 0L then nfind = 20L
       if n_elements(nsmooth) eq 0L then nsmooth = 20L

       mask = make_array(npix,/byte,value=0) ; bad pixel mask

; find emission and absorption lines as peaks.  smooth the mask by
; NSMOOTH pixels on either side of the peak
       
       epeaks = find_npeaks(normflux,wave,nfind=nfind,minsep=minsep,$
         width=10,ypeak=ypeak,xerr=ferr,npeak=npeak)
       apeaks = find_npeaks(-(normflux+min(-normflux)),wave,nfind=nfind,minsep=minsep,$
         width=10,ypeak=ypeak,xerr=ferr,npeak=npeak)

       get_element, wave, epeaks, exx
       mask[exx] = 1B

       get_element, wave, apeaks, axx
       mask[axx] = 1B

       mask = smooth(float(mask),nsmooth) gt 0B

; mask out known emission lines

;      emask = emission_mask(wave,z=z,width=30)
;      mask = (mask + (emask eq 0B)) gt 0B

       good = where(mask eq 0B,comp=bad,ncomp=nbad)
       
; iteratively reject outliers

       djs_iterstat, normflux[good], sigrej=3.0, maxiter=50, mask=imask
       mask[good] = (mask[good] + (imask eq 0B)) gt 0B

       good = where(mask eq 0B,comp=bad)

;      get_element, wave, 6800.0, ohstart ; mask out night sky emission
;      mask[ohstart<(npix-1L):npix-1L] = 0
       
    endif else begin

       mask = make_array(npix,/byte,value=0)
       good = where(mask eq 0B,ngood)
       nbad = 0L

    endelse
    
; mask out known telluric features    
    
    if keyword_set(tellmask) then begin

       tellmask = telluric_mask(wave,good=tellgood,bad=tellbad)
       mask = (mask + (tellmask eq 0B)) ge 1B
       
       good = where(mask eq 0B,comp=bad,ngood,ncomp=nbad)

;      if (good[0] eq -1L) then ngood = 0L else ngood = n_elements(good)
;      if (bad[0] eq -1L) then nbad = 0L else nbad = n_elements(bad)
       
    endif

    case method of

; b-spline fit to a "minimum" smoothed spectrum
; --------------------------------------------------
       1L: begin

          if not keyword_set(silent) then splog, 'METHOD 1: b-spline fit'

          if n_elements(xwidth) eq 0L then xwidth = 30
          if n_elements(boxstat) eq 0L then boxstat = 'median'

          if n_elements(bsorder) eq 0L then bsorder = 10.0
          bkptres = (max(wave)-min(wave))/bsorder ; break point resolution [Angstrom]
          bkpt = findgen(bsorder)*bkptres+min(wave)

          cont = dkboxstats(flux[good],xwidth=xwidth,boxstat=boxstat) ; continuum
          sset = bspline_iterfit(wave[good],flux[good],bkpt=bkpt,yfit=yfit,$
            nord=4,lower=0.1,upper=0.1,/silent)
          continuum = bspline_valu(wave,sset)
          
       end

; median box-car smoothing of the spectrum
; --------------------------------------------------
       2L: begin

          if not keyword_set(silent) then splog, 'METHOD 2: median box-car smoothing'

          if n_elements(medwidth) eq 0L then medwidth = 150
          if n_elements(smoothwidth) eq 0L then smoothwidth = 49 < (ngood-1L)
          temp = smooth(djs_median(flux[good],width=medwidth),smoothwidth,/edge_truncate)
          continuum = interpol(temp,wave[good],wave) ; interpolate over the masked points
; continuum = interpol(temp,wave[good],wave,/quadratic)

       end

; MPFIT method: find the best linear combination of polynomials
; --------------------------------------------------
       3L: begin
          
          if not keyword_set(silent) then splog, 'METHOD 3: polynomial linear combination'

          if n_elements(npoly) eq 0L then npoly = 4
          polyflux = poly_array(npix,npoly)
          
          parinfo = {value: 0.0D, fixed: 0L, limited: [0,0], limits: [0.0D,0.0D]}
          parinfo = replicate(parinfo,npoly)

          weights = mask*0.0 + 1.0
          if nbad ne 0L then weights[bad] = 1D-10

          polycoeff = mpfitfun('continuum_model',polyflux,flux,weights=weights,/quiet,$
            parinfo=parinfo,niter=niter,covar=covar,perror=perror,bestnorm=bestnorm)

          continuum = polyflux # polycoeff

       end

; legendre fit with no smoothing; fixed continuum order
; -----------------------------------------------------
       4L: begin

          if not keyword_set(silent) then splog, 'METHOD 4: Legendre polynomial fit'

          if n_elements(ncoeff) eq 0L then ncoeff = 4L

          invvar = flux[good]*0.0+1D20
          ccoeff = func_fit(wave[good],flux[good],ncoeff,invvar=invvar*(mask[good] eq 0B),$
            function_name='flegendre')
          continuum = polyleg(wave,ccoeff)
          
       end

    endcase

; generate a plot
    
    if keyword_set(doplot) then begin

       djs_plot, wave, flux, ps=10, xsty=3, ysty=3, charsize=1.2, charthick=plotthick, $
         xthick=plotthick, ythick=plotthick, ytitle='Flux', position=[0.13,0.45,0.97,0.97], $
         xtickname=replicate(' ',10), _extra=extra
       if nbad ne 0L then djs_oplot, wave[bad], flux[bad], color='green', ps=4
       djs_oplot, wave, continuum, color='red', thick=plotthick

       djs_plot, wave, flux/continuum, ps=10, xsty=3, ysty=3, charsize=1.2, charthick=plotthick, $
         xthick=plotthick, ythick=plotthick, ytitle='Normalized Flux', position=[0.13,0.13,0.97,0.45], $
         xtitle='Wavelength (\AA)', /noerase, yminor=2
       if nbad ne 0L then djs_oplot, wave[bad], flux[bad]/continuum[bad], color='green', ps=4
       djs_oplot, !x.crange, [1,1], color='red', thick=plotthick

       if n_elements(debug) ne 0L then begin
          if not keyword_set(silent) then splog, 'Press any key to continue.'
          cc = get_kbrd(1)
       endif
       
    endif

    nocflux = flux/continuum
    
return, continuum
end
