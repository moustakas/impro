;+
; NAME:
;       IM_SPECRES_LAMP()
;
; PURPOSE:
;       Measure the instrumental resolution from a 2D arc lamp image. 
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;      arcfile  - arc lamp FITS file name
;
; OPTIONAL INPUTS:
;      datapath - I/O path
;      nmed_arc - median filter NMED rows centered on the middle row
;                 of the arc lamp spectrum (default 10)
;      nline    - identify and fit NLINE emission lines in the
;                 spectrum (default 50)
;      nback    - number of background polynomial terms to include in
;                 the fit (default 2)
;      npoly    - fit an NPOLY order polynomial to the resolution as a
;                 function of wavelength (default linear)
;      psname   - name of the postscript file (default lampres.ps)
;
; KEYWORD PARAMETERS:
;      postscript - generate postscript output
;      write      - write out the results
;
; OUTPUTS:
;      resinfo  - structure with the results of the fitting
;
; OPTIONAL OUTPUTS:
;
;
; PROCEDURES USED:
;       CWD(), RD2DSPEC(), MAKE_WAVE(), DJS_MEDIAN(), IM_NORMALIZE(), 
;       FIND_NPEAKS(), POLY_ARRAY(), ILINEBACKFIT(), DFPSPLOT,
;       DJS_PLOT, DJS_OPLOT, POLY_ITER()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 February 20, U of A
;       jm04may09uofa - generalized, substantial updates, error
;                       catching, WRITE keyword added
;
; Copyright (C) 2003-2004, John Moustakas
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

function im_specres_lamp, arcfile, datapath=datapath, nmed_arc=nmed_arc, nline=nline, $
  nback=nback, npoly=npoly, psname=psname, postscript=postscript, write=write

    narc = n_elements(arcfile) 
    
    if (narc eq 0L) then begin
       print, 'Syntax - '
       return, -1L
    endif

    light = 2.99792458D5
    fwhm2sig = 2.0*sqrt(2.0*alog(2.0))

    if n_elements(datapath) eq 0L then datapath = cwd()

    if n_elements(nmed_arc) eq 0L then nmed_arc = 10L
    if n_elements(nline) eq 0L then nline = 35L
    if n_elements(nback) eq 0L then nback = 1L
    if n_elements(npoly) eq 0L then npoly = 1L
    if n_elements(psname) eq 0L then psname = 'lampres.ps'

; initialize the output structure    

    resinfo = {$
      nline:         nline, $
      nback:         nback, $
      npoly:         npoly, $
      linewave:      dblarr(nline), $
      linesigma:     dblarr(nline), $
      linewidth:     dblarr(nline), $
      coeff:         dblarr(npoly+1)}
    resinfo = replicate(resinfo,narc)

    plotsym, 8, 2.0;, /fill

    splog, 'Opening '+datapath+psname+'.'
    im_openclose, datapath+psname, postscript=postscript, /color, /silent
    
    for iarc = 0L, narc-1L do begin
    
; read in the arc lamp and extract the center row

       if file_test(datapath+arcfile[iarc],/regular) eq 0L then begin
          splog, 'Arc lamp file '+datapath+arcfile[iarc]+' not found.'
          return, -1L
       endif else splog, 'Reading '+datapath+arcfile[iarc]+'.'
       
       arc = mrdfits(datapath+arcfile[iarc],0,header,/silent)
       arcerr = sqrt(arc)
       arcinvvar = 1.0/arcerr^2.0
       wave = make_wave(header,cd1_1=cd1_1,crval1=crval1)

       arcdims = size(arc,/dimension)
       ncols = arcdims[0]
       nrows = arcdims[1]
       midrow = nrows/2L

       flux = djs_median(arc[*,midrow-nmed_arc/2L:midrow+nmed_arc/2L],2) 
       flux_stats = im_stats(flux)
       flux = im_normalize(flux,/max,const=norm)
       
       ferr = djs_median(arcerr[*,midrow-nmed_arc/2L:midrow+nmed_arc/2L],2)
       ferr = ferr/norm
       invvar = 1.0/ferr^2.0

;; find the NLINE brightest emission lines and sort them by increasing 
;; wavelength
;
;;   linewave = find_npeaks(flux,wave,nfind=nline,minsep=minsep,$
;;     width=10,ypeak=ypeak,xerr=ferr,npeak=npeak)
;;   linewave = linewave[sort(linewave)] 

; trace on the full 2D image to find good lines to fit
       
       xcen = trace_crude(arc,arcinvvar,radius=radius,thresh=3.0*flux_stats.median)
       xcen = trace_fix(xcen,minsep=10.0,ngrow=ngrow,ycen=ycen,xerr=xerr)
       nfind = (size(xcen,/dimension))[1]

       linewave = float(reform(crval1+cd1_1*xcen[midrow,*]))
       lineres = linewave*0.0
       
;      linewave = linewave[fix(randomu(seed,nline)*n_elements(linewave))]
;      linewave = linewave[uniq(linewave,sort(linewave))]
;      nline = n_elements(linewave)
       
;      plot, wave, flux, ps=10, xsty=3, ysty=3;, xrange=[5500,6200]
;      for i = 0L, nline-1 do oplot, [linewave[i],linewave[i]], [0,0.5], line=2

; fit the emission lines with multiple Gaussians;  constrain the
; Gaussians to have the same redshift and to be positive.  initialize
; the (constant) background terms
       
       background = poly_array(ncols,nback)

       zindex = lonarr(nfind)
       findex = lindgen(nfind)
       windex = lindgen(nfind)
       fvalue = replicate(0.5,nfind)
       
       splog, 'Fitting '+strn(nfind)+' arc lamp emission lines.'
       t0 = systime(1)
       result = ilinefit(flux,wave,linewave,lineres,invvar=invvar,$
         background=background,zindex=zindex,findex=findex,windex=windex,$
         fvalue=fvalue,bterms=bterms,specfit=specfit,bfit=bfit,$
         sigmax=1000.0)
       splog, 'Time to fit: '+string((systime(1)-t0),format='(F6.2)')+' seconds.'

;      good = where((result.linearea_err gt 0.0) and (result.linesigma_err gt 0.0),ngood)
;      if ngood ne 0L then begin
;         
;         linewave = linewave[good]
;         linesigma = result[good].linesigma
;         linesigma_err = result[good].linesigma_err
;
;      endif else message, 'No good lines found!'

       linesigma = result.linesigma
       linesigma_err = result.linesigma_err
       linewidth = float(linewave*linesigma/light)
       linewidth_err = float(linewave*linesigma_err/light)

; fit a polynomial to the data with iterative sigma rejection

       coeff = robust_poly_fit(linewave,linewidth,npoly)
;      poly_iter, linewave, linewidth, npoly, 2.0, coeff=coeff
;      coeff = poly_fit(linewave,linewidth,npoly)
       polyfit = poly(wave,coeff)

; generate a plot of the fitting       
       
       djs_plot, wave, flux, xsty=3, ysty=3, ps=10, charsize=2.0, charthick=5.0, $
         xthick=5.0, ythick=5.0, ytitle='Normalized Flux', position=[0.17,0.45,0.97,0.97], $
         xtickname=replicate(' ',10), thick=5.0, xrange=minmax(wave)
       djs_oplot, wave, specfit, ps=10, color='red', thick=5.0
       legend, repstr(repstr(arcfile[iarc],'.fits',''),'.gz',''), /left, $
         /top, box=0, charsize=2.0, charthick=5.0
       
       yrange = [0,max(linewidth)*1.3]
       djs_plot, linewave, linewidth, ps=8, xsty=3, ysty=3, $
;      ploterror, linewave, linewidth, linewidth_err, ps=8, xsty=3, ysty=3, $
         charsize=2.0, charthick=5.0, xrange=minmax(wave), $
         xthick=5.0, ythick=5.0, ytitle=textoidl('\sigma [\AA]'), $
         position=[0.17,0.13,0.97,0.45], xtitle='Wavelength (\AA)', $
         /noerase, yminor=2, yrange=yrange
       djs_oplot, wave, polyfit, line=0, thick=3.0, color='green'
;      legend, string(reform(coeff),format='(G0.3)'), /right, /top, box=0, $
;        charsize=1.5, charthick=5.0

; store the results; only keep the first NLINE wavelengths

       resinfo[iarc].linewave = linewave[0:nline-1L]
       resinfo[iarc].linewidth = poly(resinfo[iarc].linewave,coeff)
       resinfo[iarc].linesigma = light*resinfo[iarc].linewidth/resinfo[iarc].linewave
       resinfo[iarc].coeff = coeff

    endfor

    splog, 'Closing '+datapath+psname+'.'
    im_openclose, postscript=postscript, /close, /silent

    if keyword_set(write) then begin
       resname = repstr(psname,'.ps','.fits')
       splog, 'Writing '+datapath+resname+'.'
       mwrfits, resinfo, datapath+resname, /create
       spawn, ['gzip -f '+datapath+resname], /sh
    endif
    
return, resinfo
end
