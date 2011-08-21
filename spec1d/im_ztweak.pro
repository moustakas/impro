;+
; NAME:
;       IM_ZTWEAK()
;
; PURPOSE:
;       Compute a small radial velocity shift via cross-correlation
;       between a spectrum and a reference spectrum.
;
; CALLING SEQUENCE:
;       zans = im_ztweak(spec,wave,refspec,refwave,specivar=,mask=, $
;          vmin=,vmax=,nsamp=,minwave=,maxwave=,snrmin=,/doplot,/silent)
;
; INPUTS:
;       spec      - input spectrum [NPIX]
;       wave      - corresponding wavelength vector [NPIX]
;       refspec   - reference spectrum [NREFPIX]
;       refwave   - corresponding wavelength vector [NREFPIX]
;
; OPTIONAL INPUTS:
;       specivar  - inverse variance for SPEC [NPIX]
;       mask      - pixel mask (0 = ignore, 1 = keep) [NPIX] (not
;                   supported) 
;       vmin      - minimum velocity shift allowed (default -400 km/s) 
;       vmax      - maximum velocity shift allowed (default +400 km/s) 
;       nsamp     - oversample the data by NSAMP before
;                   cross-correlating (default 10) 
;       minwave   - only cross-correlate the data between MINWAVE and
;                   MAXWAVE (default WAVE[0])
;       maxwave   - see MINWAVE (default WAVE[NPIX-1]) 
;       snrmin    - minimum S/N in the data spectrum (not supported) 
;
; KEYWORD PARAMETERS:
;       doplot    - generate a plot of the chi2 array
;       silent    - do not print the results to STDOUT
;
; OUTPUTS:
;       zans      - output data structure
;          pshift     - pixel shift
;          pshift_err - pixel shift error
;          zshift     - redshift
;          zshift_err - redshift error
;          vshift     - velocity shift
;          vshift_err - velocity shift error
;          errflag    - error flag 
;             0 - success
;             1 - S/N < SNRMIN
;             2 - no minimum found (edge chi2 solution)
;             3 - PSHIFT/PSHIFT_ERR < 3
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: 
;       Talk about SPECIVAR.
;
; PROCEDURES USED:
;       COMBINE1FIBER, FINDCHI2MIN, DJS_PLOT, DJS_OPLOT, IM_BINSPEC(), 
;       DJS_MEDIAN(), GET_ELEMENT, LEGEND, DJS_ICOLOR(), TEXTOIDL(),
;       SPLOG 
;
; EXAMPLE:
;       IDL> zans = im_ztweak(spec,wave,refspec,refwave,vmin=-100.0,/doplot)
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 13, U of A - written
;       jm03july2uofa - added MASK and SNRMIN keywords
;       jm04mar09uofa - use ZFIND() instead of my crummy
;                       cross-correlation code; major cleanup and
;                       better diagnostic plots 
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

function im_ztweak, spec, wave, refspec, refwave, specivar=specivar, $
  mask=mask, vmin=vmin, vmax=vmax, nsamp=nsamp, minwave=minwave, $
  maxwave=maxwave, snrmin=snrmin, doplot=doplot, silent=silent

    light = 2.99792458D5        ; speed of light [km/s]

    npix = n_elements(spec)
    nwave = n_elements(wave)

    nrefpix = n_elements(refspec)
    nrefwave = n_elements(refwave)

    nspecivar = n_elements(specivar)
    nmask = n_elements(mask)
    
    if (npix eq 0L) or (nwave eq 0L) or $
      (nrefpix eq 0L) or (nrefwave eq 0L) then begin
       print, 'Syntax - zans = im_ztweak()'
       return, -1L
    endif

    if (npix ne nwave) then begin
       print, 'SPEC and WAVE must have the same number of elements.'
       return, -1L       
    endif

    if (nrefpix ne nrefwave) then begin
       print, 'REFSPEC and REFWAVE must have the same number of elements.'
       return, -1L       
    endif

    if (nspecivar eq 0L) then specivar = spec*0.0+1.0 else begin
       if nspecivar ne npix then begin
          print, 'SPECIVAR and SPEC must have the same number of elements.'
          return, -1L       
       endif
    endelse
    
    if nmask eq 0L then mask = fix(spec*0.0)+1 else begin
       if nmask ne npix then begin
          print, 'MASK and SPEC must have the same number of elements.'
          return, -1L       
       endif
    endelse
    
    if n_elements(vmin) eq 0L then vmin = -500.0 ; [km/s]
    if n_elements(vmax) eq 0L then vmax = +500.0 ; [km/s]
    if (not keyword_set(silent)) then begin
       splog, 'Minimum velocity shift: '+string(vmin,format='(G0.0)')+' km/s'
       splog, 'Maximum velocity shift: '+string(vmax,format='(G0.0)')+' km/s'
    endif
    
    zmin = vmin/light
    zmax = vmax/light

    if vmin gt vmax then begin
       print, 'VMIN must be less than VMAX.'
       return, -1L
    endif
    
    if n_elements(nsamp) eq 0L then nsamp = 10.0
    nsamp = 1.0

; constrain the minimum and maximum wavelengths but do not use any
; pixels 10% from each edge    

    if n_elements(minwave) eq 0L then minwave = min(refwave) else minwave = minwave > min(refwave)
    if n_elements(maxwave) eq 0L then maxwave = max(refwave) else maxwave = maxwave < max(refwave)
    if (not keyword_set(silent)) then begin
       splog, 'Minimum wavelength: '+string(minwave,format='(G0.0)')+' A'
       splog, 'Maximum wavelength: '+string(maxwave,format='(G0.0)')+' A'
    endif

    if n_elements(snrmin) eq 0L then snrmin = 10.0

; initialize the output data structure    
    
    zans = {$
      zshift: 0.0, zshift_err: -1.0, $
      vshift: 0.0, vshift_err: -1.0, $
      errflag: 0L} ; 0 = good
    
; resample the data such that they are constant in log-lambda (each
; pixel is the same size in km/s); also oversample by a factor of
; NSAMP

    coeff0 = alog10(minwave) ; starting logarithmic wavelength 
    coeff1 = 1D-4            ; logarithmic pixel size [log-Angstrom]

    pmin = fix(alog10(1+vmin/light)*nsamp/coeff1) ; minimum oversampled pixel shift
    pmax = fix(alog10(1+vmax/light)*nsamp/coeff1) ; maximum oversampled pixel shift

    if (pmax-pmin) eq 1L then message, 'Problem!'
    
    newaxis = long(1.0D + nsamp*(alog10(maxwave)-alog10(minwave))/coeff1)
    newloglam = coeff0 + coeff1*findgen(newaxis)/nsamp
    newlam = 10.0^newloglam

; resample the data spectrum

    combine1fiber, alog10(wave), spec, specivar, newloglam=newloglam, $
      newflux=logspec, newivar=logspecivar
;   logmask = fix(logspec*0.0+1)

; resample the reference spectrum

    combine1fiber, alog10(refwave), refspec, newloglam=newloglam, newflux=logrefspec

    logspec = logspec/median(logspec)
    logspecivar = logspec*median(logspec)^2.0
    logrefspec = logrefspec/median(logrefspec)

    mkhdr, hdr, logspec
    sxaddpar, hdr, 'COEFF0', coeff0
    sxaddpar, hdr, 'COEFF1', coeff1

    result = zfind(logspec,logspecivar,hdr=hdr,starflux=reform(logrefspec,newaxis,1),$
      starloglam0=coeff0,npoly=npoly,zguess=0.0,zmin=zmin,zmax=zmax,doplot=0)

    zans.zshift = result.z
    zans.zshift_err = result.z_err
    zans.vshift = result.z*light
    zans.vshift_err = result.z_err*light

    if (result.z_err lt 0.0) then zans.errflag = result.z_err

    if keyword_set(doplot) and (zans.zshift_err ne -1.0) then begin

       if !d.window ne 2L then window, 2, xs=600, ys=400

       djs_plot, newlam, logrefspec, ps=10, xsty=3, ysty=3, xrange=[minwave,maxwave], $
;      djs_plot, refwave, refspec, ps=10, xsty=3, ysty=3, xrange=[minwave,maxwave], $
         thick=2.0, xthick=2.0, ythick=2.0, charsize=1.2, charthick=2.0, $
         xtitle='Wavelength', ytitle='Relative Flux'
       djs_oplot, newlam*(1+zans.zshift), logspec, ps=10, color='red', thick=1.0, line=0
;      djs_oplot, wave*(1+zans.zshift), spec, ps=10, color='red', thick=2.0
       legend, ['Reference Spectrum','Data Spectrum'], linestyle=[0,2], charsize=1.0, $
         charthick=2.0, color=djs_icolor(['default','red']), /right, /top, box=0, $
         thick=2.0
       legend, textoidl('\chi^{2}')+' = '+string(result.rchi2diff,format='(F0.0)'), $
         /left, /top, box=0, charsize=1.0, charthick=2.0
       
;      djs_plot, newlam, logrefspec, ps=10, xsty=3, ysty=3, xrange=[minwave,maxwave], $
;      djs_plot, refwave, refspec, ps=10, xsty=3, ysty=3, xrange=[minwave,maxwave], $
;        thick=2.0, xthick=2.0, ythick=2.0, charsize=1.2, charthick=2.0, $
;        xtitle='Wavelength', ytitle='Relative Flux', position=[0.13,0.5,0.95,0.95]
;      djs_oplot, newlam*(1+zans.zshift), logspec, ps=10, color='red', thick=1.0, line=2
;      djs_oplot, wave*(1+zans.zshift), spec, ps=10, color='red', thick=2.0
;      legend, ['Reference Spectrum','Data Spectrum'], linestyle=[0,2], charsize=1.0, $
;        charthick=2.0, color=djs_icolor(['default','red']), /right, /top, box=0, $
;        thick=2.0
       
;      djs_plot, lags, chi2array, ps=-4, xsty=3, ysty=3, thick=2.0, xthick=2.0, $
;        ythick=2.0, charsize=1.2, charthick=2.0, xtitle='Pixel shift', $
;        ytitle='\chi^{2}', position=[0.13,0.1,0.95,0.4], /noerase
;      cc = get_kbrd(1)
       
    endif

return, zans
end
