;+
; NAME:
;   IABSLINEEW()
;
; PURPOSE:
;   Measure absorption-line equivalent widths or continuum indices
;   using pseudo-continuum side-bands. 
;
; INPUTS:
;   wave     - wavelength vector [NPIX]
;   flux     - data spectrum [NPIX]
;   linewave - line wavelengths [NLINE]; if NOLINE=1 then 
;              LINEWAVE can be an array of zeros
;
; OPTIONAL INPUTS:
;   ivar      - corresponding inverse variance spectrum [NPIX] (set
;     IVAR=0 for crummy pixels)
;   llimit   - central wavelength of the lower continuum window
;              (default LINEWAVE minus 20 Angstrom)
;   lwidth   - full width of LLIMIT (default 15 Angstrom) 
;   ulimit   - central wavelength of the upper continuum window
;              (default LINEWAVE plus 20 Angstrom)
;   uwidth   - full width of ULIMIT (default 15 Angstrom) 
;   lline    - wavelength of the lower (left) part of the line
;              (default LINEWAVE minus 5 Angstrom) 
;   uline    - wavelength of the upper (right) part of the line
;              (default LINEWAVE plus 5 Angstrom)
;   ncoeff   - number of coefficients of the fit between the two
;              continuum side-bands (default 2, linear)
;   label    - title for the plot if DEBUG=1
;
; KEYWORD PARAMETERS:
;   noline     - do not compute the absorption-line EW (just
;                compute the continuum side-bands and their ratio)  
;   debug      - generate a diagnostic plot and wait for a
;                keystroke 
;   magnitude  - compute the magnitude of the absorption-line EW 
;   fnu        - FLUX is in units of [erg/s/cm2/Hz]
;   postscript - if set then do not open a window and suppress
;                waiting for a keystroke between plots
;   keystroke  - wait for a keystroke
;   silent     - do not print warning messages to STDOUT 
;
; OUTPUTS:
;   absline   - output data structure
;
; OPTIONAL OUTPUTS:
;   absplot   - data structure with parameters used in the
;               plotting
;
; COMMENTS:
;
; PROCEDURE:
;   The continuum windows are defined by LLIMIT, ULIMIT, LWIDTH,
;   and UWIDTH.  The mean continuum in each window is computed and
;   an NCOEFF function is fitted between the points.  The
;   continnum is interpolated at the central wavelength of the
;   line. 
;
; EXAMPLE:
;
; PROCEDURES USED:
;   PARTIAL_PIXEL_SUM(), GET_ELEMENT, DJS_PLOT, DJS_OPLOT, SPLOG,
;   FLAM_UNITS(), FNU_UNITS()
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 November 26, U of A
;   jm03dec1uofa - added POSTSCRIPT keyword
;   jm04jan5uofa - additional error checking; added FNU keyword 
;   jm04apr22uofa - plotting improvements
;   jm04oct25uofa - added SILENT keyword
;   jm05jul28uofa - added ALLOW_PARTIAL keyword 
;   jm09dec08ucsd - changed FERR optional input to IVAR (inverse
;     variance), so that bad pixels could be tracked (IVAR=0)
;
; Copyright (C) 2003-2009, John Moustakas
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

pro qaplot_abslineew, wave, flux, mask=mask, absline=absline, $
  absplot=absplot, fnu=fnu, noline=noline, postscript=postscript

    npix = n_elements(wave)
    
; make a plot    
    if keyword_set(postscript) then $
      colorlist = ['black','light blue','light red','black','black','orange'] else $
      colorlist = ['yellow','blue','red','white','cyan','light green']

    plotsym, 8, 2, /fill
    medflux = abs(median(flux))
;   if (medflux le 0.0) then begin
;      splog, 'Median flux negative - no QAplot possible!'
;      return
;   endif
    power = ceil(abs(alog10(medflux)))
    scale = 10.0^power
    if keyword_set(fnu) then units = fnu_units() else units = flam_units()
    ytitle = 'Flux (10^{-'+string(power,format='(I0)')+'} '+units+')'

; compute the plotting range and plot the sub-spectrum
    wplot = 20.0
    get_element, wave, [absline.llimit,absline.ulimit]+$
      [-absline.lwidth,+absline.uwidth]/2+wplot*[-1,1], ww
    absplot.wloindx = ww[0]
    absplot.whiindx = ww[1]

    pixoff = 20
    xrange = fltarr(2)
    xrange[0] = interpol(wave,findgen(npix),ww[0]-pixoff)
    xrange[1] = interpol(wave,findgen(npix),ww[1]+pixoff)

    yrange = minmax(flux[(ww[0]-pixoff)>0L:(ww[1]+pixoff)<(npix-1)])*[1.0,1.02]

; initialize the plot       
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, charsize=1.8, $
      xtitle='Rest Wavelength (\AA)', ytitle=ytitle, $
      xrange=xrange, yrange=yrange*scale, _extra=extra
    legend, textoidl(absplot.label), /right, /top, box=0, charsize=2.0

; overplot the lower continuum window
    bb = [absplot.lloindx,absplot.lhiindx]
    polyfill, [wave[bb[0]],wave[bb[1]],wave[bb[1]],wave[bb[0]]], $ ; blue
      [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
      /line_fill, orientation=135, color=djs_icolor(colorlist[1]), spacing=0.1, linestyle=1
    polyfill, [wave[bb[0]],wave[bb[1]],wave[bb[1]],wave[bb[0]]], $ ; blue
      [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
      /line_fill, orientation=45, color=djs_icolor(colorlist[1]), spacing=0.1, linestyle=1
    
; overplot the upper continuum window
    rr = [absplot.uloindx,absplot.uhiindx]
    polyfill, [wave[rr[0]],wave[rr[1]],wave[rr[1]],wave[rr[0]]], $ ; red
      [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
      /line_fill, orientation=135, color=djs_icolor(colorlist[2]), spacing=0.1, linestyle=1
    polyfill, [wave[rr[0]],wave[rr[1]],wave[rr[1]],wave[rr[0]]], $ ; red
      [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
      /line_fill, orientation=45, color=djs_icolor(colorlist[2]), spacing=0.1, linestyle=1

; overplot the line window       
    if (keyword_set(noline) eq 0) then begin
       ll = [absplot.llindx,absplot.luindx]
       polyfill, [wave[ll[0]],wave[ll[1]],wave[ll[1]],wave[ll[0]]], $ ; line
         [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
         /line_fill, orientation=135, color=djs_icolor(colorlist[5]), $
         spacing=0.1, linestyle=1
       polyfill, [wave[ll[0]],wave[ll[1]],wave[ll[1]],wave[ll[0]]], $ ; line
         [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
         /line_fill, orientation=45, color=djs_icolor(colorlist[5]), $
         spacing=0.1, linestyle=1
    endif

; now plot the spectrum over the continuum and line windows
    djs_oplot, wave[(ww[0]-pixoff)>0L:(ww[1]+pixoff)<(npix-1)], $
      scale*flux[(ww[0]-pixoff)>0L:(ww[1]+pixoff)<(npix-1)], $
      ps=10, color=colorlist[0]
    if (mask[0] ne -1) then for jj = 0, n_elements(mask)-1 do $
      plots, wave[mask[jj]], scale*flux[mask[jj]], psym=7, symsize=2.0, $
      color=djs_icolor('dark green')

; finally overplot the continuum points and the linear fit between the
; two pseudo-continua 
    plots, absline.llimit, scale*absline.lcflux, ps=8, $
      color=djs_icolor(colorlist[3]), symsize=1.5
    plots, absline.ulimit, scale*absline.ucflux, ps=8, $
      color=djs_icolor(colorlist[3]), symsize=1.5

    cwave = wave[bb[0]:rr[1]]
    cflux = poly(cwave-absline.wave0,absline.ccoeff)
    djs_oplot, cwave, scale*cflux, line=5, color=colorlist[4]

    if keyword_set(keystroke) then cc = get_kbrd(1)

return
end

function init_abslineew, linewave=linewave, llimit=llimit, $
  lwidth=lwidth, ulimit=ulimit, uwidth=uwidth, lline=lline, $
  uline=uline, ncoeff=ncoeff, label=label, forplot=forplot
; initialize the output and plotting data structures
    if keyword_set(forplot) then begin
       info = {$
         label:  label,$
         lloindx:    0,$        ; lower continuum lo
         lhiindx:    0,$        ; lower continuum hi
         uloindx:    0,$        ; upper continuum lo
         uhiindx:    0,$        ; upper continuum hi
         llindx:     0,$        ; line lo
         luindx:     0,$        ; line hi
         wloindx:    0,$
         whiindx:    0}
    endif else begin
       info = {$
         linewave:      linewave,       $
         lineew:        0.0D,           $ ; line EW
         lineew_err:   -2.0D,           $ ; line EW error
         lineflux:      0.0D,           $ ; line flux
         lineflux_err: -2.0D,           $ ; line flux error
         linec:         0.0D,           $ ; interpolated continuum at line center
         linec_err:    -2.0D,           $ ; interpolated continuum error at line center
         llimit:        llimit,         $
         lwidth:        lwidth,         $
         ulimit:        ulimit,         $
         uwidth:        uwidth,         $
         lline:         lline,          $
         uline:         uline,          $
;        lcpartial:        0,           $ ; 1 means part of the band is outside of the wavelength range
         lcflux:        0.0D,           $ ; lower continuum flux
         lcflux_err:   -2.0D,           $ ; lower continuum flux error
;        ucpartial:        0,           $ ; 1 means part of the band is outside of the wavelength range
         ucflux:        0.0D,           $ ; upper continuum flux
         ucflux_err:   -2.0D,           $ ; upper continuum flux error
         cratio:        0.0D,           $ ; ratio of the UPPER to LOWER continuum flux
         cratio_err:   -2.0D,           $ ; error in CRATIO
         magnitude:        0,           $ ; Yes/No
         wave0:         0.0D,           $ ; mean wavelength
         ccoeff:        dblarr(ncoeff), $ ; coefficients of the fit between LCFLUX and UCFLUX
         ccoeff_err:    dblarr(ncoeff)-2} ; error in CCOEFF (assumed zero)
    endelse
return, info
end

function iabslineew, wave, flux, linewave, ivar=ivar1, llimit=llimit, $
  lwidth=lwidth, ulimit=ulimit, uwidth=uwidth, lline=lline, uline=uline, $
  ncoeff=ncoeff, wplot=wplot, label=label, absplot=absplot, noline=noline, $
  debug=debug, magnitude=magnitude, fnu=fnu, allow_partial=allow_partial, $
  postscript=postscript, keystroke=keystroke, silent=silent, _extra=extra

    npix = n_elements(wave)
    nflux = n_elements(flux)
    nline = n_elements(linewave)

    if (npix eq 0) or (nflux eq 0) or (nline eq 0) then begin
       doc_library, 'iabslineew'
       return, -1
    endif

; prepare the IVAR array    
    if (n_elements(ivar1) eq 0) then begin
       sig = djsig(flux)
       if (sig le 0) then sig = 1.0
       ivar = flux*0.0+1.0/sig^2
    endif else ivar = ivar1

    neg = where(ivar lt 0,nneg)
    if (nneg ne 0L) then message, 'IVAR array contains negative values'
    nivar = n_elements(ivar)

    if (npix ne nflux) or (npix ne nivar) then begin
       print, 'Dimensions of WAVE, FLUX, and IVAR do not agree'
       return, -1
    endif

    if keyword_set(postscript) then debug = 0

; set some respectable defaults    
    if (nline eq 1) then linewave = linewave[0]
    if n_elements(llimit) eq 0 then llimit = linewave-20.0
    if n_elements(ulimit) eq 0 then ulimit = linewave+20.0
    if n_elements(lwidth) eq 0 then lwidth = linewave*0.0+15
    if n_elements(uwidth) eq 0 then uwidth = linewave*0.0+15
    if n_elements(lline) eq 0 then lline = linewave-5.0
    if n_elements(uline) eq 0 then uline = linewave+5.0
    if n_elements(label) eq 0 then label = replicate('',nline)
    if n_elements(ncoeff) eq 0 then ncoeff = 2

    if (nline gt 1) then begin
       for k = 0, nline-1 do begin
          absline1 = iabslineew(wave,flux,linewave[k],ivar=ivar,llimit=llimit[k],$
            lwidth=lwidth[k],ulimit=ulimit[k],uwidth=uwidth[k],lline=lline[k],$
            uline=uline[k],ncoeff=ncoeff,label=label[k],_extra=extra,$
            absplot=absplot1,noline=noline,debug=debug,magnitude=magnitude,$
            fnu=fnu,allow_partial=allow_partial,postscript=postscript,silent=silent)
          if (k eq 0) then begin
             absline = absline1 
             absplot = absplot1
             if keyword_set(debug) and (absline1.cratio_err gt 0.0) then $
               splog, 'Press any key to continue'
          endif else begin 
             absline = [[absline],[absline1]]
             absplot = [[absplot],[absplot1]]
          endelse
          if keyword_set(debug) and (absline1.cratio_err gt 0.0) then $
            cc = get_kbrd(1)
       endfor
       absplot = reform(absplot)
       return, reform(absline)
    endif else begin
       llimit = llimit[0]
       ulimit = ulimit[0]
       lwidth = lwidth[0]
       uwidth = uwidth[0]
       lline = lline[0]
       uline = uline[0]
    endelse

; initialize the output structures
    absline = init_abslineew(linewave=linewave,llimit=llimit,$
      lwidth=lwidth,ulimit=ulimit,uwidth=uwidth,lline=lline, $
      uline=uline,ncoeff=ncoeff)
    absplot = init_abslineew(label=label,/forplot)
    absline.wave0 = (djs_mean([absline.llimit,absline.ulimit])>min(wave))<max(wave)

; build the error array and deal with masked pixels
    ferr = sqrt(1.0/(ivar+(ivar eq 0))*(ivar ne 0))
    good = where(ivar gt 0,ngood,comp=mask,ncomp=nmask)
    if (ngood eq 0) then message, 'All pixels masked!'
    
; compute the lower continuum
    lowave = absline.llimit - absline.lwidth/2.0
    hiwave = absline.llimit + absline.lwidth/2.0
    if (lowave ge min(wave)) and (hiwave le max(wave)) then begin
       absplot.lloindx = findex(wave,lowave)
       absplot.lhiindx = findex(wave,hiwave)
       absline.lcflux = im_integral(wave[good],flux[good],lowave,hiwave)/(hiwave-lowave)
       absline.lcflux_err = sqrt(im_integral(wave[good],ferr[good]^2,lowave,hiwave))/(hiwave-lowave)
    endif else absline.lcflux = -2.0

; compute the upper continuum
    lowave = absline.ulimit - absline.uwidth/2.0
    hiwave = absline.ulimit + absline.uwidth/2.0
    if (lowave ge min(wave)) and (hiwave le max(wave)) then begin
       absplot.uloindx = findex(wave,lowave)
       absplot.uhiindx = findex(wave,hiwave)
       absline.ucflux = im_integral(wave[good],flux[good],lowave,hiwave)/(hiwave-lowave)
       absline.ucflux_err = sqrt(im_integral(wave[good],ferr[good]^2,lowave,hiwave))/(hiwave-lowave)
       absplot.uloindx = findex(wave,lowave) ; for the plot
       absplot.uhiindx = findex(wave,hiwave)
    endif else absline.ucflux = -2.0

; if both are defined, fit a polynomial between the two
; pseudo-continuum points; to get the errors in the coefficients right
; subtract off the starting wavelength

    if (absline.lcflux_err gt 0.0) and (absline.ucflux_err gt 0.0) then begin
       absline.ccoeff = poly_fit([absline.llimit,absline.ulimit]-absline.wave0,$
         [absline.lcflux,absline.ucflux],ncoeff-1,measure_errors=$
         [absline.lcflux_err,absline.ucflux_err],sigma=ccoeff_err,$
         /double,chisq=chisq,covar=covar,yerror=yerror,yband=yband,$
         yfit=yfit,status=status)
       absline.ccoeff_err = ccoeff_err
    endif
    
;; RELEGATED (do not allow)! if just one of the two continuum bands is
;; defined (because one of the two bands falls outside the wavelength
;; range of the spectrum), then set the missing continuum to the
;; existing continuum level, effectively assuming a flat slope over the
;; line
;    if keyword_set(allow_partial) then begin
;       if (absline.lcflux ne -2.0) and (absline.ucflux eq -2.0) then begin
;          absline.ucflux     = absline.lcflux
;          absline.ucflux_err = absline.lcflux_err
;          absline.ccoeff     = [absline.lcflux,0.0] ; slope of zero
;          absline.ccoeff_err = [absline.lcflux_err,0.0]
;       endif
;       if (absline.lcflux eq -2.0) and (absline.ucflux ne -2.0) then begin
;          absline.lcflux     = absline.ucflux
;          absline.lcflux_err = absline.ucflux_err
;          absline.ccoeff     = [absline.ucflux,0.0] ; slope of zero
;          absline.ccoeff_err = [absline.ucflux_err,0.0]
;       endif
;    endif
       
; compute the ratio of the continuum side-bands and the error
    if (absline.lcflux_err gt 0.0) and (absline.ucflux_err gt 0.0) then begin
       absline.cratio = (absline.lwidth/absline.uwidth)*(absline.ucflux/absline.lcflux)
       absline.cratio_err = (absline.lwidth/absline.uwidth)*$
         sqrt((absline.ucflux_err/absline.lcflux)^2 + $
         (absline.lcflux_err*absline.ucflux/absline.lcflux^2)^2)
    endif

; compute the line EW; evaluate the continuum over the absorption line
; and continuum points; see Trager et al. (1998)
    lowave = absline.lline
    hiwave = absline.uline
    if (keyword_set(noline) eq 0) and (absline.lcflux_err gt 0.0) and $
      (absline.ucflux_err gt 0.0) and (lowave ge min(wave)) and $
      (hiwave le max(wave)) then begin
       absplot.llindx = findex(wave,lowave)
       absplot.luindx = findex(wave,hiwave)
       
       clineflux = poly(wave-absline.wave0,absline.ccoeff)
       clineflux_err = sqrt(poly((wave-absline.wave0)^2,absline.ccoeff_err^2)) ; continuum error
       
       absline.lineflux = im_integral(wave[good],(clineflux-flux)[good],lowave,hiwave)/(hiwave-lowave)
       absline.lineflux_err = sqrt(im_integral(wave[good],ferr[good]^2,lowave,hiwave))/(hiwave-lowave)
       absline.linec = interpol(clineflux,wave,linewave)
       absline.linec_err = interpol(clineflux_err,wave,linewave)

       index = 0.0
       index_err = -2.0
       if keyword_set(magnitude) then begin 
          absline.magnitude = 1
          lineflux = flux/clineflux
       endif else lineflux = 1.0 - flux/clineflux
       lineflux_err = ferr/clineflux
          
       index = im_integral(wave[good],lineflux[good],lowave,hiwave)
       index_err = sqrt(im_integral(wave[good],lineflux_err[good]^2,lowave,hiwave))

       if keyword_set(magnitude) and (index gt 0.0) then begin
          index_err = (2.5/alog(10.0))*index_err/index
          index = -2.5*alog10(index)
       endif

       absline.lineew = index
       absline.lineew_err = index_err
    endif 

; build a QAplot if requested
    if keyword_set(debug) and (absline.cratio_err gt 0.0) then begin
       qaplot_abslineew, wave, flux, mask=mask, absline=absline, $
         absplot=absplot, fnu=fnu, noline=noline, postscript=postscript
    endif
    
return, absline
end
