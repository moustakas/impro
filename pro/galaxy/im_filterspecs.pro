;+
; NAME:
;       IM_FILTERSPECS()
;
; PURPOSE:
;       Compute and return useful filter specifications. 
;
; INPUTS:
;       None
;
; OPTIONAL INPUTS:
;       filterlist - list of filters (see KCORRECT)
;       extra      - extra keywords for K_VEGA2AB()
;
; KEYWORD PARAMETERS:
;       verbose - print the parameters to the screen
;       doplot  - show the filter curve
;
; OUTPUTS:
;       specs - output data structure
;
; OPTIONAL OUTPUTS:
;       None
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Oct 12, U of A
;       jm05mar03uofa - added VEGA_FNU, VEGA_FLAM, and VEGA_JY tags;
;                       added VERBOSE keyword
;       jm05jun01uofa - added _EXTRA keyword for K_VEGA2AB()
;       jm05jun07uofa - added DOPLOT keyword
;       jm08mar12nyu  - documentation cleanup
;
; Copyright (C) 2004-2005, 2008, John Moustakas
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

function im_filterspecs, filterlist=filterlist, verbose=verbose, $
  doplot=doplot, _extra=extra

    light = 2.99792458D18       ; speed of light [A/s]

    if (n_elements(filterlist) eq 0L) then $
      filterlist = [$
        'bessell_'+['U','B','V','R','I'],$
        'sdss_'+['u0','g0','r0','i0','z0'],$
        'twomass_'+['J','H','Ks']]+'.par'
    nfilt = n_elements(filterlist)

    vega2ab = k_vega2ab(filterlist=filterlist,_extra=extra,/kurucz,/silent)
    solarmags = k_solar_magnitudes(filterlist=filterlist,/silent)
    weff = k_lambda_eff(filterlist=filterlist)

;   Umsun = +5.66 ; from Cox (2000) and Worthy 1994
;   Bmsun = +5.47
;   Vmsun = +4.82
;   Rmsun = +4.46
;   Imsun = +4.14

;   umsun_sdss = +6.41 ; all from Bell et al. (2003)
;   gmsun_sdss = +5.15
;   rmsun_sdss = +4.67
;   imsun_sdss = +4.56
;   zmsun_sdss = +4.53

;   Jmsun = +3.70  ; 2MASS, Johnson
;   Hmsun = +3.37  ; 2MASS, Johnson
;   Ksmsun = +3.32 ; 2MASS, from Bell et al. (2003)

; load the filter curves

    k_load_filters, filterlist, filter_nlambda, filter_lambda, filter_pass, _extra=extra
    nfilterpts = max(filter_nlambda)

; initialize the output data structure    
    
    specs = {$
      filter:     '', $
      vega_fnu:  0.0, $ ; Vega flux [erg/s/cm2/Hz] 
      vega_flam: 0.0, $ ; Vega flux [erg/s/cm2/A] 
      vega_Jy:   0.0, $ ; Vega flux [Jy] 
      sdssflag:   0L, $ ; SDSS filter true/false
      vega2ab:   0.0, $
      solarmags: 0.0, $
      weff:      0.0, $
      fwhm:      0.0, $
      filtw:     fltarr(nfilterpts), $
      filtf:     fltarr(nfilterpts), $
      filtn:       0L}
    specs = replicate(specs,nfilt)

; convert the solar magnitudes in non-SDSS filters to Vega magnitudes 

    sdssflag = lonarr(nfilt)+1L
    notsdss = where(strmatch(filterlist,'*sdss*',/fold) eq 0B,nnotsdss)
    if (nnotsdss ne 0L) then begin
       solarmags[notsdss] = solarmags[notsdss] - vega2ab[notsdss]
       sdssflag[notsdss] = 0L
    endif
    
    specs.filter = filterlist
    specs.sdssflag = sdssflag
    specs.vega2ab = vega2ab
    specs.solarmags = solarmags
    specs.weff = weff
    specs.vega_fnu = 3.631D-20*10^(-0.4*specs.vega2ab)  ; [erg/s/cm2/Hz]
    specs.vega_flam = specs.vega_fnu*light/specs.weff^2 ; [erg/s/cm2/A]
    specs.vega_Jy = specs.vega_fnu*1D23
    
; store the filter curves and compute the bandpass width 

    specs.filtw = filter_lambda
    specs.filtf = filter_pass
    specs.filtn = filter_nlambda

    if keyword_set(doplot) then window, 0, xs=500, ys=500
    
    for ifilt = 0L, nfilt-1L do begin

       good = lindgen(specs[ifilt].filtn)

       filtweff = specs[ifilt].weff
       filtw = specs[ifilt].filtw[good]
       filtf = specs[ifilt].filtf[good]
       
; refine the wavelength grid

       minw = min(filtw)
       maxw = max(filtw)
       dw = 1.0D                ; 1 A pixels

       ww = minw+dw*findgen((maxw-minw+dw)/dw)
       ff = interpol(filtf,filtw,ww)
       filtarea = int_tabulated(ww,ff,/sort) ; filter area
       junk = check_math()
       
       specs[ifilt].fwhm = 2.0*sqrt(int_tabulated(ww,ff*(ww-filtweff)^2.0)/filtarea) ; [Angstrom]
       
       if keyword_set(doplot) then begin
          if specs[ifilt].weff gt 1D4 then norm = 1D4 else norm = 1D
          plot, ww/norm, ff, xsty=3, ysty=3, xthick=2, ythick=2, title=strupcase(specs[ifilt].filter), $
            thick=2, charsize=1.4, charthick=2.0, xtitle='Wavelength', $
            ytitle='Response'
          cc = get_kbrd(1)

       endif
       
    endfor

    if keyword_set(verbose) then struct_print, struct_trimtags(specs,except=['FILTW','FILTF','FILTN'])
    
return, specs
end
    
