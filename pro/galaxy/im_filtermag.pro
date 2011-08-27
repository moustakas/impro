;+
; NAME:
;   IM_FILTERMAG()
;
; PURPOSE:
;   Compute the AB mag of a spectrum in a given filter. 
;
; INPUTS:
;   wave - wavelength array [Angstrom]
;   flux_flam - flux array [ergs/s/cm2/A]
;
; OPTIONAL INPUTS:
;   filterlist - filter or list of filters to convolve with the
;     spectrum (analogous to KCORRECT)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   maggies - flux density in erg/s/cm2/Hz; the AB magnitude is given
;     by -2.5*alog10(maggies) 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The spectrum can be on an irregular wavelength grid.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2006 July 14, U of A, stolen from C. Tremonti
;      and modified to be more like K_PROJECT_FILTERS()
;   jm08mar12nyu - documented 
;   jm10feb17ucsd - cleaned up
;
; Copyright (C) 2006, 2008, 2010, John Moustakas
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

function im_filtermag, wave, flux_flam, filterlist=filterlist

    npix = n_elements(wave)
    nfilt = n_elements(filterlist)
    if (npix eq 0) or (nfilt eq 0) then begin
       doc_library, 'im_filtermag'
       return, -1
    endif
    if (npix ne n_elements(flux_flam)) then begin
       splog, 'Dimensions of WAVE and FLUX_FLAM must match'
       return, -1
    endif

; load the filters    
    k_load_filters, filterlist, k_filter_nlambda, $
      k_filter_lambda, k_filter_pass
    filtinfo = {filtw: k_filter_lambda, filtf: k_filter_pass, $
      filtn: k_filter_nlambda}

; call this routine recursively    
    if (nfilt gt 1) then begin
       maggies = fltarr(nfilt)
       for i = 0L, nfilt-1L do maggies[i] = im_filtermag(wave,flux_flam,$
         filtinfo=filtinfo[i])
       return, maggies
    endif

    filter_wave = filtinfo.filtw[0:filtinfo.filtn-1]
    filter_trans = filtinfo.filtf[0:filtinfo.filtn-1]

; determine width of each pixel in log-lambda units
    pix = findgen(npix)
    logwave = alog10(wave)
    logdiff = abs(logwave[1:*] - logwave) ; difference from 1 pixel center to the next
    linterp, pix[0:npix-2]+0.5, logdiff, pix, pixwidth 

; interpolate filter transmission to match the wavelength array of the
; input spectrum
    filter_sort = sort(filter_wave)
    linterp, filter_wave[filter_sort], filter_trans[filter_sort], $
      wave, filter_img, missing=0.0

; convert the input SED from [erg/s/cm^2/A] to [erg/s/cm^2/Hz]
    flux_fnu = flux_flam*wave^2/2.99792458D18

; add up the flux transmitted by the filter; NOTE!  We are really
; summing f_nu * dnu / nu * filter_trans; the trick is that
; dlog-lambda = dnu / nu / ln(10); the reason that each pixel is
; divided by nu is that we want to integrate in units of photons
; rather than units of energy
    filter_norm = total(filter_img*pixwidth)
    filter_flux = total(flux_fnu*filter_img*pixwidth)/filter_norm

    maggies = filter_flux*10.0^(0.4*48.6)
;   filter_abmag = -2.5 * alog10(filter_flux) - 48.6

return, maggies
end

