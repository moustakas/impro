;+
; NAME:
;   IM_GALEX_TO_MAGGIES
;
; PURPOSE:
;   Convert an input GALEX photometric catalog to Galactic
;   extinction-corrected AB maggies. 
;
; INPUTS: 
;   galex - input catalog [NGAL]
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;   allow_nuv_nondetect - allow photometry of objects not detected in
;     the NUV, but with measured FUV photometry (in general objects
;     without solid NUV detections should not be used) (has no effect
;     if /PSF_FLUX)
;   psf_flux - use the PSF flux from the Schiminovich catalogs (default is
;     to use the AUTO flux) 
;
; OUTPUTS: 
;   maggies - [2,NGAL] output FUV/NUV maggies 
;   ivarmaggies - [2,NGAL] corresponding inverse variance
;
; OPTIONAL OUTPUTS:
;   filterlist - GALEX filterlist names
;
; COMMENTS:
;   A minimum photometric error of [0.052,0.026] for [FUV,NUV] is
;   applied, as recommended by Morrissey+07.  Uses the MAG_AUTO
;   photometry by default.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Apr 30, UCSD - based loosely on
;     M. Blanton's GALEX_TO_MAGGIES
;   jm11aug30ucsd - added ALLOW_NUV_NONDETECT keyword
;
; Copyright (C) 2010-2011, John Moustakas
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

pro im_galex_to_maggies, galex, maggies, ivarmaggies, psf_flux=psf_flux, $
  filterlist=filterlist, allow_nuv_nondetect=allow_nuv_nondetect

    ngal = n_elements(galex)
    if (ngal eq 0L) then begin
       doc_library, 'im_galex_to_maggies'
       return
    endif

    filterlist = galex_filterlist()
    nbands = n_elements(filterlist)

; correct for Galactic extinction; use the standard reddening values
; from Wyder+07; ignore the quadratic reddening term in the FUV
; channel, and also be sure to check for non-detections
    kl = [8.24,8.20] ; [FUV,NUV]
    ebv = fltarr(ngal)
    if tag_exist(galex,'alpha_j2000') then begin
       good = where(galex.alpha_j2000 gt -900.0,ngood)
       ra = galex[good].alpha_j2000
       dec = galex[good].delta_j2000
    endif else begin
       good = where(galex.ra gt -900.0,ngood)
       ra = galex[good].ra
       dec = galex[good].dec
    endelse
    if (ngood ne 0L) then begin
       glactc, ra, dec, 2000.0, gl, gb, 1, /deg
       ebv[good] = dust_getval(gl,gb,/interp,/noloop)
    endif

; magnitude zeropoints: see http://galex.stsci.edu/doc/GI_Doc_Ops7.pdf 
    fuv_zpt = 18.82
    nuv_zpt = 20.08

; convert to maggies; ignore artifacts; require an NUV detection for
; the FUV photometry, unless /allow_nuv_nondetect
    obsmaggies = fltarr(2,ngal)-999.0
    obsmaggieserr = fltarr(2,ngal)-999.0

    if keyword_set(psf_flux) then begin
       fuvfluxtag = tag_indx(galex,'fuv_flux')
       fuvfluxerrtag = tag_indx(galex,'fuv_fluxerr')
       nuvfluxtag = tag_indx(galex,'nuv_flux')
       nuvfluxerrtag = tag_indx(galex,'nuv_fluxerr')
    endif else begin
       fuvfluxtag = tag_indx(galex,'fuv_flux_auto')
       fuvfluxerrtag = tag_indx(galex,'fuv_fluxerr_auto')
       nuvfluxtag = tag_indx(galex,'nuv_flux_auto')
       nuvfluxerrtag = tag_indx(galex,'nuv_fluxerr_auto')
    endelse

    if keyword_set(psf_flux) then begin
; FUV
       good = where((galex.(fuvfluxtag) gt -90.0),ngood)
       if (ngood ne 0L) then begin
          obsmaggies[0,good] = galex[good].(fuvfluxtag)*10^(-0.4*23.9) ; microJy-->maggies
          obsmaggieserr[0,good] = galex[good].(fuvfluxerrtag)*10^(-0.4*23.9)
       endif
; NUV
       good = where((galex.(nuvfluxtag) gt -90.0),ngood)
       if (ngood ne 0L) then begin
          obsmaggies[1,good] = galex[good].(nuvfluxtag)*10^(-0.4*23.9) ; microJy-->maggies
          obsmaggieserr[1,good] = galex[good].(nuvfluxerrtag)*10^(-0.4*23.9)
       endif
    endif else begin
       good = where((galex.(nuvfluxtag) gt -90.0),ngood)
       if (ngood ne 0L) then begin
          obsmaggies[1,good] = galex[good].(nuvfluxtag)*10^(-0.4*nuv_zpt) ; galex flux-->maggies
          obsmaggieserr[1,good] = galex[good].(nuvfluxerrtag)*10^(-0.4*nuv_zpt)
; use the aperture-matched FUV photometry, if possible
          good_fuv = where(galex[good].fuv_ncat_flux gt -900.0,ngood_fuv)
          obsmaggies[0,good[good_fuv]] =  galex[good[good_fuv]].fuv_ncat_flux*10^(-0.4*23.9) ; microJy-->maggies
          obsmaggieserr[0,good[good_fuv]] =  galex[good[good_fuv]].fuv_ncat_fluxerr*10^(-0.4*23.9)
       endif

; optionally allow NUV non-detections       
       if keyword_set(allow_nuv_nondetect) then begin
          good = where((galex.(nuvfluxtag) lt -90.0) and (galex.(fuvfluxtag) gt -90.0),ngood)
          if (ngood ne 0L) then begin
             obsmaggies[0,good] = galex[good].(fuvfluxtag)*10^(-0.4*fuv_zpt) ; galex flux-->maggies
             obsmaggieserr[0,good] = galex[good].(fuvfluxerrtag)*10^(-0.4*fuv_zpt)
; use the aperture-matched NUV photometry, if possible
             good_nuv = where(galex[good].nuv_fcat_flux gt -900.0,ngood_nuv)
             obsmaggies[1,good[good_nuv]] =  galex[good[good_nuv]].nuv_fcat_flux*10^(-0.4*23.9) ; microJy-->maggies
             obsmaggieserr[1,good[good_nuv]] =  galex[good[good_nuv]].nuv_fcat_fluxerr*10^(-0.4*23.9)
          endif 
       endif 
    endelse
    
; now correct for extinction, convert to ivarmaggies, and return
    maggies = dblarr(2,ngal)
    ivarmaggies = dblarr(2,ngal)
    
    ngood = where((obsmaggies[1,*] gt -900.0),nngood)
    if (nngood ne 0L) then begin
       factor = 10^(+0.4*kl[1]*ebv[ngood])
       maggies[1,ngood] = obsmaggies[1,ngood]*factor
       ivarmaggies[1,ngood] = 1.0/(obsmaggieserr[1,ngood]*factor)^2.0
    endif

    fgood = where((obsmaggies[0,*] gt -900.0),nfgood)
    if (nfgood ne 0L) then begin
       factor = 10^(+0.4*kl[0]*ebv[fgood])
       maggies[0,fgood] = obsmaggies[0,fgood]*factor
       ivarmaggies[0,fgood] = 1.0/(obsmaggieserr[0,fgood]*factor)^2.0
    endif
    
; minimum photometric error from Morrissey+07, plus a little
    minerr = [0.052,0.026]
;   minerr = sqrt([0.052,0.026]^2 + [0.05,0.05]^2)
    k_minerror, maggies, ivarmaggies, minerr

return    
end
