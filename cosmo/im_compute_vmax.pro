;+
; NAME:
;       IM_COMPUTE_VMAX
;
; PURPOSE:
;       Given the limits of a survey and the results of running
;       k-correct, compute ZMIN, ZMAX, VMAX, and V/VMAX.
;
; INPUTS:
;       zobj          - galaxy redshifts [NGAL]
;       mabs          - absolute magnitude in FILTER [NGAL]
;       kcorr_coeffs  - best-fit coefficients from k-correct [5,NGAL]
;       survey_filter - filter bandpasses in which the survey is limited 
;       survey_mlimit - magnitude limits of the survey (e.g.,
;                       [14.5,17.77] in r for the SDSS)
;
; OPTIONAL INPUTS:
;       survey_zmin - minimum sample redshift [default min(ZOBJ)]
;       survey_zmax - maximum sample redshift [default max(ZOBJ)]
;       survey_area - survey area [arcsec^2]
;       vname       - k-corrected templates (default
;                     'default.nolines') 
;
; KEYWORD PARAMETERS:
;       debug - generate a visualization plot for debugging purposes
;
; OUTPUTS:
;       result - output data structure with everything you need
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Jun 12, NYU - written, based on discussions
;         with L. Moustakas and M. Blanton
;
; Copyright (C) 2007, John Moustakas
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

function im_compute_vmax, zobj, mabs, mapp, kcorr_coeffs, survey_filter, survey_mlimit, $
  survey_zmin=survey_zmin, survey_zmax=survey_zmax, survey_area=survey_area, $
  vname=vname, debug=debug

    ngalaxy = n_elements(zobj)
    if (ngalaxy eq 0L) or (n_elements(kcorr_coeffs) eq 0L) or $
      (n_elements(survey_filter) eq 0L) or (n_elements(survey_mlimit) eq 0L) then begin
       doc_library, 'im_compute_vmax'
       return, -1L
    endif

    if (ngalaxy ne (size(kcorr_coeffs,/dim))[1]) then begin
       splog, 'ZOBJ and KCORR_COEFFS must have the same number of objects'
       return, -1L
    endif

    if (n_elements(survey_mlimit) ne 2L) then begin
       splog, 'SURVEY_MLIMIT must be a two-element array.'
       return, -1L
    endif
    
    if (n_elements(survey_zmin) eq 0L) then survey_zmin = min(zobj)
    if (n_elements(survey_zmax) eq 0L) then survey_zmax = max(zobj)
    if (n_elements(vname) eq 0L) then vname = 'default.nolines'

; load the filter look-up table
    
    survey_nz = 500L

    k_load_vmatrix, vmatrix, lambda, vname=vname
    k_projection_table, rmatrix, vmatrix, lambda, survey_zvals, survey_filter, $ 
      zmin=survey_zmin, zmax=survey_zmax, nz=survey_nz

    result = {zobj: zobj, mabs: mabs, kcorr_coeffs: kcorr_coeffs, survey_filter: survey_filter, $
      survey_mlimit: survey_mlimit, survey_zmin: survey_zmin, survey_zmax: survey_zmax, $
      survey_zvals: survey_zvals, survey_rmatrix: rmatrix, $
      vol: 0.0D, zmin: fltarr(ngalaxy), zmax: fltarr(ngalaxy), vmax: dblarr(ngalaxy)}

    result.vol = qpint1d('dvcomoving',survey_zmin,survey_zmax,functargs={mpc:1L}) ; total volume [Mpc^3]

    for igal = 0L, ngalaxy-1L do begin

       print, format='("Working on ",I0,"/",I0,".",A10,$)', igal+1, ngalaxy, string(13b)

; compute K(z)       
       
       k_reconstruct_maggies, kcorr_coeffs[*,igal], 0.0, restmaggies, $
         vmatrix=vmatrix, lambda=lambda, filterlist=survey_filter, /silent
       k_reconstruct_maggies, rebin(kcorr_coeffs[*,igal],5.0,survey_nz), survey_zvals, $
         survey_obsmaggies, vmatrix=vmatrix, lambda=lambda, filterlist=survey_filter, /silent
       kcorr = 2.5*alog10(restmaggies[0]/survey_obsmaggies) ; K(z)

; compute m(z); see eq. (17) in Blanton et al. (2003)       
       
       survey_appmag = mabs[igal] + dmodulus(survey_zvals) + kcorr
       
       result.zmin[igal] = interpol(survey_zvals,survey_appmag,survey_mlimit[0]) > survey_zmin
       result.zmax[igal] = interpol(survey_zvals,survey_appmag,survey_mlimit[1]) < survey_zmax  

       result.vmax[igal] = qpint1d('dvcomoving',result.zmin[igal],result.zmax[igal],functargs={mpc:1L})
       
       if keyword_set(debug) then begin

          plot, survey_zvals, survey_appmag, ps=-4, xsty=3, ysty=3, charsize=1.8, charthick=2.0, $
            xtitle='Redshift', ytitle='Apparent Magnitude ['+survey_filter+']'
          oplot, !x.crange, survey_mlimit[0]*[1,1], thick=2, color=djs_icolor('red')
          oplot, !x.crange, survey_mlimit[1]*[1,1], thick=2, color=djs_icolor('red')
          oplot, !x.crange, mapp[igal]*[1,1]
          oplot, result.zmax[igal]*[1,1], !y.crange, color=djs_icolor('blue')
          oplot, result.zmin[igal]*[1,1], !y.crange, color=djs_icolor('blue')
          oplot, zobj[igal]*[1,1], !y.crange

          print, mabs[igal], mapp[igal]-interpol(survey_appmag,survey_zvals,zobj[igal]), $
            zobj[igal], result.zmin[igal], result.zmax[igal]
          cc = get_kbrd(1)

       endif
          
    endfor

return, result
end 
