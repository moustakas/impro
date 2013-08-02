;+
; NAME:
;   IM_ZMINZMAX()
;
; PURPOSE:
;   Compute the minimum and maximum redshift at which an object
;   could have been observed in a survey of specified flux limits,
;   given precomputed K-corrections.
;
; INPUTS: 
;   redshift - galaxy redshift [NGAL]
;   mag - apparent magnitude in the selection band [NGAL]
;   coeffs - coefficients calculated by KCORRECT [5,NGAL]
;   bright - bright magnitude cut-off of the survey
;   faint - faint magnitude cut-off of the survey
;   filter - filter function name corresponding to MAG
;
; OPTIONAL INPUTS: 
;   vname - basis set used to calculate K-corrections
;   q0, q1, qz0 - see K_EVOLVE()
;   h100, omega0, omegalambda - desired cosmology for RED; default to
;     h100=0.7, omega0=0.3, omegalambda=0.7
;   sample_zmin - absolute minimum redshift of the sample under
;     consideration (default 1E-3; see COMMENTS)
;   sample_zmax - absolute maximum redshift of the sample under
;     consideration (default 7.0; see COMMENTS)
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   Data structure with the following quantities:
;     REDSHIFT - input redshift
;     MAG - input apparent magnitude
;     ZMIN - redshift at which each object's brightness equals BRIGHT 
;     ZMAX - redshift at which each object's brightness equals FAINT 
;
; OPTIONAL OUTPUTS:
;   offset - magnitude offset applied to the model magnitudes to match
;     MAG at the redshift of the galaxy, presumably caused by
;     scatter/zeropoint differences between the selection magnitude
;     and the magnitude(s) used for K-corrections
;
; COMMENTS:
;   The apparent magnitude of a source at redshift, ztrue, is:
;     (1) m(obs) = M + DM(ztrue) + K(ztrue) - E(ztrue)
;   
;   where M is the absolute magnitude, DM(z) is the distance modulus,
;   K(z) is the K-correction, and E(z) is the evolutionary
;   correction (which is positive for q0>0).
; 
;   Similarly, the hypothetical apparent magnitude of the same source
;   at a different redshift would be:
;     (2) m(z) = M + DM(z) + K(z) - E(z)
;
;   By definition, BRIGHT and FAINT are the apparent magnitude at ZMIN
;   and ZMAX, respectively.  Therefore, subtracting equation (2) from
;   (1) we get: 
;     BRIGHT-m(obs) = [DM(zmin)-DM(ztrue)]+[K(zmin)-K(ztrue)]-[E(zmin)-E(ztrue)]
;     FAINT -m(obs) = [DM(zmax)-DM(ztrue)]+[K(zmax)-K(ztrue)]-[E(zmax)-E(ztrue)]
;
;   So we just need to compute DM(z), K(z), and E(z) on a grid of
;   redshifts and at the true redshift of each galaxy and then
;   interpolate to get ZMIN and ZMAX.  Note that one nice thing about
;   these equations is that all the zeropoints drop out as long as
;   BRIGHT, FAINT, and m(obs) are defined the same way (e.g., Vega,
;   uncorrected for reddening, etc.).  
;
;   Note that if evolution is allowed, then the delta-mag vs redshift
;   relationship may not be monotonic, and therefore ZMIN and ZMAX are
;   not defined.  In this case, try setting SAMPLE_ZMAX to a redshift
;   that is much smaller than the default (7.0), closer to the
;   absolute maximum redshift of the sample under consideration. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 07, UCSD
;   jm10sep09ucsd - SAMPLE_ZMAX optional input added
;
; Copyright (C) 2010, John Moustakas
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

function im_zminzmax, redshift, mag, coeffs, bright=bright, $
  faint=faint, filter=filter, offset=offset, vname=vname, $
  sample_zmin=sample_zmin, sample_zmax=sample_zmax, debug=debug, $
  q0=q0, q1=q1, qz0=qz0, h100=h100, omega0=omega0, $
  omegalambda=omegalambda

    ngal = n_elements(redshift)
    if (ngal eq 0) or (n_elements(mag) eq 0) or $
      (n_elements(coeffs) eq 0) then begin
       doc_library, 'im_zminzmax'
       return, -1
    endif

    if (ngal ne n_elements(mag)) then begin
       splog, 'Dimensions of REDSHIFT and NGAL must match'
       return, -1
    endif
    if (n_elements(filter) eq 0) then begin
       splog, 'FILTER input required'
       return, -1
    endif
    if (n_elements(bright) eq 0) or (n_elements(faint) eq 0) then begin
       splog, 'BRIGHT and FAINT inputs required'
       return, -1
    endif
    toofaint = where(mag gt faint,nfaint)
    if (nfaint ne 0) then begin
       splog, 'WARNING: Some objects are too faint given FAINT'
    endif
    toobright = where(mag lt bright,nbright)
    if (nbright ne 0) then begin
       splog, 'WARNING: Some objects are too bright given BRIGHT'
    endif

; cosmology    
    if (n_elements(h100) eq 0) then h100 = 0.7
    if (n_elements(omega0) eq 0) then omega0 = 0.3
    if (n_elements(omegalambda) eq 0) then omegalambda = 0.7
;   red, h100=h100, omega0=omega0, omegalambda=omegalambda

; evolutionary parameters    
    if (n_elements(q0) eq 0) then q0 = 0.0
    if (n_elements(q1) eq 0) then q1 = 0.0
    if (n_elements(qz0) eq 0) then qz0 = 0.0
    
; initialize the output structure    
    result = replicate({redshift: 0.0, mag: 0.0, $
      zmin: 0.0, zmax: 0.0},ngal)
    result.redshift = redshift
    result.mag = mag

; define the fiducial redshift grid    
    if (n_elements(sample_zmin) eq 0) then $
      zrefmin = 1E-3 else zrefmin = sample_zmin
    if (n_elements(sample_zmax) eq 0) then $
      zrefmax = 7.0 else zrefmax = sample_zmax
    nzref = 500
;   zref = range(zrefmin,zrefmax,nzref)
    zref = findgen(nzref)*(zrefmax-zrefmin)/(nzref-1)+zrefmin

; distance modulus and evolution
    dm_gal = lf_distmod(redshift,omega0=omega0,omegal0=omegalambda)
    dm_ref = lf_distmod(zref,omega0=omega0,omegal0=omegalambda)
;   dm_gal = dmodulus(redshift)
;   dm_ref = dmodulus(zref)
    evol_gal = k_evolve(redshift*0.0,redshift,q0,q1,qz0)
    evol_ref = k_evolve(zref*0.0,zref,q0,q1,qz0)

; compute the absolute magnitude of each source in the selection band
; on the fiducial redshift grid, allowing for luminosity evolution
    k_load_vmatrix, vmatrix, lambda, vname=vname
    k_projection_table, rmatrix, vmatrix, lambda, zvals, $
      filter, zmin=zrefmin, zmax=zrefmax, nz=nzref, /silent

    k_reconstruct_maggies, coeffs, redshift*0.0, $
      restmaggies, rmatrix=rmatrix, zvals=zvals
    k_reconstruct_maggies, coeffs, redshift, $
      maggies_gal, rmatrix=rmatrix, zvals=zvals
    kcorr_gal = +2.5*alog10(reform(restmaggies/maggies_gal))

    kcorr_ref = fltarr(nzref,ngal)
    for ii = 0L, nzref-1L do begin
       k_reconstruct_maggies, coeffs, redshift*0.0+zref[ii], $
         maggies, rmatrix=rmatrix, zvals=zvals
       kcorr_ref[ii,*] = +2.5*alog10(reform(restmaggies/maggies))
    endfor
    
; now loop through each galaxy and compute zmin and zmax; force the
; apparent magnitude from the model at the redshift of the galaxy to
; be equal to the measured magnitude in the selection band
    offset = fltarr(ngal)
    for jj = 0L, ngal-1 do begin
       appmag = mag[jj] + (dm_ref-dm_gal[jj]) + (kcorr_ref[*,jj]-kcorr_gal[jj]) - $
         (evol_ref-evol_gal[jj])
       offset[jj] = mag[jj] - interpol(appmag,zref,redshift[jj])
;      splog, redshift[jj], mag[jj], offset[jj]
       appmag = appmag + offset[jj]
       result[jj].zmin = interpol(zref,appmag,bright)>zrefmin
       result[jj].zmax = interpol(zref,appmag,faint)<zrefmax
; debugging plot
       if keyword_set(debug) then begin
          djs_plot, zref, appmag, psym=-6, sym=0.5, xsty=3, ysty=3, $
            xrange=[zrefmin,zrefmax], yrange=[bright-0.05,faint+0.05], $
            xtitle='Redshift', ytitle='Apparent AB magnitude'
          djs_oplot, redshift[jj]*[1,1], !y.crange, color='red'
          djs_oplot, !x.crange, mag[jj]*[1,1], color='red'
          djs_oplot, !x.crange, bright*[1,1], color='blue', line=3
          djs_oplot, !x.crange, faint*[1,1], color='blue', line=3
          djs_oplot, result[jj].zmin*[1,1], !y.crange, color='orange', line=5
          djs_oplot, result[jj].zmax*[1,1], !y.crange, color='orange', line=5
          cc = get_kbrd(1)
       endif
    endfor
    
return, result
end
    
