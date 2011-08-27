;+
; NAME:
;   IM_CALC_VMAX
;
; PURPOSE:
;   Calculate Vmax, given the flux and redshift limits of the sample. 
;
; INPUTS: 
;   appm - apparent magnitude of each object in the selection band [N] 
;   vmax_coeffs - [NT, N] K-correction coeffs for each object
;   select_filter - selection filter name (e.g., sdss_r0.par)
;   marea  - area of mag limit region [NA]
;   mmin - minimum apparent mag for each mag limit region [NA]
;   mmax - maximum apparent mag for each mag limit region [NA]
;   sample_zmin - minimum redshift of the sample
;   sample_zmax - maximum redshift of the sample
;   actual_z - actual redshift of each object [N]
;
; OPTIONAL INPUTS:
;   q0, q1, qz0 - evolutionary parameters (see K_EVOLVE)
;   omega0 - omega_matter (default: 0.3)
;   omegal0 - omega_lambda (default: 0.7)
;   h100 - Hubble constant divided by 100; Vmax is scaled to this
;     value (default 0.7)
;   im - index of magnitude limit region (default 0) [N]
;   vname - templates used when computing K-corrections (default
;     'default') 
;   magoffset - 
; 
; KEYWORD PARAMETERS:
;   silent - suppress messages to STDOUT
;   smoothhr - 
;
; OPTIONAL OUTPUTS:
;   zmin - local zmin of each [N]
;   zmax - local zmax of each [N]
;   vmax - Vmax of each integrated over all mag limit regions [N]
;
; COMMENTS:
;   VMAX is returned in h^{-3} Mpc^3 comoving
;
; COMMENTS:
;   Shamelessly hacked version of M. Blanton's LF_CALC_VMAX! 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009-Nov-11, UCSD - written
;
; Copyright (C) 2009, John Moustakas
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

function im_calc_vmax, appm, vmax_coeffs, in_select_filter, marea, mmin, mmax,  $
  sample_zmin, sample_zmax, actual_z=actual_z, q0=in_q0, q1=in_q1, $
  qz0=in_qz0, omega0=in_omega0, omegal0=in_omegal0, h100=in_h100, im=im, $
  vname=vname, magoffset=in_magoffset, rmatrix=rmatrix, zvals=zvals, $
  smoothr=smoothr, silent=silent

    nobj = n_elements(appm)
    if (nobj eq 0) or (n_elements(vmax_coeffs) eq 0) or $
      (n_elements(in_select_filter) eq 0) or $
      (n_elements(marea) eq 0) or (n_elements(mmin) eq 0) or $
      (n_elements(mmax) eq 0) or (n_elements(sample_zmin) eq 0) or $
      (n_elements(sample_zmax) eq 0) then begin
       doc_library, 'im_calc_vmax'
       return, -1
    endif
    
; defaults
    if (n_elements(in_omega0) eq 0) then in_omega0 = 0.3
    if (n_elements(in_omegal0) eq 0) then in_omegal0 = 0.7
    if (n_elements(in_h100) eq 0) then in_h100 = 0.7
    if (n_elements(in_q0) eq 0) then in_q0 = 0.0
    if (n_elements(in_q1) eq 0) then in_q1 = 0.0
    if (n_elements(in_qz0) eq 0) then in_qz0 = 0.0
    if (n_elements(in_magoffset) eq 0) then in_magoffset = 0.0
    if (n_elements(im) eq 0) then im = lonarr(nobj)
    if (n_elements(vname) eq 0) then vname = 'default'

    omega0 = in_omega0
    omegal0 = in_omegal0
    h100 = in_h100
    select_filter = in_select_filter
    band_shift = 0.0 ; not needed, but set to zero here
    q0 = in_q0
    q1 = in_q1
    qz0 = in_qz0
    magoffset = in_magoffset
    nv = n_elements(vmax_coeffs)/n_elements(appm)
    nk = 1

; calculate the absolute magnitude in the selection band
    if (n_elements(zvals) ne 0) then this_zvals = zvals
    k_load_vmatrix, this_vmatrix, this_lambda, vname=vname
    k_projection_table, this_rmatrix, this_vmatrix, this_lambda, $
      this_zvals, select_filter, silent=silent
    k_reconstruct_maggies, vmax_coeffs, replicate(band_shift,nobj), recm, $
      rmatrix=this_rmatrix, zvals=this_zvals, filterlist=select_filter, $
      silent=silent
    recm = recm/(1.0+band_shift)
    k_reconstruct_maggies, vmax_coeffs, actual_z, recm0, rmatrix=this_rmatrix, $
      zvals=this_zvals, filterlist=select_filter, silent=silent
    select_kcorr = reform(2.5*alog10(recm/recm0))
    select_absmag = appm-lf_distmod(actual_z,omega0=omega0,omegal0=omegal0)-select_kcorr ; h=1

; run kcorrect to get rmatrix and zvals
    if (n_elements(rmatrix) eq 0) or (n_elements(zvals) eq 0) then $
      kcorrect, dummaggies, dummaggies_ivar, 0.0, dumk, band_shift=band_shift, $
      filterlist=[select_filter], rmatrix=rmatrix, zvals=zvals, coeffs=vmax_coeffs[*,0],$
      vname=vname, silent=silent
    nz = n_elements(rmatrix)/nv/nk
    
    if (keyword_set(smoothr)) then begin
       old_rmatrix=rmatrix
       for i=0L,  nv-1L do $
         rmatrix[*,i]=smooth(old_rmatrix[*,i], 4., /edge)
    endif

; for each object...
    vmax=fltarr(nobj)
    zmin=fltarr(nobj)
    zmax=fltarr(nobj)
    for i=0L, nobj-1L do begin
       if (keyword_set(silent) eq 0) then if ((i mod 1000) eq 0) then $
         splog,' galaxy '+string(i)
       curr_select_absmag=k_evolve(select_absmag[i], actual_z[i], q0, q1, qz0)
       curr_coeffs=vmax_coeffs[*,i]
       if (keyword_set(smoothr)) then begin
          off=-2.5*alog10(interpol(rmatrix#abs(curr_coeffs), zvals, $
            actual_z[i]))+ $
            2.5*alog10(interpol(old_rmatrix#abs(curr_coeffs), zvals, $
            actual_z[i]))
          curr_select_absmag=curr_select_absmag-off
       endif
       
       for j=0L, n_elements(marea)-1L do begin
          curr_zmin=-1.0
          curr_zmax=-1.0
          soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
            root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')
          curr_mmin=mmin[j]
          curr_mmax=mmax[j]
          retval=call_external(soname, 'idl_lf_calc_vmax', float(curr_select_absmag), $
            float(curr_coeffs),long(nv),float(zvals), $
            long(nz),float(rmatrix),long(nk), $
            float(sample_zmin),float(sample_zmax), $
            float(curr_mmin),float(curr_mmax), $
            float(q0),float(q1),float(qz0), $
            float(band_shift),float(magoffset), $
            float(omega0),float(omegal0), $
            float(curr_zmin),float(curr_zmax))
          
; calculate vmax; don't allow it to be less than zero
          if (im[i] eq j) then begin
             zmin[i]=curr_zmin
             zmax[i]=curr_zmax
             if (zmin[i] gt actual_z[i]) then begin
                if (abs(curr_mmin-appm[i]) gt 0.005 AND $
                  abs(zmin[i]-actual_z[i]) gt 1.e-4) then begin
                   splog, i, actual_z[i], zmin[i], zmax[i], appm[i], $
                     select_absmag[i] + lf_distmod(actual_z[i]) + select_kcorr[i], $
                     curr_mmin, curr_mmax, curr_select_absmag
                   stop
                endif
             endif
             if (zmax[i] lt actual_z[i]) then begin
                if (abs(curr_mmax-appm[i]) gt 0.005 AND $
                  abs(zmax[i]-actual_z[i]) gt 1.e-4)then begin
;                  message,'galaxy outside allowed putative selection limits', /info
                   splog, i, actual_z[i], zmin[i], zmax[i], appm[i], $
                     select_absmag[i] + lf_distmod(actual_z[i]) + select_kcorr[i], $
                     curr_mmin, curr_mmax, curr_select_absmag
                endif
             endif
          endif
          vmax[i] = vmax[i]+((marea[j]*(lf_comvol(curr_zmax,omega0=omega0,omegal0=omegal0)- $
            lf_comvol(curr_zmin,omega0=omega0,omegal0=omegal0))) > 0.0)/3.0
       endfor
    endfor 

    vol = (marea/3.0)*(lf_comvol(actual_z)-(lf_comvol(sample_zmin))[0]) ; comoving volume

; scale to h100; note that K-correct *always* uses h100=1;
; don't scale SELECT_ABSMAG (which is for h100=1), which has no
; use outside of this routine
    vmax = vmax/h100^3.0
    vol = vol/h100^3.0
    
; pack everything into a structure and return
    out = {$
      omega0:               omega0,$
      omegal0:             omegal0,$
      h100:                   h100,$
      q0:                       q0,$
      q1:                       q1,$
      qz0:                     qz0,$
      vname:                 vname,$
      marea:                 marea,$ 
      mmin:                   mmin,$
      mmax:                   mmax,$
      sample_zmin:     sample_zmin,$
      sample_zmax:     sample_zmax,$
      magoffset:         magoffset,$
      select_filter: select_filter,$
      select_absmag: select_absmag,$ ; note! h100=1
      select_kcorr:   select_kcorr,$
      vmax_coeffs:     vmax_coeffs,$ 
      appm:                   appm,$
      actual_z:           actual_z,$
      im:                       im,$
      zmin:                   zmin,$
      zmax:                   zmax,$
      vmax:                   vmax,$ ; 
      vol:                     vol}  ; 
    
return, out
end
