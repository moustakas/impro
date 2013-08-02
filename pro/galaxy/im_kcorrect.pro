;+
; NAME:
;   IM_KCORRECT()
;
; PURPOSE:
;   Wrapper on Blanton's KCORRECT to compute K-corrections from
;   and to an arbitrary set of bandpasses; this is a (shamelessly
;   hacked) generalized version of Blanton's DEEP_KCORRECT.
;
; INPUTS: 
;   redshift       - redshifts [NGAL]
;   maggies        - maggies, Galactic-reddening corrected
;                    [IN_NBAND,NGAL] 
;   ivarmaggies    - inverse variances corresponding to MAGGIES
;                    [IN_NBAND,NGAL]  
;   in_filterlist  - list of filter functions corresponding to MAGGIES
;                    [IN_NBAND] 
;   out_filterlist - list of output filters [OUT_NBAND]
;
; OPTIONAL INPUTS:
;   use_coeffs - use *these* coefficients (do not recompute
;     K-corrections and therefore COEFFS); 
;   band_shift - blueshift of output bandpasses (to get ^{z}b type
;     bands) (default 0.0)  
;   vname - name of fit to use (defaults to 'default')  
;   h100 - Hubble constant divided by 100; the stellar mass estimates
;     and absolute magnitudes are scaled to this value (default 0.7) 
;   omega0 - Omega_matter (default 0.3) 
;   omegal0 - Omega_Lambda (default 0.7)
;   extra - additional keywords for KCORRECT
;
; KEYWORD PARAMETERS: 
;   vega - convert the output magnitudes to Vega (default is AB)   
;   not_closest - disable the default, which is to compute the
;     K-correction using the closest input bandpass to the desired
;     output bandpass  
;   silent - do not print messages to STDOUT
;
; OUTPUTS: 
;   kcorrect - K-corrections satisfying based on the best fit sum of
;     templates: m_R = M_Q + DM(z) + K_QR(z) (all in AB)
;     [OUT_NBAND,NGAL]  
;
; OPTIONAL OUTPUTS:
;   chi2 - chi^2 of fit
;   coeffs - coefficients of fit [5,NGAL]
;   mass - total current stellar mass from model [NGAL]
;   mtol - mass-to-light ratios from model in each output band
;     [OUT_NBAND,NGAL]  
;   obands - which input bands the K-corrections refer to
;     [IN_NBAND,NGAL]  
;   bestmaggies - synthesized *observed-frame* maggies corresponding
;     to IN_FILTERLIST (equivalent to RMAGGIES returned by K-correct)
;     [IN_NBAND,NGAL]  
;   synth_outmaggies_obs - synthesized *observed-frame* maggies
;     corresponding to OUT_FILTERLIST [OUT_NBAND,NGAL]   
;   synth_outmaggies_rest - synthesized *rest-frame* maggies
;     corresponding to OUT_FILTERLIST [OUT_NBAND,NGAL]    
;   absmag - absolute magnitude (for missing data the code uses the
;     synthesized rest-frame photometry in each output band
;     [OUT_NBAND,NGAL]   
;   ivarabsmag - inverse variance of absolute magnitude (for missing
;     data = 0) in each output band [OUT_NBAND,NGAL]   
;   synth_absmag - absolute magnitudes synthesized from the
;     best-fitting model [OUT_NBAND,NGAL]   
;   clineflux - continuum flux at the wavelengths of strong emission
;     lines: [OII], H-beta, [OIII], and H-alpha [erg/s/cm^2/A]
;   uvflux - continuum flux at 1500 and 2800 A from the best-fitting
;     model, useful for computing SFRs [erg/s/cm^2/A]
;
; COMMENTS:
;   This routine uses the bandpass that is closest in effective
;   wavelength to the desired output bandpass to compute the
;   K-correction (identical to always using /CLOSEST in
;   Blanton's DEEP_KCORRECT).
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Aug 07, NYU - shamelessly hacked version
;     of M. Blanton's DEEP_KCORRECT
;   jm08dec21nyu - added NOT_CLOSEST keyword
;   jm09apr30nyu - renamed RMAGGIES to SYNTH_INMAGGIES and added
;     SYNTH_OUTMAGGIES as an optional output; bug fix having to do
;     with OUT_FILTERLIST
;   jm09aug18ucsd - lots of tweaks
;   jm09nov24ucsd - added UVFLUX output
;   jm11aug29ucsd - added SYNTH_ABSMAG optional output 
;
; Copyright (C) 2008-2009, 2011, John Moustakas
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

function im_kcorrect, redshift, maggies, ivarmaggies, in_filterlist, $
  new_out_filterlist, use_coeffs=use_coeffs, band_shift=in_band_shift, vname=in_vname, $
  h100=h100, omega0=omega0, omegal0=omegal0, chi2=chi2, coeffs=coeffs, mass=mass, $
  mtol=mtol, intsfh=intsfh, obands=obands, bestmaggies=bestmaggies, $
  synth_outmaggies_obs=synth_outmaggies_obs, synth_outmaggies_rest=synth_outmaggies_rest, $
  absmag=absmag, ivarabsmag=ivarabsmag, synth_absmag=synth_absmag, clineflux=clineflux, $
  uvflux=uvflux, vega=vega, not_closest=not_closest, silent=silent, reset_rmatrix=reset_rmatrix, $
  psfile=psfile, _extra=extra

    common com_im_kcorrect, out_rmatrix, out_zvals, out_filterlist, band_shift, vname

    nredshift = n_elements(redshift)
    in_nband = n_elements(in_filterlist)
    out_nband = n_elements(new_out_filterlist)
    if (nredshift eq 0L) or (in_nband eq 0L) or (out_nband eq 0L) then begin
       doc_library, 'im_kcorrect'
       return, -1
    endif

    if (n_elements(maggies) eq 0L) or $
      (n_elements(ivarmaggies) eq 0L) then begin
       splog, 'MAGGIES and IVARMAGGIES input required'
       return, -1
    endif
    
    if (n_elements(omega0) eq 0L) then omega0 = 0.3
    if (n_elements(omegal0) eq 0L) then omegal0 = 0.7
    if (n_elements(h100) eq 0L) then h100 = 0.7
    
; need to reset rmatrix if vname or band_shift change or if USE_COEFFS
; was passed
    if (n_elements(use_coeffs) ne 0) then rmatrix = 0
    
    if (n_elements(in_vname) gt 0) then begin
       use_vname=in_vname
    endif else begin
       use_vname='default.nolines'
    endelse
    if (n_elements(vname) gt 0) then begin
       if(vname ne use_vname) then begin
          rmatrix=0
          zvals=0
          out_rmatrix=0
          out_zvals=0
       endif
    endif
    vname = use_vname

    if (n_elements(in_band_shift) gt 0) then $
      use_band_shift = in_band_shift else use_band_shift = 0.0
    if (n_elements(band_shift) gt 0) then begin
       if(band_shift ne use_band_shift) then begin
          out_rmatrix = 0
          out_zvals = 0
       endif
    endif
    band_shift = use_band_shift

; need to reset rmatrix if out_filterlist changes
    if (n_elements(out_filterlist) gt 0) then begin
       if (n_elements(out_filterlist) eq n_elements(new_out_filterlist)) then begin
          for i = 0L, n_elements(out_filterlist)-1L do begin
             if (strtrim(new_out_filterlist[i],2) ne strtrim(out_filterlist[i],2)) then begin
                out_rmatrix = 0
                out_zvals = 0
             endif
          endfor
       endif else begin
          out_rmatrix = 0
          out_zvals = 0
       endelse
    endif 
    out_filterlist = new_out_filterlist

    if keyword_set(reset_rmatrix) then begin
       out_rmatrix = 0
       out_zvals = 0
    endif

; call kcorrect unless USE_COEFFS was passed
    if (n_elements(use_coeffs) eq 0) then begin
       kcorrect, maggies, ivarmaggies, redshift, kcdum, band_shift=band_shift, $
         coeffs=coeffs, rmaggies=bestmaggies, vname=vname, mass=mass, $
         mtol=mtol, absmag=absmag, amivar=ivarabsmag, filterlist=in_filterlist, $
         silent=silent, chi2=chi2, _extra=extra;, rmatrix=im_rmatrix, $
;        zvals=im_zvals
    endif else begin
       dim = size(use_coeffs,/dim)
       if (dim[0] ne 5) or (dim[1] ne nredshift) then $
         message, 'USE_COEFFS must be a [5,NGAL] array!'
       coeffs = use_coeffs

       k_load_vmatrix, vm, lam, vfile=vfile, lfile=lfile, vpath=vpath, vname=vname
       k_projection_table, in_rmatrix, vm, lam, in_zvals, $
         in_filterlist, /silent
       k_reconstruct_maggies, coeffs, redshift, bestmaggies, $
         rmatrix=in_rmatrix, zvals=in_zvals, /silent
    endelse

; calculate the preliminaries; force kcorrect to recalculate RMATRIX
; and ZVALS every time
    if (keyword_set(out_rmatrix) eq 0) or (keyword_set(out_zvals) eq 0) then begin
       if (keyword_set(vmatrix) eq 0) or (keyword_set(lambda) eq 0) then begin
          k_load_vmatrix, vmatrix, lambda, vfile=vfile, $
            lfile=lfile, vpath=vpath, vname=vname
          k_projection_table, out_rmatrix, vmatrix, lambda, out_zvals, $
            out_filterlist, _extra=extra
       endif
    endif

; synthesized observed- and rest-frame output maggies
    k_reconstruct_maggies, coeffs, replicate(band_shift,nredshift), $
      synth_outmaggies_rest, rmatrix=out_rmatrix, zvals=out_zvals, $
      silent=silent
    synth_outmaggies_rest = synth_outmaggies_rest/(1.0+band_shift)
    k_reconstruct_maggies, coeffs, redshift, synth_outmaggies_obs, $
      rmatrix=out_rmatrix, zvals=out_zvals, silent=silent

    obands = lindgen(n_elements(out_filterlist))#replicate(1L,nredshift)

; default: use the magnitude in the closest filter bandpass with
; *measured* photometry (but upper limits are OK)
    if (keyword_set(not_closest) eq 0) then begin
       lambda_in = k_lambda_eff(filterlist=in_filterlist)
       lambda_out = k_lambda_eff(filterlist=out_filterlist,band_shift=band_shift)
       for ii = 0L, nredshift-1 do begin
          for jj = 0L, n_elements(lambda_out)-1 do begin
             lambdadist = abs(lambda_in/(1.0+redshift[ii])-lambda_out[jj])
             dmin = min(lambdadist + (ivarmaggies[*,ii] eq 0)*1D6,imin)
;            dmin = min(lambdadist,imin) ; old code
             obands[jj,ii] = imin
          endfor
       endfor
    endif

    kcorrect = fltarr(n_elements(out_filterlist), nredshift)
    for i=0L, nredshift-1L do $
      for j=0L, n_elements(out_filterlist)-1L do $
        kcorrect[j,i]=synth_outmaggies_rest[j,i]/bestmaggies[obands[j,i],i]
    kcorrect = +2.5*alog10(kcorrect)

; convert to Vega magnitudes if requested (code stolen from
; SDSS2BESSELL) 
    if keyword_set(vega) then begin
       for i = 0L, n_elements(out_filterlist)-1L do begin
          v2ab = k_vega2ab(filterlist=out_filterlist[i], /kurucz, $
            band_shift=band_shift,silent=silent)
          ;; you need to convert FROM AB TO VEGA. the K-correction is
          ;; M=m-DM-K, so adding v2ab to K makes M=m-DM-K-v2ab,
          ;; equivalent to M=m-DM-K+ab2v
          kcorrect[i,*] = kcorrect[i,*]+v2ab[0]
;;        mtol[i,*]=mtol[i,*]*10.^(-0.4*v2ab[0])
       endfor
    endif

; calculate absolute magnitudes    
    synth_absmag = fltarr(n_elements(out_filterlist), nredshift)
    ivarabsmag = fltarr(n_elements(out_filterlist), nredshift)
    if arg_present(absmag) or arg_present(synth_absmag) then begin
       for i = 0L, n_elements(out_filterlist)-1L do $
         synth_absmag[i,*] = -2.5*alog10(synth_outmaggies_rest[i,*])- $
         lf_distmod(redshift, omega0=omega0, omegal0=omegal0)
       absmag = synth_absmag
       for j = 0L, nredshift-1L do begin
          igood = where(ivarmaggies[obands[*,j],j] gt 0. and $
            maggies[obands[*,j],j] gt 0., ngood)
          if(ngood gt 0) then begin
             absmag[igood,j] = -2.5*alog10(maggies[obands[igood,j],j])- $
               lf_distmod(redshift[j],omega0=omega0,omegal0=omegal0)- $
               kcorrect[igood,j]
             ivarabsmag[igood,j] = maggies[obands[igood,j],j]^2* $
               ivarmaggies[obands[igood,j],j]*(0.4*alog(10.))^2
          endif
       endfor
    endif

; if requested, compute the continuum flux corresponding to the
; wavelengths of strong emission lines
    if arg_present(clineflux) then begin
       k_load_vmatrix, vmatrix, lambda, vfile=vfile, $
         lfile=lfile, vpath=vpath, vname=vname
       linewave = [3727.420,4861.325,4958.911,5006.843,6562.800]
       lineindx = findex(k_lambda_to_centers(lambda),linewave)
       clineflux = fltarr(n_elements(linewave),nredshift)
       for ii = 0L, nredshift-1L do clineflux[*,ii] = $ ; [erg/s/cm2/A]
         interpolate(smooth(reform(vmatrix#coeffs[*,ii]),10),lineindx) 
    endif
       
; if requested, compute the continuum flux at 1500 and 2800 A
    if arg_present(uvflux) then begin
       k_load_vmatrix, vmatrix, lambda, vfile=vfile, $
         lfile=lfile, vpath=vpath, vname=vname
       wave = k_lambda_to_centers(lambda)
       uvwave = [1500.0,2800.0]
       uvflux = fltarr(n_elements(uvwave),nredshift)
;      uvindx = findex(wave,uvwave)
;      for ii = 0L, nredshift-1L do uvflux[*,ii] = $ ; [erg/s/cm2/A]
;        interpolate(reform(vmatrix#coeffs[*,ii]),uvindx)
       ww1500 = where((wave gt uvwave[0]-30) and (wave lt uvwave[0]+30))
       ww2800 = where((wave gt uvwave[1]-30) and (wave lt uvwave[1]+30))
       for ii = 0L, nredshift-1L do begin
          sed = reform(vmatrix#coeffs[*,ii])
          uvflux[0,ii] = djs_mean(sed[ww1500])
          uvflux[1,ii] = djs_mean(sed[ww2800])
       endfor
    endif
       
;   rmatrix = out_rmatrix

; scale everything to h100; note that K-correct *always* uses h100=1
    if arg_present(mass) then mass = mass/h100^2.0
    if arg_present(absmag) or arg_present(synth_absmag) then begin
       absmag = absmag - 5.0*alog10(1.0/h100)
       synth_absmag = synth_absmag - 5.0*alog10(1.0/h100)
    endif

; if requested, build a QAplot    
    if (n_elements(psfile) ne 0) then begin
       str = replicate({z: 0.0, chi2: 0.0, mass: 0.0, coeffs: fltarr(5), $
         maggies: fltarr(in_nband), $
         ivarmaggies: fltarr(in_nband), bestmaggies: fltarr(in_nband), $
         out_filterlist: out_filterlist, synth_outmaggies_obs: fltarr(out_nband)},nredshift)
       str.z = redshift
       str.chi2 = chi2
       str.mass = alog10(mass)
       str.coeffs = coeffs
       str.maggies = maggies
       str.ivarmaggies = ivarmaggies
       str.bestmaggies = bestmaggies
       str.synth_outmaggies_obs = synth_outmaggies_obs
       str.out_filterlist = out_filterlist
       kcorrect_qaplot, str[0:100<(nredshift-1L)], in_filterlist=in_filterlist, $
         vname=vname, psfile=psfile, /clobber
    endif
    
return, kcorrect
end
