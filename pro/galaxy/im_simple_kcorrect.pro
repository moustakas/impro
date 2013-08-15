;+
; NAME:
;   IM_SIMPLE_KCORRECT()
;
; PURPOSE:
;   Compute K-corrections from and to an arbitrary set of bandpasses,
;   but using a given SED for each object (see also IM_KCORRECT).
;
; INPUTS: 
;   redshift      - redshifts [NGAL]
;   maggies       - maggies, Galactic-reddening corrected
;                   [IN_NBAND,NGAL] 
;   ivarmaggies   - inverse variances corresponding to MAGGIES
;                   [IN_NBAND,NGAL]  
;   in_filterlist - list of filter functions corresponding to MAGGIES
;                   [IN_NBAND] 
;   sedwave       - rest-frame wavelength vector
;   sedflux       - rest-frame spectral energy distribution
;
; OPTIONAL INPUTS: 
;   band_shift - blueshift of output bandpasses (to get ^{z}b type
;                bands) (default 0.0) 
;   vname      - name of fit to use (defaults to 'default.nolines') 
;   omega0     - Omega_matter (default 0.3)
;   omegal0    - Omega_Lambda (default 0.7)
;   h100       - Hubble constant divided by 100; the stellar mass
;                estimates and absolute magnitudes are scaled to this
;                value (default 0.7)
;
; KEYWORD PARAMETERS: 
;   sdss   - see OUT_FILTERLIST
;   vega   - convert the output magnitudes to Vega (default is AB)  
;   silent - do not print messages to STDOUT
;
; OUTPUTS: 
;   kcorrect   - K-corrections satisfying based on the best fit sum of
;                templates: m_R = M_Q + DM(z) + K_QR(z) (all in AB)
;                [OUT_NBAND,NGAL] 
;
; OPTIONAL OUTPUTS:
;   out_filterlist - list of output filters; default UBVRI (Bessell),
;                    unless /SDSS, in which case the 'ugriz' filters
;                    are used [OUT_NBAND]
;   chi2       - chi^2 of fit
;   rmaggies   - reconstructed maggies from the fit [IN_NBAND,NGAL] 
;   obands     - which input bands the K-corrections refer to
;                [IN_NBAND,NGAL] 
;   scale      - derived normalization factor used to scale the input
;                SED to the observed photometry; note that if the
;                input SED is normalized per unit mass, then SCALE is
;                trivially related to the stellar mass 
;   absmag     - absolute magnitude (for missing data, substitutes 
;                model fit) in each output band [OUT_NBAND,NGAL] 
;   ivarabsmag - inverse variance of absolute magnitude (for
;                missing data = 0) in each output band
;                [OUT_NBAND,NGAL]  
;   synth_absmag - absolute magnitudes synthesized from the
;     best-fitting model [OUT_NBAND,NGAL]   
;   clineflux  - continuum flux at the wavelengths of strong emission 
;                lines: [OII], H-beta, [OIII], and H-alpha [erg/s/cm2/A] 
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
;   J. Moustakas, 2009 Jan 29, NYU - simplified version of
;     IM_KCORRECT() 
;   jm09aug18ucsd - lots of tweaks
;   jm10jan07ucsd - added UVFLUX output
;   jm11aug29ucsd - added SYNTH_ABSMAG optional output 
;
; Copyright (C) 2009-2011, John Moustakas
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

function im_simple_kcorrect, redshift, maggies, ivarmaggies, in_filterlist, $
  out_filterlist, sedwave1, sedflux1, band_shift=band_shift, h100=h100, $
  omega0=omega0, omegal0=omegal0, chi2=chi2, scale=scale, obands=obands, $
  bestmaggies=bestmaggies, synth_outmaggies_obs=synth_outmaggies_obs, $
  synth_outmaggies_rest=synth_outmaggies_rest, absmag=absmag, $
  ivarabsmag=ivarabsmag, synth_absmag=synth_absmag, clineflux=clineflux, $
  uvflux=uvflux, vega=vega, silent=silent

    nredshift = n_elements(redshift)
    in_nband = n_elements(in_filterlist)
    out_nband = n_elements(out_filterlist)
    if (nredshift eq 0L) or (in_nband eq 0L) or (out_nband eq 0L) then begin
       doc_library, 'im_simple_kcorrect'
       return, -1
    endif

    if (n_elements(maggies) eq 0L) or $
      (n_elements(ivarmaggies) eq 0L) then begin
       splog, 'MAGGIES and IVARMAGGIES input required'
       return, -1
    endif
    
; handle the input SED    
    ndim_wave = size(sedwave1,/n_dim)
    ndim_flux = size(sedflux1,/n_dim)
    dims_wave = size(sedwave1,/dim)
    dims_flux = size(sedflux1,/dim)

    if (dims_wave[0] eq 0L) or (dims_flux[0] eq 0L) then begin
       splog, 'SEDWAVE and SEDFLUX input required'
       return, -1
    endif
    if (dims_wave[0] ne dims_flux[0]) then begin
       splog, 'Dimensions of SEDWAVE and SEDFLUX must match'
       return, -1
    endif
    if (ndim_wave eq 2L) then if (dims_wave[1] ne nredshift) then begin
       splog, 'Dimensions of REDSHIFT and SEDWAVE must match'
       return, -1
    endif
    if (ndim_flux eq 2L) then if (dims_flux[1] ne nredshift) then begin
       splog, 'Dimensions of REDSHIFT and SEDFLUX must match'
       return, -1
    endif

    if (ndim_wave eq 2L) then sedwave = sedwave1 else $
      sedwave = rebin(sedwave1,dims_wave[0],nredshift)
    if (ndim_flux eq 2L) then sedflux = sedflux1 else $
      sedflux = rebin(sedflux1,dims_flux[0],nredshift)

; other preliminaries     
    if (n_elements(band_shift) eq 0L) then band_shift = 0.0
    if (n_elements(omega0) eq 0L) then omega0 = 0.3
    if (n_elements(omegal0) eq 0L) then omegal0 = 0.7
    if (n_elements(h100) eq 0L) then h100 = 0.7

; compute the relevant observed- and rest-frame magnitudes    
    chi2 = fltarr(nredshift)
    scale = fltarr(nredshift)+1.0
    bestmaggies = fltarr(n_elements(in_filterlist),nredshift)
    synth_outmaggies_obs = fltarr(out_nband,nredshift)
    synth_outmaggies_rest = fltarr(out_nband,nredshift)

    for jj = 0L, nredshift-1L do begin
       restlambda = reform(k_lambda_to_edges(sedwave[*,jj]))
       restflux = reform(sedflux[*,jj])
       
       lambda = reform(k_lambda_to_edges(sedwave[*,jj]*(1.0+redshift[jj])))
       flux = reform(sedflux[*,jj]/(1.0+redshift[jj]))
       
       bestmaggies[*,jj] = k_project_filters(lambda,flux,$          ; in_filterlist, observed
         silent=silent,filterlist=in_filterlist)
       synth_outmaggies_rest[*,jj] = k_project_filters(restlambda,$ ; out_filterlist, rest
         restflux,silent=silent,filterlist=out_filterlist,$
         band_shift=band_shift)
       synth_outmaggies_obs[*,jj] = k_project_filters(lambda,flux,$ ; out_filterlist, observed
         silent=silent,filterlist=out_filterlist)

; note the absolute value, to protect against negative fluxes       
       scale[jj] = total(ivarmaggies[*,jj]*abs(maggies[*,jj])*bestmaggies[*,jj],/double)/$
         total(ivarmaggies[*,jj]*bestmaggies[*,jj]^2.0,/double)
       chi2[jj] = total(ivarmaggies[*,jj]*(abs(maggies[*,jj])-scale[jj]*bestmaggies[*,jj])^2.0,/double)
       scale[jj] = abs(scale[jj])
       
       bestmaggies[*,jj] = scale[jj]*bestmaggies[*,jj]
       synth_outmaggies_rest[*,jj] = scale[jj]*synth_outmaggies_rest[*,jj]
       synth_outmaggies_obs[*,jj] = scale[jj]*synth_outmaggies_obs[*,jj]
    endfor

; choose the closest filter bandpass
    obands = lindgen(out_nband)#replicate(1L,nredshift)
    lambda_in = k_lambda_eff(filterlist=in_filterlist)
    lambda_out = k_lambda_eff(filterlist=out_filterlist,band_shift=band_shift)
    for ii = 0L, nredshift-1 do begin
       for jj = 0L, n_elements(lambda_out)-1 do begin
          lambdadist = abs(lambda_in/(1.0+redshift[ii])-lambda_out[jj])
          dmin = min(lambdadist + (ivarmaggies[*,ii] eq 0)*1D6,imin)
;         dmin = min(lambdadist,imin) ; old code
          obands[jj,ii]= imin
       endfor
    endfor

; compute the K-correction
    kcorrect = fltarr(out_nband,nredshift)
    for i = 0L, nredshift-1L do $
      for j = 0L, out_nband-1L do $
        kcorrect[j,i]=synth_outmaggies_rest[j,i]/bestmaggies[obands[j,i],i]
    kcorrect = +2.5*alog10(kcorrect)
; convert to Vega magnitudes if requested (code stolen from
; SDSS2BESSELL) 
    if keyword_set(vega) then begin
       for i = 0L, out_nband-1L do begin
          v2ab = k_vega2ab(filterlist=out_filterlist[i], /kurucz, $
            band_shift=band_shift,silent=silent)
; you need to convert FROM AB TO VEGA. the K-correction is M=m-DM-K,
; so adding v2ab to K makes M=m-DM-K-v2ab, equivalent to M=m-DM-K+ab2v
          kcorrect[i,*] = kcorrect[i,*]+v2ab[0]
       endfor
    endif

; calculate absolute magnitudes    
    synth_absmag = fltarr(out_nband,nredshift)
    ivarabsmag = fltarr(out_nband,nredshift)
    if arg_present(absmag) or arg_present(synth_absmag) then begin
       for i = 0L, out_nband-1L do $
         synth_absmag[i,*]=-2.5*alog10(synth_outmaggies_rest[i,*])- $
         lf_distmod(redshift, omega0=omega0, omegal0=omegal0)
       absmag = synth_absmag
       for j=0L, nredshift-1L do begin
          igood = where(ivarmaggies[obands[*,j],j] gt 0.0 and $
            maggies[obands[*,j],j] gt 0.0, ngood)
          if (ngood gt 0) then begin
             absmag[igood,j] = -2.5*alog10(maggies[obands[igood,j],j])- $
               lf_distmod(redshift[j],omega0=omega0,omegal0=omegal0)- $
               kcorrect[igood,j]
             ivarabsmag[igood,j]=maggies[obands[igood,j],j]^2* $
               ivarmaggies[obands[igood,j],j]*(0.4*alog(10.))^2
          endif
       endfor 
    endif

; if requested, compute the continuum flux corresponding to the
; wavelengths of strong emission lines
    if arg_present(clineflux) then begin
       linewave = [3727.420,4861.325,4958.911,5006.843,6562.800]
       clineflux = fltarr(n_elements(linewave),nredshift)
       for ii = 0L, nredshift-1L do clineflux[*,ii] = $ ; [erg/s/cm2/A]
         interpol(scale[ii]*sedflux[*,ii],sedwave[*,ii],linewave)
    endif
       
; if requested, compute the continuum flux at 1500 and 2800 A
    if arg_present(uvflux) then begin
       uvwave = [1500.0,2800.0]
       uvflux = fltarr(n_elements(uvwave),nredshift)
       for ii = 0L, nredshift-1L do uvflux[*,ii] = $ ; [erg/s/cm2/A]
         interpol(scale[ii]*sedflux[*,ii],sedwave[*,ii],uvwave)
    endif
       
; scale everything to h100; note that K-correct *always* uses h100=1
    if arg_present(absmag) or arg_present(synth_absmag) then begin
       absmag = absmag - 5.0*alog10(1.0/h100)
       synth_absmag = synth_absmag - 5.0*alog10(1.0/h100)
    endif

return, kcorrect
end
