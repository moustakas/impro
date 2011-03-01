;+
; THIS CODE IS RIDICULOUSLY SLOW, BUT CONTAINS SOME USEFUL
; CONVENTIONS!  USE IM_KCORRECT INSTEAD!
;
; NAME:
;       IM_PHOTO_KCORR()
;
; PURPOSE:
;       This is a wrapper on M. Blanton's k-correct to compute
;       k-corrections and many other useful quantities.
;
; CALLING SEQUENCE:
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED: 
;
; COMMENTS:
;      From Blanton in response to "How can I overplot the observed
;      fluxes on the best-fitting template?" 
; 
;      "I claim in the documentation that the units are erg/S/cm2/A  
; 
;         fluxfactor= 3.631D-29*(3D18)/(central)^2 ; AB nMgy to f_lambda
;         flambda= nmgys*fluxfactor
; 
;      Takes you from nanomaggies to f_lambda. 
; 
;      Then remember that the vmatrix is in the restframe, so you have
;      to multiply the wavelengths by (1+z) and divide the spectrum by
;      (1+z) to compare to the observed fluxes.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Oct 26, NYU - written, based on lots of
;         older code
;
; Copyright (C) 2006, John Moustakas
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

function im_photo_kcorr, kcorrphot, obsfilters, restfilters=restfilters, restvega=restvega, $
  restbandpasses=restbandpasses, filters_path=filters_path, kcorr_method=kcorr_method, $
  debug=debug

    ngalaxy = n_elements(kcorrphot)
    nobsfilter = n_elements(obsfilters)

    if (ngalaxy eq 0L) or (nobsfilter eq 0L) then begin
       print, 'Syntax - result = im_photo_kcorr(kcorrphot,obsfilters,restfilters=,$'
       print, '   restvega=,restbandpasses=,filters_path=,kcorr_method=,/debug)'
       return, -1L
    endif

; define the cosmology and some constants

    red, omega0=0.3, omegalambda=0.7, h100=0.70
    omega0 = redomega0() & omegal = redomegal() & h100 = redh100()

    imf_chabrier_to_salpeter = +0.25
    light = 2.99792458D18       ; speed of light [A/s]

    if keyword_set(debug) then begin
       plotsym, 0, 2.0, /fill, color=djs_icolor('red')
       im_window, 0, xratio=0.5, yratio=0.4
    endif
    
    if (n_elements(filters_path) eq 0L) then filters_path = getenv('KCORRECT_DIR')+'/data/filters/'
    if (n_elements(kcorr_method) eq 0L) then kcorr_method = 2L

; load the rest filters of interest

    if (arg_present(restfilters) eq 0L) and (arg_present(restvega) or arg_present(restbandpasses)) then begin
       splog, 'RESTFILTERS, RESTVEGA, and RESTBANDPASSES must be defined together.'
       return, -1L
    endif
    
    splog, 'Loading rest-frame filters:'
    if (n_elements(restfilters) eq 0L) then begin
       restfilters = ['bessell_'+['U','B','V','R','I'],'sdss_'+['u0','g0','r0'],'twomass_'+['J','H']]+'.par'
       restvega = [1,1,1,1,1,0,0,0,1,1]
       restbandpasses = ['U','B','V','R','I','sdss_u','sdss_g','sdss_r','twomass_j','twomass_h']
    endif
    nrestfilter = n_elements(restfilters)
    restfiltinfo = im_filterspecs(filterlist=restfilters,/verbose)

    notvega = where((restvega eq 0L),nnotvega)
    if (nnotvega ne 0L) then restfiltinfo[notvega].vega2ab = 0.0

; observed-frame photometric filter information

    splog, 'Reading observed-frame filters:'
    obsbandpasses = repstr(obsfilters,'.par','')
    obsfiltinfo = im_filterspecs(filterlist=obsfilters,filterpath=filters_path,/verbose)
    
; initialize the output data structure

    result = {$

      galaxy:                                     '', $
      z:                                      -999.0, $
      kcorr_chi2:                             -999.0, $ ; reduced chi2
      kcorr_mass:                             -999.0, $
      kcorr_mets:                             -999.0, $
      kcorr_intsfh:                           -999.0, $
      kcorr_zzmax:                            -999.0, $ ; compute z/zmax (see below)
      kcorr_kcorr:         fltarr(nrestfilter)-999.0, $ ; k-correction
      kcorr_filter:     replicate('...',nrestfilter), $ ; reference filter
      kcorr_mrest:         fltarr(nrestfilter)-999.0, $ ; k-corrected, rest frame, using "actual observed photometry"
      kcorr_mrest_err:     fltarr(nrestfilter)-999.0, $ ; error in MREST
      kcorr_mrest_synth:         fltarr(nrestfilter), $ ; synthesized, rest frame
;     kcorr_mrest_combined:      fltarr(nrestfilter), $ ; rest frame, using combined "actual" and "synthesized" observed photometry
      kcorr_mobs_ab:        fltarr(nobsfilter)-999.0, $ ; "actual observed photometry"
      kcorr_mobs_ab_err:    fltarr(nobsfilter)-999.0, $ ; error in MOBS_AB
      kcorr_mobs_ab_synth:        fltarr(nobsfilter), $ ; "synthesized observed photometry"
      kcorr_coeffs:                        fltarr(5)}   ; eigentemplate coefficients

    result = replicate(result,ngalaxy)

    result.galaxy            = kcorrphot.galaxy
    result.z                 = kcorrphot.z
    result.kcorr_mobs_ab     = kcorrphot.mobs_ab
    result.kcorr_mobs_ab_err = kcorrphot.mobs_ab_err

; compute the k-corrections      

    splog, 'Fitting k-correct templates.'
    t0 = systime(1)
    kcorrect, kcorrphot.maggies, kcorrphot.maggies_ivar, kcorrphot.z, $
      kcorrect, band_shift=band_shift, filterlist=obsfilters, $
      filterpath=filters_path, rmatrix=rmatrix, zvals=zvals, $
      lambda=lambda, vmatrix=vmatrix, coeffs=coeffs, chi2=chi2, $
      maxiter=maxiter, omega0=omega0, omegal0=omegal, mass=mass, $
      intsfh=intsfh, mets=mets, b300=b300, b1000=b1000, mtol=mtol, $
      vname='default.nolines' ; NOTE!
    splog, 'Total time = '+string((systime(1)-t0)/60.0,format='(G0.0)')+' minutes.'

; store the basic results; convert the stellar mass to my IMF and
; Hubble constant; store the reduced chi2 below

    good = where((mass gt 1E6) and finite(mass),ngood)
    if (ngood ne 0L) then begin

       result[good].kcorr_mass = alog10(mass[good]) - alog10(h100) + imf_chabrier_to_salpeter
       result[good].kcorr_mets = mets[good]
       result[good].kcorr_intsfh = alog10(intsfh[good])
       result[good].kcorr_coeffs = coeffs[*,good]

    endif

; compute rest-frame magnitudes    
    
    k_reconstruct_maggies, coeffs, kcorrphot.z*0.0, reconstruct_maggies, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, vmatrix=vmatrix, $
      lambda=lambda, filterpath=filters_path, filterlist=restfilters, $
      zmin=zmin, zmax=zmax, nz=nz, silent=silent

stop    
    
; for each object, reconstruct the best-fitting SED and compute
; k-corrections and rest-frame apparent magnitudes

    restwave = lambda
       
    splog, 'Computing k-corrections.'
    t0 = systime(1)
    for i = 0L, ngalaxy-1L do begin
;   for i = 0L, 45L do begin

       print, format='("Galaxy ",I0,"/",I0,".    ",A1,$)', i+1, ngalaxy, string(13b)

       if (result[i].kcorr_mass gt -999.0) and (result[i].z gt 0.0) then begin
   
          good = where((result[i].kcorr_mobs_ab gt 0.0),ngood)
          if (ngood ne 0L) then result[i].kcorr_chi2 = chi2[i] / float(ngood)

          restwave = k_lambda_to_centers(lambda)
          wave = restwave*(1.0+kcorrphot[i].z)
          restflux = vmatrix#double(coeffs[*,i]) ; f_lambda
          flux = restflux/(1.0+kcorrphot[i].z)
          flux_fnu = flux*wave^2.0/light*10.0^(0.4*48.6) ; f_nu
          flux_ab = -2.5*alog10(flux_fnu>1D-50)   ; AB
          
; synthesize magnitudes in the observed frame

          obsmaggies_synth = im_filtermag(wave,flux,filtinfo=obsfiltinfo)
          result[i].kcorr_mobs_ab_synth = -2.5*alog10(obsmaggies_synth) ; AB

; compute synthesized and k-corrected rest-frame apparent magnitudes;
; remember: m_R = m_Q + K_QR; for each output _rest_ bandpass,
; identify the _observed_ bandpass that minimizes the k-correction

          restmaggies_synth = im_filtermag(restwave,restflux,filtinfo=restfiltinfo)
;         restmaggies_synth2 = reform(k_project_filters(k_lambda_to_edges(restwave),restflux,filterlist=restfilters))
          result[i].kcorr_mrest_synth = -2.5*alog10(restmaggies_synth) - restfiltinfo.vega2ab

          for ifilt = 0L, nrestfilter-1L do begin

; we have two options to compute k-corrections, with the first one
; being preferred because it minimizes the k-correction, although for
; some objects the magnitude will be defined but the magnitude error
; won't be: (1) take the nearest (redder?) observed bandpass that
; minimizes the k-correction, given the redshift; if no photometry has
; been measured, use the synthesized observed-frame magnitude; (2)
; search for the nearest bandpass with measured photometry, at the
; expense of increasing the k-correction

; -------------------------
; Method 1:             
; -------------------------

             case kcorr_method of

                1L: begin
                   
                   get_element, obsfiltinfo.weff/(1.0+kcorrphot[i].z), restfiltinfo[ifilt].weff, refindx

                   kcorr = -2.5*alog10(im_filtermag(restwave,restflux,filtinfo=restfiltinfo[ifilt])/$
                     im_filtermag(wave,flux,filtinfo=obsfiltinfo[refindx]))
;                  kcorr = -2.5*alog10(im_filtermag(restwave,restflux,filterlist=restfilters[ifilt])/$
;                    im_filtermag(wave,flux,filterlist=obsfilters[refindx]))
                   result[i].kcorr_kcorr[ifilt] = kcorr
                   result[i].kcorr_filter[ifilt] = obsbandpasses[refindx]

                   mobs = result[i].kcorr_mobs_ab[refindx]
                   if (mobs lt -900.0) then mobs = result[i].kcorr_mobs_ab_synth[refindx] ; <-- the big difference!
                   result[i].kcorr_mrest[ifilt] = mobs + kcorr - restfiltinfo[ifilt].vega2ab
;                  print, result[i].kcorr_mrest[ifilt], mobs, kcorr, $
;                    ' '+result[i].kcorr_filter[ifilt]+'-->'+restbandpasses[ifilt]
                   result[i].kcorr_mrest_err[ifilt] = result[i].kcorr_mobs_ab_err[refindx]

                end

; -------------------------
; Method 2:
; -------------------------

                2L: begin
                   
                   allrefindx = sort(abs(obsfiltinfo.weff/(1.0+kcorrphot[i].z)-restfiltinfo[ifilt].weff))
                   pickfilt = where(result[i].kcorr_mobs_ab[allrefindx] gt -900.0,nrefindx)
                   if (nrefindx eq 0L) then message, 'Problem here!'
                   refindx = allrefindx[pickfilt[0]]
                   
;                  get_element, filtinfo.weff/(1.0+kcorrphot[i].z), restfiltinfo[ifilt].weff, test
;                  if (test ne refindx) then print, result[i].kcorr_mobs_ab[test]

;                  niceprint, obsbandpasses, filtinfo.weff/(1.0+kcorrphot[i].z) & print, restfiltinfo[ifilt].weff

                   kcorr = -2.5*alog10(im_filtermag(restwave,restflux,filtinfo=restfiltinfo[ifilt])/$
                     im_filtermag(wave,flux,filtinfo=obsfiltinfo[refindx]))
;                  kcorr = -2.5*alog10(im_filtermag(restwave,restflux,filterlist=restfilters[ifilt])/$
;                    im_filtermag(wave,flux,filterlist=obsfilters[refindx]))
                   result[i].kcorr_kcorr[ifilt] = kcorr
                   result[i].kcorr_filter[ifilt] = obsbandpasses[refindx]

                   result[i].kcorr_mrest[ifilt] = result[i].kcorr_mobs_ab[refindx] + kcorr - restfiltinfo[ifilt].vega2ab
                   result[i].kcorr_mrest_err[ifilt] = result[i].kcorr_mobs_ab_err[refindx]

                end

             endcase
                
          endfor
;         print & cc = get_kbrd(1)

;         niceprint, restfilters, result[i].kcorr_kcorr, result[i].kcorr_filter

; debug plot          
          
          if keyword_set(debug) then begin

             if (ngood ne 0L) then begin

                xrange = [min(obsfiltinfo[good].weff-1.5*obsfiltinfo[good].fwhm),$
                  max(obsfiltinfo[good].weff+1.5*obsfiltinfo[good].fwhm)] ;/(1.0+zobj)
                xrange[0] = xrange[0]>90.0
                get_element, wave, xrange, xx

                if keyword_set(fnu) then begin
                   
                   yrange = [min(kcorrphot[i].maggies[good]-kcorrphot[i].maggies_err[good]),$
                     max(kcorrphot[i].maggies[good]+kcorrphot[i].maggies_err[good])]
                   yrange[0] = yrange[0]<min(flux_fnu[xx[0]:xx[1]])
                   yrange[1] = yrange[1]>max(flux_fnu[xx[0]:xx[1]])
                   
                   plot, [0], [0], /nodata, xsty=3, ysty=3, charsize=1.5, charthick=2.0, xthick=2.0, $
                     ythick=2.0, xtitle='Observed Wavelength [\AA]', ytitle=textoidl('f_{\nu}'), yrange=yrange, $
                     xrange=xrange, title=string(kcorrphot[i].galaxy,format='(A0)')+' z = '+$
                     strtrim(string(kcorrphot[i].z,format='(F12.3)'),2)
                   oploterror, obsfiltinfo[good].weff, kcorrphot[i].maggies[good], obsfiltinfo[good].fwhm, $
                     kcorrphot[i].maggies_err[good], ps=8
                   oplot, wave, flux_fnu, line=0.1

                endif else begin

                   yrange = fltarr(2)
                   yrange[0] = max(-2.5*alog10(kcorrphot[i].maggies[good]))>max(flux_ab[xx[0]:xx[1]]*1.00)
                   yrange[1] = min(-2.5*alog10(kcorrphot[i].maggies[good]))<min(flux_ab[xx[0]:xx[1]]*0.98)

                   plot, [0], [0], /nodata, xsty=3, ysty=3, charsize=1.5, charthick=2.0, xthick=2.0, $
                     ythick=2.0, xtitle='Observed Wavelength [\AA]', ytitle=textoidl('m_{AB}'), $
                     yrange=yrange, xrange=xrange, title=string(kcorrphot[i].galaxy,format='(A0)')+' z = '+$
                     strtrim(string(kcorrphot[i].z,format='(F12.3)'),2)
                   oploterror, obsfiltinfo[good].weff, -2.5*alog10(kcorrphot[i].maggies[good]), obsfiltinfo[good].fwhm, $
                     2.5*kcorrphot[i].maggies_err[good]/kcorrphot[i].maggies[good]/alog(10.0), ps=8
                   oplot, wave, flux_ab, line=0.1

                endelse
                   
                label = textoidl('\chi^{2}_{\nu} = '+strtrim(string(result[i].kcorr_chi2,format='(F12.2)'),2))
                legend, label, /right, /bottom, box=0, charsize=1.5, charthick=2.0

               cc = get_kbrd(1)

             endif
         
          endif

      endif 

    endfor
    splog, 'Total time = '+string((systime(1)-t0)/60.0,format='(G0.0)')+' minutes.'

; now expand the rest-frame photometry into individual structure tags

    morephot = create_struct(restbandpasses[0],-999.0,restbandpasses[0]+'_err',-999.0,'M_'+restbandpasses[0],-999.0,$
      'M_'+restbandpasses[0]+'_err',-999.0,restbandpasses[0]+'_lum',-999.0,restbandpasses[0]+'_lum_err',-999.0)
    for i = 1L, nrestfilter-1L do morephot = create_struct(morephot,restbandpasses[i],-999.0,restbandpasses[i]+'_err',-999.0,$
      'M_'+restbandpasses[i],-999.0,'M_'+restbandpasses[i]+'_err',-999.0,restbandpasses[i]+'_lum',-999.0,$
      restbandpasses[i]+'_lum_err',-999.0)

    result = struct_addtags(result,replicate(morephot,ngalaxy))
    
; compute absolute magnitudes and luminosities; temporarily append a
; distance tag; DISTANCE ERR should be zero
    
    splog, 'Computing absolute magnitudes and luminosities.'
    dtags = replicate({distance: -999.0, distance_err: 0.0},ngalaxy)
    good = where(result.z gt 0.0,ngood)
    if (ngood ne 0L) then dtags[good].distance = dluminosity(result[good].z,/mpc)
    result = struct_addtags(result,dtags)

    for iband = 0L, nrestfilter-1L do begin

       true = tag_exist(result,restbandpasses[iband],index=photindx)
       true = tag_exist(result,restbandpasses[iband]+'_err',index=photindx_err)
       true = tag_exist(result,'M_'+restbandpasses[iband],index=absphotindx)
       true = tag_exist(result,'M_'+restbandpasses[iband]+'_err',index=absphotindx_err)
       true = tag_exist(result,restbandpasses[iband]+'_lum',index=lumphotindx)
       true = tag_exist(result,restbandpasses[iband]+'_lum_err',index=lumphotindx_err)

       result.(photindx) = result.kcorr_mrest[iband]
       result.(photindx_err) = result.kcorr_mrest_err[iband]

; NOTE! absolute magnitudes are not computed for all objects if
; KCORR_METHOD = 1 since the magnitude errors for these objects will
; be -999.0, even though the magnitudes are not -999.0
       
       good = where((result.distance gt -900.0) and (result.(photindx) gt -900.0) and (result.(photindx_err) gt -900.0),ngood)
       if (ngood ne 0L) then begin
          mags = im_absolute_magnitudes(restbandpasses[iband],result[good].(photindx),$
            result[good].distance,mag_err=result[good].(photindx_err),$
            distance_err=result[good].distance_err)
          if (size(mags,/type) ne 8L) then message, 'Problem here.'
          result[good].(absphotindx)     = mags.absmag
          result[good].(absphotindx_err) = mags.absmag_err
          result[good].(lumphotindx)     = mags.lum
          result[good].(lumphotindx_err) = mags.lum_err
       endif

    endfor 

    result = struct_trimtags(result,except=['DISTANCE','DISTANCE_ERR'])

; compute color-based stellar masses and append

    splog, 'Computing color-based stellar masses.'
    mass = compute_stellar_mass(result)
    result = struct_addtags(result,mass)

return, result
end
