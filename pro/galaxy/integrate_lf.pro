;+
; NAME:
;   INTEGRATE_LF()
;
; PURPOSE:
;   Integrate a 1/Vmax-weighted luminosity function and, optionally,
;   the best-fitting Schechter model. 
;
; INPUTS: 
;   lf_vmax - 1/Vmax weighted luminosity function in the style of
;     IM_MF_VMAX() 
;
; OPTIONAL INPUTS: 
;   schechter - input IM_LF_SCHECHTER-style data structure
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   result - output data structure with lots of goodies
;
; COMMENTS:
;   Error checking isn't great.
;
;   Here are the analytic formulae for some of the integrals: 
;     rhotot = Phistar*Mstar*Gamma(alpha+2)
;     numtot = Phistar*Gamma(alpha+1)
;     rho(>M) = Phistar*Mstar*IGamma(alpha+2,M/Mstar)
;     num(>M) = Phistar*IGamma(alpha+1,M/Mstar)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Aug 15, UCSD
;
; Copyright (C) 2011, John Moustakas
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

function integrate_lf, lf_vmax, absmagfaint=absmagfaint1, $
  absmagbright=absmagbright1, schechter=schechter, nmonte=nmonte

    if (n_elements(lf_vmax) eq 0) and (n_elements(schechter) eq 0) then begin
       doc_library, 'integrate_lf'
       return, -1
    endif

    if (n_elements(nmonte) eq 0) then nmonte = 300

; initialize the output data structure
    result = {absmagfaint: -1.0, absmagbright: -1.0, absmagfaint_data: -1.0, absmagbright_data: -1.0, $
      absmagfaint_model: -1.0, absmagbright_model: -1.0, absmagfaint_tot: -1.0, absmagbright_tot: -1.0, $
      rho: -1.0, rho_cor: -1.0, rho_model: -1.0, rhotot_model: -1.0, rhoerr: -1.0, $
      rhoerr_model: -1.0, rhototerr_model: -1.0, num: -1.0, num_cor: -1.0, num_model: -1.0, $
      numerr: -1.0, numerr_model: -1.0}

; integrate the LF over the measured range of absolute magnitudes,
; using the model, if provided, to extrapolate to ABSMAG=[-10,-30];
; get the errors by Monte Carlo
    if (n_elements(lf_vmax) ne 0) then begin
       gd = where(lf_vmax.limit eq 1,ngd)
       absmagfaint_data = min(lf_vmax.absmag[gd])
       absmagbright_data = max(lf_vmax.absmag[gd])
       result.absmagfaint_data = absmagfaint_data
       result.absmagbright_data = absmagbright_data

       if (n_elements(absmagfaint1) eq 0) then absmagfaint = absmagfaint_data else absmagfaint = absmagfaint1
       if (n_elements(absmagbright1) eq 0) then absmagbright = absmagbright_data else absmagbright = absmagbright1
       if (absmagfaint_data le absmagfaint1) and (absmagbright_data ge absmagfaint1) then begin
          result.absmagfaint = absmagfaint
          result.absmagbright = absmagbright
          
          these = where((lf_vmax.limit eq 1) and (lf_vmax.absmag ge absmagfaint) and $
            (lf_vmax.absmag le absmagbright),nthese)
          if (nthese gt 1) then begin ; need at least two points
             absmag = lf_vmax.absmag[these]
             phi = lf_vmax.phi[these]
             phierr = lf_vmax.phierr[these]

             result.rho = im_integral(absmag,10D^absmag*phi,absmagfaint,absmagbright) ; [Msun Mpc^-3]
             result.num = im_integral(absmag,phi,absmagfaint,absmagbright)          ; [Mpc^-3]

             nbin = n_elements(absmag)
             rhomonte = fltarr(nmonte)
             nummonte = fltarr(nmonte)
             phimonte = randomn(seed,nbin,nmonte)*rebin(reform(phierr,nbin,1),nbin,nmonte)+$
               rebin(reform(phi,nbin,1),nbin,nmonte)
             for ii = 0L, nmonte-1 do begin
                rhomonte[ii] = im_integral(absmag,10D^absmag*phimonte[*,ii],absmagfaint,absmagbright)
                nummonte[ii] = im_integral(absmag,phimonte[*,ii],absmagfaint,absmagbright)
             endfor
             result.rhoerr = djsig(rhomonte)
             result.numerr = djsig(nummonte)
          endif else begin
             result.rho = 10D^lf_vmax.absmag[these]*lf_vmax.phi[these]*lf_vmax.binsize
             result.num = lf_vmax.phi[these]*lf_vmax.binsize
             result.rhoerr = 10D^lf_vmax.absmag[these]*lf_vmax.phierr[these]*lf_vmax.binsize
             result.numerr = lf_vmax.phierr[these]*lf_vmax.binsize
          endelse
; if necessary, "correct" the measured luminosity densities and number
; densities using the model to extrapolate to ABSMAGBRIGHT
          if (n_elements(schechter) ne 0) and (result.absmagbright_data lt absmagbright1) then begin
             modelabsmag = range(result.absmagbright_data,absmagbright1,500)
             modelphi = im_lf_schechter(modelabsmag,schechter=schechter)
             rho_add = im_integral(modelabsmag,10D^modelabsmag*modelphi,result.absmagbright_data,absmagbright1)
             num_add = im_integral(modelabsmag,modelphi,result.absmagbright_data,absmagbright1)
             
             if (result.rho le 0) then message, 'Tenemos problema, chico!'
             result.rho_cor = result.rho + rho_add
             result.num_cor = result.num + num_add
          endif
       endif
    endif 
    
; integrate the model    
    if (n_elements(schechter) ne 0) then begin
       if (n_elements(absmagfaint1) eq 0) then absmagfaint = -10D else absmagfaint = absmagfaint1
       if (n_elements(absmagbright1) eq 0) then absmagbright = -30D else absmagbright = absmagbright1

       result.absmagfaint_model = absmagfaint
       result.absmagbright_model = absmagbright       
       modelabsmag = range(absmagfaint,absmagbright,1000)

; need to convert to luminosity       
       modellum = -0.4*(modelabsmag+20.0)
       minlum = -0.4*(absmagfaint+20)
       maxlum = -0.4*(absmagbright+20)

; use empirical integration because we can integrated the double
; Schechter function easily; the empirical and analytic results match 
       modelphi = im_lf_schechter(modelabsmag,schechter=schechter)
       result.rho_model = -2.5*alog10(im_integral(modellum,10^modellum*modelphi,minlum,maxlum))-20.0

stop       
       
       result.num_model = im_integral(modellum,modelphi,absmagbright,absmagfaint)

; also integrate the full Schechter function       
       absmagfainttot = 7D
       absmagbrighttot = 15D
       result.absmagfaint_tot = absmagfainttot
       result.absmagbright_tot = absmagbrighttot
       modelabsmagtot = range(absmagfainttot,absmagbrighttot,1000)

       modelphitot = im_lf_schechter(modelabsmagtot,schechter=schechter)
       result.rhotot_model = im_integral(modelabsmagtot,10D^modelabsmagtot*modelphitot,absmagfainttot,absmagbrighttot)

; analytically:       
;      result.rho_model = schechter.phistar*10D^schechter.logmstar*gamma(schechter.alpha+2)*$
;        (1-igamma(schechter.alpha+2,10D^(absmagfaint-schechter.logmstar)))
;      result.num_model = schechter.phistar*gamma(schechter.alpha+2)*$
;        (1-igamma(schechter.alpha+2,10D^(absmagfaint-schechter.logmstar)))
;      result.rhotot_model = schechter.phistar*10D^schechter.logmstar*gamma(schechter.alpha+2)
    endif 
    
return, result
end
