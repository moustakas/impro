;+
; NAME:
;   INTEGRATE_MF()
;
; PURPOSE:
;   Integrate a 1/Vmax-weighted stellar mass function and, optionally,
;   the best-fitting Schechter model.
;
; INPUTS: 
;   mf_vmax - 1/Vmax weighted stellar mass function in the style of
;     IM_MF_VMAX()
;
; OPTIONAL INPUTS: 
;   schechter - input MF_SCHECHTER-style data structure
;
; KEYWORD PARAMETERS: 
;   double - assume that SCHECHTER describes the double Schechter
;     function used by MF_DOUBLE_SCHECHTER
;   plus - assume that SCHECHTER describes the modified Schechter
;     function used by MF_SCHECHTER_PLUS
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

function integrate_mf, mf_vmax, minmass=minmass1, maxmass=maxmass1, $
  schechter=schechter, nmonte=nmonte, plus=plus, double=double, debug=debug

    if (n_elements(mf_vmax) eq 0) and (n_elements(schechter) eq 0) then begin
       doc_library, 'integrate_mf'
       return, -1
    endif

    if (n_elements(nmonte) eq 0) then nmonte = 300

; initialize the output data structure
    result = {minmass: -1.0, maxmass: -1.0, minmass_data: -1.0, maxmass_data: -1.0, $
      minmass_model: -1.0, maxmass_model: -1.0, minmass_tot: -1.0, maxmass_tot: -1.0, $
      rho: -1.0, rho_cor: -1.0, rho_model: -1.0, rhotot_model: -1.0, rhoerr: -1.0, $
      rhoerr_model: -1.0, rhototerr_model: -1.0, num: -1.0, num_cor: -1.0, num_model: -1.0, $
      numtot_model: -1.0, numerr: -1.0, numerr_model: -1.0, numtoterr_model: -1.0, $
      mass_cumuphi50: -1.0, mass_cumuphi45: -1.0, mass_cumuphi40: -1.0, mass_cumuphi35: -1.0, $
      mass_cumuphi30: -1.0, mass_cumuphi25: -1.0}

; integrate the MF over the measured range of stellar masses, using
; the model, if provided, to extrapolate to from log(M)=[5,15]; get
; the errors by Monte Carlo
    if (n_elements(mf_vmax) ne 0) then begin
       gd = where(mf_vmax.limit eq 1 and mf_vmax.number gt 3,ngd)
;      gd = where(mf_vmax.limit eq 1,ngd)
       minmass_data = min(mf_vmax.meanmass[gd])
       maxmass_data = max(mf_vmax.meanmass[gd])
       result.minmass_data = minmass_data
       result.maxmass_data = maxmass_data

       if (n_elements(minmass1) eq 0) then minmass = minmass_data else minmass = minmass1
       if (n_elements(maxmass1) eq 0) then maxmass = maxmass_data else maxmass = maxmass1
       if (minmass_data le minmass) and (maxmass_data ge minmass) then begin
          result.minmass = minmass
          result.maxmass = maxmass
          
          these = where((mf_vmax.limit eq 1) and (mf_vmax.meanmass ge minmass) and $
            (mf_vmax.meanmass le maxmass) and (mf_vmax.number gt 3),nthese)
          if (nthese eq 0) then begin ; special case
             result.rho = 0D
             result.num = 0D
             result.rhoerr = 0D
             result.numerr = 0D
          endif else begin
             if (nthese gt 1) then begin ; need at least two points
                mass = mf_vmax.meanmass[these]
                phi = mf_vmax.phi[these]
                phierr = mf_vmax.phierr[these]

                result.rho = im_integral(mass,10D^mass*phi,minmass,maxmass) ; [Msun Mpc^-3]
                result.num = im_integral(mass,phi,minmass,maxmass)          ; [Mpc^-3]

                nbin = n_elements(mass)
                rhomonte = fltarr(nmonte)
                nummonte = fltarr(nmonte)
                phimonte = randomn(seed,nbin,nmonte)*rebin(reform(phierr,nbin,1),nbin,nmonte)+$
                  rebin(reform(phi,nbin,1),nbin,nmonte)
                for ii = 0L, nmonte-1 do begin
                   rhomonte[ii] = im_integral(mass,10D^mass*phimonte[*,ii],minmass,maxmass)
                   nummonte[ii] = im_integral(mass,phimonte[*,ii],minmass,maxmass)
                endfor
                result.rhoerr = djsig(rhomonte)
                result.numerr = djsig(nummonte)
             endif else begin
                result.rho = 10D^mf_vmax.meanmass[these]*mf_vmax.phi[these]*mf_vmax.binsize
                result.num = mf_vmax.phi[these]*mf_vmax.binsize
                result.rhoerr = 10D^mf_vmax.meanmass[these]*mf_vmax.phierr[these]*mf_vmax.binsize
                result.numerr = mf_vmax.phierr[these]*mf_vmax.binsize
             endelse
          endelse
; if necessary, "correct" the measured mass densities and number
; densities using the model to extrapolate to MAXMASS
          if (n_elements(schechter) ne 0) then begin
             if (result.maxmass_data lt maxmass1) then begin
                modelmass = range(result.maxmass_data,maxmass1,100)
                if (keyword_set(double) eq 0) and (keyword_set(plus) eq 0) then begin
                   modelphi = mf_schechter(modelmass,schechter=schechter)
                endif else begin
                   if keyword_set(double) then modelphi = mf_double_schechter(modelmass,schechter=schechter)
                   if keyword_set(plus) then modelphi = mf_schechter_plus(modelmass,schechter=schechter)
                endelse
                rho_add = im_integral(modelmass,10D^modelmass*modelphi,result.maxmass_data,maxmass1)
                num_add = im_integral(modelmass,modelphi,result.maxmass_data,maxmass1)
                
                if (result.rho lt 0) or (result.num lt 0) then message, 'Tenemos problema, chico!'
                result.rho_cor = result.rho + rho_add
                result.num_cor = result.num + num_add                               
             endif else begin
                result.rho_cor = result.rho
                result.num_cor = result.num
             endelse
; compute the stellar mass at which the *cumulative* number density
; equals various thresholds
             nnmodel = 100
             modelmass = range(7D,14D,nnmodel)
;            modelmass = range(result.maxmass_data,14D,100)
             if (keyword_set(double) eq 0) and (keyword_set(plus) eq 0) then begin
                modelphi = mf_schechter(modelmass,schechter=schechter)
             endif else begin
                if keyword_set(double) then modelphi = mf_double_schechter(modelmass,schechter=schechter)
                if keyword_set(plus) then modelphi = mf_schechter_plus(modelmass,schechter=schechter)
             endelse

; integrate the model             
             model_cumuphi = modelphi*0
             for ii = nnmodel-2, 0, -1 do model_cumuphi[ii] = im_integral(modelmass,$
               modelphi,modelmass[ii],modelmass[nnmodel-1])

; integrate the data             
             keep = where(modelmass gt max(mf_vmax.meanmass[gd]) or $
               modelmass lt min(mf_vmax.meanmass[gd]))
             mass1 = [mf_vmax.meanmass[gd],modelmass[keep]]
             phi1 = [mf_vmax.phi[gd],modelphi[keep]]
             uu = uniq(mass1,sort(mass1))
             mass1 = mass1[uu]
             phi1 = phi1[uu]
             cumuphi = phi1*0
             nn = n_elements(mass1)
             for ii = nn-2, 0, -1 do cumuphi[ii] = im_integral(mass1,phi1,mass1[ii],mass1[nn-1])
             mm = interpol(mass1,alog10(cumuphi>1D-30),[-5.0,-4.5,-4.0,-3.5,-3.0,-2.5])

             if keyword_set(debug) then begin
                djs_plot, alog10(cumuphi>1D-30), mass1, psym=6, xsty=3, ysty=3, xrange=[-1.5,-6], yrange=[7,12.2]
                djs_oplot, alog10(model_cumuphi>1D-30), modelmass, psym=6, color='orange'
                djs_oplot, alog10(model_cumuphi[keep]>1D-30), modelmass[keep], psym=6, color='blue'
                djs_oplot, -5*[1,1], !y.crange, line=0
                djs_oplot, -4*[1,1], !y.crange, line=5
                djs_oplot, -3*[1,1], !y.crange, line=3
                cc = get_kbrd(1)
             endif

             result.mass_cumuphi50 = mm[0]
             result.mass_cumuphi45 = mm[1]
             result.mass_cumuphi40 = mm[2]
             result.mass_cumuphi35 = mm[3]
             result.mass_cumuphi30 = mm[4]
             result.mass_cumuphi25 = mm[5]
          endif  
       endif
    endif 
    
; integrate the model    
    if (n_elements(schechter) ne 0) then begin
       if (n_elements(minmass1) eq 0) then minmass = 7D else minmass = minmass1
       if (n_elements(maxmass1) eq 0) then maxmass = 14D else maxmass = maxmass1

       result.minmass_model = minmass
       result.maxmass_model = maxmass       
       modelmass = range(minmass,maxmass,500)

; use empirical integration because we can integrated the double
; Schechter function easily; the empirical and analytic results match 
       if (keyword_set(double) eq 0) and (keyword_set(plus) eq 0) then begin
          modelphi = mf_schechter(modelmass,schechter=schechter)
       endif else begin
          if keyword_set(double) then modelphi = mf_double_schechter(modelmass,schechter=schechter)
          if keyword_set(plus) then modelphi = mf_schechter_plus(modelmass,schechter=schechter)
       endelse
       result.rho_model = im_integral(modelmass,10D^modelmass*modelphi,minmass,maxmass)
       result.num_model = im_integral(modelmass,modelphi,minmass,maxmass)

; also integrate the full Schechter function       
       minmasstot = 7D
       maxmasstot = 14D
       result.minmass_tot = minmasstot
       result.maxmass_tot = maxmasstot
       modelmasstot = range(minmasstot,maxmasstot,500)

       if (keyword_set(double) eq 0) and (keyword_set(plus) eq 0) then begin
          modelphitot = mf_schechter(modelmasstot,schechter=schechter)
       endif else begin
          if keyword_set(double) then modelphitot = mf_double_schechter(modelmasstot,schechter=schechter)
          if keyword_set(plus) then modelphitot = mf_schechter_plus(modelmasstot,schechter=schechter)
       endelse
       result.rhotot_model = im_integral(modelmasstot,10D^modelmasstot*modelphitot,minmasstot,maxmasstot)
       result.numtot_model = im_integral(modelmasstot,modelphitot,minmasstot,maxmasstot)

; analytically:       
;      result.rho_model = schechter.phistar*10D^schechter.logmstar*gamma(schechter.alpha+2)*$
;        (1-igamma(schechter.alpha+2,10D^(minmass-schechter.logmstar)))
;      result.num_model = schechter.phistar*gamma(schechter.alpha+2)*$
;        (1-igamma(schechter.alpha+2,10D^(minmass-schechter.logmstar)))
;      result.rhotot_model = schechter.phistar*10D^schechter.logmstar*gamma(schechter.alpha+2)
    endif 
    
return, result
end
