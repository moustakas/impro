;+
; NAME:
;   ISEDFIT_RECONSTRUCT_SFH()
;
; PURPOSE:
;   Reconstruct an iSEDfit-compatible star formation history. 
;
; INPUTS: 
;   info - iSEDfit information structure with the following mandatory
;     tags (see BUILD_ISEDFIT_SFHGRID) for definitions and details 
;     NBURST
;     FBURST
;     TBURST
;     DTBURST
;     MINTBURST
;     MAXTBURST
; 
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
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

function isedfit_reconstruct_sfh, info, useage=useage, outage=outage, mtau=mtau, $
  aburst=aburst, mburst=mburst, mgalaxy=outmgalaxy, sfr100=outsfr100, $
  b100=outb100, notruncate=notruncate, sfhtau=outsfhtau, sfhburst=outsfhburst, $
  sfrage=outsfrage, linear=linear, debug=debug, _extra=extra, nagemax=nagemax, $
  dage=dage
; jm10dec22ucsd - given an iSEDfit structure, reconstruct the star
; formation history, allowing for multiple bursts

; recommend you only use AGE on output, unless the array was
; constructed using BUILD_ISEDFIT_AGEGRID

; MTAU is the normalization of the model; default to 1 Msun    

; notruncate allows you to override a truncated burst (used by
; BUILD_ISEDFIT_AGEGRID)
    
; note that if a truncated burst is requested then SFHTAU and SFHBURST
; outputs are not modified
    
    if (n_elements(info) eq 0) then begin
       doc_library, 'isedfit_reconstruct_sfh'
       return, -1
    endif

; burst preliminaries (need this info up here)
    nb = info.nburst 
    if (nb gt 0) then begin
       fburst = im_double(info.fburst[0:nb-1])
       tburst = im_double(info.tburst[0:nb-1])
       dtburst = im_double(info.dtburst[0:nb-1])
    endif
    
; need a highly sampled age grid to get the integrations right; this
; grid is used internally and then the output quantities are
; interpolated at OUTAGE
    if (n_elements(nagemax) eq 0) then nagemax = 300L
    if (n_elements(dage) eq 0) then dage = 0.02D ; default 20 Myr time interval

    minage = 0D
    if (nb eq 0) then maxage = im_double(info.maxage) else $
      maxage = im_double(info.maxage>info.maxtburst)
    nage = long(((maxage-minage)/dage)<nagemax)
    
    if (n_elements(useage) eq 0) then begin
       age = build_isedfit_agegrid(info,debug=0,nage=nage,$
         minage=minage,maxage=maxage,linear=linear)
    endif else age = useage

    if (n_elements(outage) eq 0) then outage = age
    nage = n_elements(age)

; deal with the simple tau model, including delayed tau models, and
; with bursty SFHs     
    if (n_elements(mtau) eq 0) then mtau = 1D ; normalization
    if (info.tau eq 0D) then sfhtau = dblarr(nage) else begin
       if info.delayed then sfhtau = mtau*age*1D9*exp(-age/info.tau)/(info.tau*1D9)^2 else $
         sfhtau = mtau*exp(-age/info.tau)/(info.tau*1D9)
    endelse

; type of burst: 0=step function (default); 1=gaussian; 2=step
; function with exponential wings 
    if (nb gt 0) then begin
       sfhburst1 = reform(dblarr(nage,nb),nage,nb)
; compute the amplitude, and then build each burst in turn
       aburst = dblarr(nb)
       mburst = dblarr(nb)
       for ib = 0, nb-1 do begin
          if (info.tau eq 0D) then $
            aburst[ib] = fburst[ib]*mtau/(dtburst[ib]*1D9) else $
              aburst[ib] = fburst[ib]*mtau*(1.0-exp(-tburst[ib]/info.tau))/(dtburst[ib]*1D9)
; step-function burst (default)
          if info.bursttype eq 0 then begin
             if (max(age) ge tburst[ib]) then begin
                t1 = (findex(age,tburst[ib]))>0
                t2 = (findex(age,tburst[ib]+dtburst[ib]))<(nage-1)
                sfhburst1[t1:t2,ib] = aburst[ib]
                mburst[ib] = aburst[ib]*(interpolate(age,t2)-interpolate(age,t1))*1D9 ; =dtburst[ib] except on the edges
             endif
          endif
; Gaussian burst; could do the integral analytically, but no need to
          if info.bursttype eq 1 then begin
             sfhburst1[*,ib] = aburst[ib]*exp(-0.5*((age-tburst[ib])/dtburst[ib])^2)/sqrt(2.0*!pi) ; [Msun/yr]
             mburst[ib] = im_integral(age*1D9,sfhburst1[*,ib]) ; [Msun]
          endif
; step-function burst with exponential wings; the MBURST() integral
; won't be quite right because of the discontinuity
          if info.bursttype eq 2 then begin
             during = where((age ge tburst[ib]) and (age le tburst[ib]+dtburst[ib]),nduring)
             before = where(age lt tburst[ib],nbefore)
             after = where(age gt tburst[ib]+dtburst[ib],nafter)
             if (nbefore ne 0) then sfhburst1[before,ib] += aburst[ib]*exp(-(tburst[ib]-age[before])/0.01D)
             if (nafter ne 0) then sfhburst1[after,ib] += aburst[ib]*exp(-(age[after]-(tburst[ib]+dtburst[ib]))/0.01D)
             if (nduring ne 0) then sfhburst1[during,ib] += aburst[ib]
             mburst[ib] = im_integral(age*1D9,sfhburst1[*,ib]) ; [Msun]
          endif
       endfor 
       sfhburst = total(sfhburst1,2,/double) ; add 'em up!
    endif else sfhburst = sfhtau*0D

    sfh = sfhtau + sfhburst ; final star formation history

; truncate the last burst?
    dotruncate = 0
    ilast = -1
    if (nb gt 0) and (info.trunctau gt 0D) then begin
       if (keyword_set(notruncate) eq 0) then begin
          dotruncate = 1
          ilast = (long(findex(age,tburst[nb-1])))>0
          if ilast lt 0 then message, 'This should not happen'
          if info.bursttype eq 1 then $ ; Gaussian burst shape
            post = where(age ge tburst[nb-1],npost) else $ ; after the peak of the burst
              post = where(age ge tburst[nb-1]+dtburst[nb-1],npost) ; after the full width of the burst
          if (npost gt 0) then begin
; truncate the full SFH
             sfrpostburst = interpol(sfh,age,tburst[nb-1])
             sfh[post] = sfrpostburst*exp(-(age[post]-tburst[nb-1])/info.trunctau)
             
; code below truncates SFHBURST and SFHTAU, and adjusts MBURST for the
; truncation; in general we don't use these quantities, so skip
; the extra work              
;             sfrpostburst1 = interpol(sfhburst1[*,nb-1],age,tburst[nb-1])
;             sfhburst1[post,nb-1] = sfrpostburst1*exp(-(age[post]-tburst[nb-1])/info.trunctau)
;             mburst[nb-1] = im_integral(age*1D9,sfhburst1[*,nb-1]) ; [Msun]
;;            mburst[nb-1] = 0.5D*mburst[nb-1] + sfrpostburst*info.trunctau[nb-1]*1D9*$
;;              exp(-tburst[nb-1]/info.trunctau[nb-1]) ; [Msun]                
;             sfhtau[post] = 0
;             sfhburst = total(sfhburst1,2,/double)
          endif
       endif
    endif

; if requested, compute the <SFR> over the previous 100 Myr, the
; birthrate parameter, and the SFH-weighted age
    if arg_present(outmgalaxy) or arg_present(outsfr100) or $
      arg_present(outb100) or arg_present(outsfrage) then begin
       dt = 0.1D ; [100 Myr]
       mgalaxy = sfh*0D
       sfr100 = sfh*0D
       b100 = sfh*0D
       sfrage = sfh*0D
       for iage = 0L, nage-1 do begin
; if the SFH has been truncated then do the integrals, otherwise
; calculate the burst masses with integrals and the tau-model
; analytically
          if dotruncate and (iage gt ilast) then begin
             mtot100 = im_integral(age*1D9,sfh,1D9*(age[iage]-dt)>0,1D9*age[iage])
;if n_elements(uniq(age,sort(age))) ne 500 then stop
;if age[ilast] eq age[iage] then stop
;stop
             mgalaxy[iage] = mgalaxy[ilast] + im_integral(age*1D9,sfh,1D9*age[ilast],1D9*age[iage])
          endif else begin
             if (nb eq 0) then begin
                m100burst = 0D
                mtotburst = 0D
             endif else begin
                m100burst = im_integral(age*1D9,sfhburst,1D9*(age[iage]-dt)>0,1D9*age[iage])
                mtotburst = im_integral(age*1D9,sfhburst,0D,1D9*age[iage])
             endelse
             if (info.tau eq 0D) then begin
                mtot100 = 0D + m100burst
                mgalaxy[iage] = mtau + mtotburst
             endif else begin
                mtot100 = mtau*(exp(-((age[iage]-dt)>0)/info.tau)-exp(-age[iage]/info.tau)) + m100burst
                mgalaxy[iage] = mtau*(1D0-exp(-age[iage]/info.tau)) + mtotburst
             endelse
          endelse
          if (age[iage] gt 0D) then begin
             sfr100[iage] = mtot100/(dt*1D9)
             b100[iage] = sfr100[iage]/(mgalaxy[iage]/(1D9*age[iage]))
          endif
; compute the SFR-weighted age
          if (info.tau eq 0D) then sfrage[iage] = age[iage] else begin
             norm = im_integral(age*1D9,sfh,0D,1D9*age[iage])
             if (norm eq 0D) then sfrage[iage] = age[iage] else $
               sfrage[iage] = im_integral(age*1D9,sfh*(age[iage]-age)*1D9,0D,1D9*age[iage])/norm
          endelse
          
;         print, age[iage], mtot100, mtot100/(dt*1D9), sfh[iage]
;         plot, age, sfhburst, xr=[1.4,2.4], psym=-6, xsty=3, ysty=3
;         djs_oplot, age[iage]*[1,1], !y.crange, color='red'
;         djs_oplot, (age[iage]-dt)*[1,1], !y.crange, color='blue'
;         djs_oplot, age[0:iage], sfr100[0:iage], psym=-6, color='cyan'
;         if (outage[iage] ge 1.65) then cc = get_kbrd(1) ; stop
       endfor
    endif

; interpolate onto OUTAGE
    findx = findex(age,outage)
    outsfh = interpolate(sfh,findx)
;   outsfh = exp(interpolate(alog(sfh),findx,cubic=-1))
    outsfhtau = interpolate(sfhtau,findx)
    outsfhburst = interpolate(sfhburst,findx)

    if arg_present(outmgalaxy) then outmgalaxy = interpolate(mgalaxy,findx)
    if arg_present(outsfr100) then outsfr100 = interpolate(sfr100,findx)
    if arg_present(outb100) then outb100 = interpolate(b100,findx)
    if arg_present(outsfrage) then outsfrage = interpolate(sfrage,findx)/1D9 ; [Gyr]

; QAplot    
    if keyword_set(debug) then begin
       djs_plot, age, sfh, xlog=0, xsty=3, ysty=3, psym=-6, _extra=extra;, xr=[min(age)>0.01,max(age)]
;      djs_oplot, age, sfr100, psym=6, color='cyan'
       djs_oplot, outage, outsfh, psym=6, color='orange'
;      djs_oplot, outage, outsfr100, psym=6, color='blue'
;      djs_oplot, outage, outsfhtau, color='blue', psym=-6, sym=0.5
;      djs_oplot, outage, outsfhburst, color='red', psym=-6, sym=0.5
       if dotruncate then djs_oplot, tburst[nb-1]+info.trunctau*[1,1], !y.crange, color='yellow'
;      for ib = 0, nb-1 do djs_oplot, tburst[ib]*[1,1], !y.crange, color='yellow'
    endif

return, outsfh
end
