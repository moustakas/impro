;+
; NAME:
;   ISEDFIT_RECONSTRUCT_SFH()
;
; PURPOSE:
;   
;
; INPUTS: 
;
;
; OPTIONAL INPUTS: 
;
;
; KEYWORD PARAMETERS: 
;
;
; OUTPUTS: 
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
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

function isedfit_reconstruct_sfh, info, outage=outage, mtau=mtau, $
  aburst=aburst, mburst=mburst, mgalaxy=outmgalaxy, sfr100=outsfr100, $
  b100=outb100, notruncate=notruncate, sfhtau=outsfhtau, sfhburst=outsfhburst, $
  nooversample=nooversample, debug=debug, gaussburst=gaussburst, $
  stepburst=stepburst, _extra=extra
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
    
; need a highly sampled age grid to get the integrations right
;   if (n_elements(age) eq 0) then age = build_isedfit_agegrid(info,$
;     debug=0,_extra=extra)
    if (n_elements(outage) ne 0) then begin
       if (nb gt 0) then begin
          minage = min(outage)<min(tburst) ; burst can technically start before MINAGE
          maxage = max(outage)>max(tburst)
       endif else begin
          minage = min(outage)
          maxage = max(outage)
       endelse
       nage = 500>n_elements(outage)
    endif else nage = 500

    if keyword_set(nooversample) then begin
       if (n_elements(outage) eq 0) then message, 'NOOVERSAMPLE '+$
         'keyword requires OUTAGE!'
       age = outage
    endif else begin
       age = build_isedfit_agegrid(info,debug=0,nage=500,$
         minage=minage,maxage=maxage,linear=linear)
    endelse
    if (n_elements(outage) eq 0) then outage = age
    nage = n_elements(age)

; deal with the simple tau model and with bursty SFHs    
    if (n_elements(mtau) eq 0) then mtau = 1D ; normalization
    if (info.tau eq 0.0) then sfhtau = dblarr(nage) else $
      sfhtau = mtau*exp(-age/info.tau)/(info.tau*1D9) 

    if (nb gt 0) then begin
       sfhburst1 = reform(dblarr(nage,nb),nage,nb)
; compute the amplitude, and then build each burst in turn
       aburst = dblarr(nb)
       mburst = dblarr(nb)
       for ib = 0, nb-1 do begin
          if (info.tau eq 0.0) then $
            aburst[ib] = fburst[ib]*mtau/(dtburst[ib]*1D9) else $
              aburst[ib] = fburst[ib]*mtau*(1.0-exp(-tburst[ib]/info.tau))/(dtburst[ib]*1D9)

; step-function burst with exponential wings (default)
          if (keyword_set(gaussburst) eq 0) and (keyword_set(stepburst) eq 0) then begin
             during = where((age ge tburst[ib]) and (age le tburst[ib]+dtburst[ib]),nduring)
             before = where(age lt tburst[ib],nbefore)
             after = where(age gt tburst[ib]+dtburst[ib],nafter)
             if (nbefore ne 0) then sfhburst1[before,ib] += aburst[ib]*exp(-(tburst[ib]-age[before])/0.01D)
             if (nafter ne 0) then sfhburst1[after,ib] += aburst[ib]*exp(-(age[after]-(tburst[ib]+dtburst[ib]))/0.01D)
             if (nduring ne 0) then sfhburst1[during,ib] += aburst[ib]
             mburst[ib] = im_integral(age*1D9,sfhburst1[*,ib]) ; [Msun]
          endif
; step-function burst
          if keyword_set(stepburst) then begin
             if (max(age) ge tburst[ib]) then begin
                t1 = (findex(age,tburst[ib]))>0
                t2 = (findex(age,tburst[ib]+dtburst[ib]))<(nage-1)
                sfhburst1[t1:t2,ib] = aburst[ib]
                mburst[ib] = aburst[ib]*(interpolate(age,t2)-interpolate(age,t1))*1D9 ; =dtburst[ib] except on the edges
             endif
          endif
; Gaussian burst
          if keyword_set(gaussburst) then begin
             sfhburst1[*,ib] = aburst[ib]*exp(-0.5*((age-tburst[ib])/dtburst[ib])^2)/sqrt(2.0*!pi) ; [Msun/yr]
             mburst[ib] = im_integral(age*1D9,sfhburst1[*,ib]) ; [Msun]
          endif
       endfor 
       sfhburst = total(sfhburst1,2,/double)
    endif else sfhburst = sfhtau*0D

    sfh = sfhtau + sfhburst ; combine

; truncate the last burst? if so then we need to "correct" MBURST; we
; do not correct MTAU, although maybe we should 
    dotruncate = 0
    ilast = -1
    if tag_exist(info,'tauburst') then begin
       if (nb gt 0) and (info.tauburst gt 0.0) then begin
          if (keyword_set(notruncate) eq 0) then begin
             dotruncate = 1
             ilast = long(findex(age,tburst[nb-1]))
             post = where(age ge tburst[nb-1],npost) ; after the burst
             if (npost gt 0) then begin
; adjust the full SFH
                sfrpostburst = interpol(sfh,age,tburst[nb-1])
                sfh[post] = sfrpostburst*exp(-(age[post]-tburst[nb-1])/info.tauburst)
; now adjust SFHBURST and SFHTAU                
                sfrpostburst1 = interpol(sfhburst1[*,nb-1],age,tburst[nb-1])
                sfhburst1[post,nb-1] = sfrpostburst1*exp(-(age[post]-tburst[nb-1])/info.tauburst)
                mburst[nb-1] = im_integral(age*1D9,sfhburst1[*,nb-1]) ; [Msun]
;               mburst[nb-1] = 0.5D*mburst[nb-1] + sfrpostburst*info.tauburst[nb-1]*1D9*$
;                 exp(-tburst[nb-1]/info.tauburst[nb-1]) ; [Msun]                
                sfhtau[post] = 0
                sfhburst = total(sfhburst1,2,/double)
             endif
          endif
       endif
    endif

; if requested, compute the <SFR> over the previous 100 Myr and the
; birthrate parameter
    if arg_present(outmgalaxy) or arg_present(outsfr100) or $
      arg_present(outb100) then begin
       dt = 0.1D ; [100 Myr]
       mgalaxy = sfh*0D
       sfr100 = sfh*0D
       b100 = sfh*0D
       for iage = 0, nage-1 do begin
; if the SFH has been truncated then do the integrals, otherwise
; calculate the burst masses with integrals and the tau-model
; analytically
          if dotruncate and (iage gt ilast) then begin
             mtot100 = im_integral(age*1D9,sfh,1D9*(age[iage]-dt)>0,1D9*age[iage])
             mgalaxy[iage] = mgalaxy[ilast] + im_integral(age*1D9,sfh,1D9*age[ilast],1D9*age[iage])
          endif else begin
             if (nb eq 0) then begin
                m100burst = 0D
                mtotburst = 0D
             endif else begin
                m100burst = im_integral(age*1D9,sfhburst,1D9*(age[iage]-dt)>0,1D9*age[iage])
                mtotburst = im_integral(age*1D9,sfhburst,0D,1D9*age[iage])
             endelse
             if (info.tau eq 0.0) then begin
                mtot100 = 0D + m100burst
                mgalaxy[iage] = mtau + mtotburst
             endif else begin
                mtot100 = mtau*(exp(-((age[iage]-dt)>0)/info.tau)-exp(-age[iage]/info.tau)) + m100burst
                mgalaxy[iage] = mtau*(1D0-exp(-age[iage]/info.tau)) + mtotburst
             endelse
          endelse
          sfr100[iage] = mtot100/(dt*1D9)
          b100[iage] = sfr100[iage]/(mgalaxy[iage]/(1D9*age[iage]))

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
    outsfhtau = interpolate(sfhtau,findx)
    outsfhburst = interpolate(sfhburst,findx)
    if arg_present(outmgalaxy) then outmgalaxy = interpolate(mgalaxy,findx)
    if arg_present(outsfr100) then outsfr100 = interpolate(sfr100,findx)
    if arg_present(outb100) then outb100 = interpolate(b100,findx)

; QAplot    
    if keyword_set(debug) then begin
       djs_plot, age, sfh, xlog=0, xsty=3, ysty=3, psym=-6, _extra=extra;, xr=[min(age)>0.01,max(age)]
       djs_oplot, age, sfr100, psym=6, color='cyan'
       djs_oplot, outage, outsfh, psym=6, color='orange'
       djs_oplot, outage, outsfr100, psym=6, color='blue'
;      djs_oplot, outage, outsfhtau, color='blue', psym=-6, sym=0.5
;      djs_oplot, outage, outsfhburst, color='red', psym=-6, sym=0.5
       if dotruncate then djs_oplot, tburst[nb-1]+info.tauburst*[1,1], !y.crange, color='yellow'
;      for ib = 0, nb-1 do djs_oplot, tburst[ib]*[1,1], !y.crange, color='yellow'
    endif

return, outsfh
end
