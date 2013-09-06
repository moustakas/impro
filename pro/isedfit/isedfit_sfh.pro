;+
; NAME:
;   ISEDFIT_SFH()
;
; PURPOSE:
;   Generate an iSEDfit-compatible star formation history including
;   allowing for multiple bursts.
;
; INPUTS: 
;   sfhinfo - iSEDfit-style star formation history structure with the
;     following tags (note, only .TAU is mandatory, but see also the
;     TAU and MGAL optional inputs):
;       .TAU - characteristic e-folding time [Gyr]
;       .TRUNCTAU - characteristic e-folding time over which the last
;         burst should be truncated [Gyr]
;       .NBURST - number of bursts
;       .TBURST - time/age at which each burst begins [Gyr]
;       .DTBURST - burst duration (exact meaning depends on BURSTTTYPE) [Gyr]
;       .FBURST - total stellar mass fraction formed in each burst
;       .TOTALMASS - if this tag is present, then normalize everything
;         (SFH, SFR100, MGALAXY, etc.) to this value at the final
;         age/time bin (i.e., the last element of OUTAGE); this
;         parameter corresponds to the total (time-integrated) stellar
;         mass formed (ignoring mass lost through winds, SNe, etc.)
;         over the age of the "galaxy" (the default is to normalize to
;         the  value of MGALAXY at the last time/age bin) 
;   tau - in lieu of passing SFHINFO, a basic star-formation history
;     (no bursts, truncated or otherwise) be constructed by simply
;     passing the characteristic e-folding time through this variable
;     [Gyr] 
;   totalmass - drop-in replacement for SFHINFO.TOTALMASS
; 
; OPTIONAL INPUTS: 
;   nagemax - maximum number of age/time points to use internally
;     (default 250); probably shouldn't change this!
;   dage - minimum time resolution at which the SFH should be sampled
;     internally (default 20 Myr); probably shouldn't change
;     this! [Gyr]
;   maxage - compute the SFH between 0 and MAXAGE [Gyr] (this optional
;     input is trumped by OUTAGE)
;   outage - time/age corresponding to SFH [Gyr, NOUTAGE]
;   useage - optionally force this routine to internally use this
;     age/time vector [NAGE, GYR]; not generally recommended!
;   mtau - normalization of the TAU model after infinite time (default
;     1.0) [Msun]  
;   
; KEYWORD PARAMETERS:
;   delayed - construct a "delayed" rather a "simple" tau model (see
;     iSEDfit documentation) (recommended)
;   bursttype - type of burst (see iSEDfit documentation)
;     0 - step function
;     1 - Gaussian (default)
;     2 - step function with exponential wings
;   notruncate - override any truncated bursts specified in SFHINFO
;   linear - use a linearly spaced age/time vector in ISEDFIT_AGEGRID
;     (default is to use logarithmic spacing) (not generally
;     recommended!) 
;   debug - make a very simple debugging plot
;   _extra - additional plotting keywords (only for /DEBUG) 
; 
; OUTPUTS: 
;   sfh - star formation rate as a function of time/age [Msun/yr, NOUTAGE] 
;
; OPTIONAL OUTPUTS:
;   sfhtau - underlying smooth star formation history (i.e., just the
;     TAU component of SFH) [Msun/yr, NOUTAGE]
;   sfhburst - bursty star formation history (i.e., just the bursty
;     parts of SFH) [Msun/yr, NOUTAGE]
;   aburst - burst amplitude [Msun/yr, NBURST]
;   mburst - total mass formed in each burst [Msun, NBURST]
;   mgalaxy - total mass in stars formed (i.e., galaxy mass, ignoring
;     mass loss) [Msun, NOUTAGE]
;   sfr100 - SFR averaged over the previous 100 Myr [Msun/yr, NOUTAGE] 
;   b100 - birthrate parameter averaged over the previous 100 Myr [NOUTAGE] 
;   b_1000 - birthrate parameter averaged over the previous 1000 Myr
;     (note the underscore, needed to circumvent IDL's keyword
;     rules!) [NOUTAGE]  
;   sfrage - SFR-weighted age [Gyr, NOUTAGE]
;
; COMMENTS:
;   If a truncated burst is requested then the SFHTAU and SFHBURST
;   outputs are not modified.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Dec 22, UCSD - developed
;   jm13aug12siena - updated to the latest data model; fully
;     documented 
;
; Copyright (C) 2011, 2013, John Moustakas
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

function isedfit_sfh, sfhinfo, tau=tau, outage=outage, totalmass=totalmass, $
  mtau=mtau, useage=useage, nagemax=nagemax, dage=dage, maxage=maxage, $
  sfhtau=outsfhtau, sfhburst=outsfhburst, aburst=aburst, mburst=mburst, $
  mgalaxy=outmgalaxy, sfr100=outsfr100, b100=outb100, b_1000=outb1000, $
  sfrage=outsfrage, delayed=delayed, bursttype=bursttype, notruncate=notruncate, $
  linear=linear, debug=debug, _extra=extra
    
    if n_elements(sfhinfo) eq 0 and n_elements(tau) eq 0 then begin
       doc_library, 'isedfit_sfh'
       return, -1
    endif

; construct the variables we will need
    if n_elements(sfhinfo) eq 0 then begin
       if n_elements(tau) ne 1 then begin
          splog, 'Either pass a scalar (positive) TAU value or SFHINFO!'
          return, -1
       endif
       nb = 0
       trunctau = 0D
    endif else begin
       if tag_exist(sfhinfo,'TAU') eq 0 then begin
          splog, 'SFHINFO structure requires TAU tag!'
          return, -1
       endif else tau = sfhinfo.tau
       if tag_exist(sfhinfo,'NBURST') then nb = sfhinfo.nburst else nb = 0
       if (nb gt 0) then begin
          if tag_exist(sfhinfo,'FBURST') eq 0 or tag_exist(sfhinfo,'FBURST') eq 0 or $
            tag_exist(sfhinfo,'FBURST') eq 0 then begin
             splog, 'SFHINFO structure requires FBURST, TBURST, and DTBURST tags!'
             return, -1
          endif
          fburst = im_double(sfhinfo.fburst[0:nb-1])
          tburst = im_double(sfhinfo.tburst[0:nb-1])
          dtburst = im_double(sfhinfo.dtburst[0:nb-1])
       endif
       if tag_exist(sfhinfo,'TRUNCTAU') then $
         trunctau = sfhinfo.trunctau else trunctau = 0.0D
       if tag_exist(sfhinfo,'TOTALMASS') then totalmass = sfhinfo.totalmass
    endelse
    
; need a highly sampled age grid to get the integrations right; this
; grid is used internally and then the output quantities are
; interpolated at OUTAGE
    if n_elements(nagemin) eq 0 then nagemin = 10
    if n_elements(nagemax) eq 0 then nagemax = 250
    if n_elements(dage) eq 0 then dage = 0.02D ; default 20 Myr time interval
    if n_elements(bursttype) eq 0 then bursttype = 1 ; Gaussian burst

    minage = 0D ; never change this!
    if n_elements(maxage) eq 0 then begin
       if (n_elements(outage) eq 0) then maxage = 13.8D else $
         maxage = im_double(max(outage))
    endif

    nage = long((((maxage-minage)/dage)<nagemax)>nagemin)
    if (n_elements(useage) eq 0) then begin
       age = isedfit_agegrid(sfhinfo,tau=tau,debug=0,nage=nage,$
         minage=minage,maxage=maxage,linear=linear,$
         delayed=delayed,bursttype=bursttype)
    endif else age = useage

    if (n_elements(outage) eq 0) then outage = age
    nage = n_elements(age)

; deal with the simple tau model, including delayed tau models, and
; with bursty SFHs     
    if (n_elements(mtau) eq 0) then mtau = 1D ; normalization [Msun]
    if (tau eq 0D) then sfhtau = dblarr(nage) else begin
       if keyword_set(delayed) then $
         sfhtau = mtau*age*1D9*exp(-age/tau)/(tau*1D9)^2 else $
           sfhtau = mtau*exp(-age/tau)/(tau*1D9)
    endelse

; type of burst: 0=step function (default); 1=gaussian; 2=step
; function with exponential wings 
    if (nb gt 0) then begin
       sfhburst1 = reform(dblarr(nage,nb),nage,nb)
; compute the amplitude, and then build each burst in turn
       aburst = dblarr(nb)
       mburst = dblarr(nb)
       for ib = 0, nb-1 do begin
          if (tau eq 0D) then $
            aburst[ib] = fburst[ib]*mtau/(dtburst[ib]*1D9) else $
              aburst[ib] = fburst[ib]*mtau*(1.0-exp(-tburst[ib]/tau))/(dtburst[ib]*1D9)
; step-function burst (default)
          if bursttype eq 0 then begin
             if (max(age) ge tburst[ib]) then begin
                t1 = (findex(age,tburst[ib]))>0
                t2 = (findex(age,tburst[ib]+dtburst[ib]))<(nage-1)
                sfhburst1[t1:t2,ib] = aburst[ib]
; this should equal dtburst[ib] except on the edges
                mburst[ib] = aburst[ib]*(interpolate(age,t2)-interpolate(age,t1))*1D9 
             endif
          endif
; Gaussian burst; could do the integral analytically, but no need to
          if bursttype eq 1 then begin
             sfhburst1[*,ib] = aburst[ib]*exp(-0.5*((age-tburst[ib])/$
               dtburst[ib])^2)/sqrt(2.0*!pi)                   ; [Msun/yr]
             mburst[ib] = im_integral(age*1D9,sfhburst1[*,ib]) ; [Msun]
          endif
; step-function burst with exponential wings; the MBURST() integral
; won't be quite right because of the discontinuity
          if bursttype eq 2 then begin
             during = where((age ge tburst[ib]) and (age le tburst[ib]+dtburst[ib]),nduring)
             before = where(age lt tburst[ib],nbefore)
             after = where(age gt tburst[ib]+dtburst[ib],nafter)
             if (nbefore ne 0) then sfhburst1[before,ib] += $
               aburst[ib]*exp(-(tburst[ib]-age[before])/0.01D)
             if (nafter ne 0) then sfhburst1[after,ib] += $
               aburst[ib]*exp(-(age[after]-(tburst[ib]+dtburst[ib]))/0.01D)
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
    if (nb gt 0) and (trunctau gt 0D) then begin
       if (keyword_set(notruncate) eq 0) then begin
          dotruncate = 1
          ilast = (long(findex(age,tburst[nb-1])))>0
          if ilast lt 0 then message, 'This should not happen'
          if bursttype eq 1 then $ ; Gaussian burst shape
            post = where(age ge tburst[nb-1],npost) else $ ; after the peak of the burst
              post = where(age ge tburst[nb-1]+dtburst[nb-1],npost) ; after the full width of the burst
          if (npost gt 0) then begin
; truncate the full SFH
             sfrpostburst = interpol(sfh,age,tburst[nb-1])
             sfh[post] = sfrpostburst*exp(-(age[post]-tburst[nb-1])/trunctau)

; code below truncates SFHBURST and SFHTAU, and adjusts MBURST for the
; truncation; in general we don't use these quantities, so skip
; the extra work              
;             sfrpostburst1 = interpol(sfhburst1[*,nb-1],age,tburst[nb-1])
;             sfhburst1[post,nb-1] = sfrpostburst1*exp(-(age[post]-tburst[nb-1])/trunctau)
;             mburst[nb-1] = im_integral(age*1D9,sfhburst1[*,nb-1]) ; [Msun]
;;            mburst[nb-1] = 0.5D*mburst[nb-1] + sfrpostburst*trunctau[nb-1]*1D9*$
;;              exp(-tburst[nb-1]/trunctau[nb-1]) ; [Msun]                
;             sfhtau[post] = 0
;             sfhburst = total(sfhburst1,2,/double)
          endif
       endif
    endif

; if requested, compute the <SFR> over the previous 100 Myr, the
; birthrate parameter, and the SFH-weighted age
;   t0 = systime(1)
    if arg_present(outmgalaxy) or arg_present(outsfr100) or $
      arg_present(outb100) or arg_present(outb1000) or $
      arg_present(outsfrage) then begin
       dt100 = 0.1D ; [100 Myr]
       dt1000 = 1D ; [1000 Myr = 1 Gyr]
       mgalaxy = sfh*0D
       sfr100 = sfh*0D
       sfr1000 = sfh*0D
       b100 = sfh*0D
       b1000 = sfh*0D
       sfrage = sfh*0D
       for iage = 0L, nage-1 do begin
; if the SFH has been truncated then do the integrals, otherwise
; calculate the burst masses with integrals and the tau-model
; analytically
          if dotruncate and (iage gt ilast) then begin
             mtot100 = im_integral(age*1D9,sfh,1D9*(age[iage]-dt100)>0,1D9*age[iage])
             mtot1000 = im_integral(age*1D9,sfh,1D9*(age[iage]-dt1000)>0,1D9*age[iage])
             mgalaxy[iage] = mgalaxy[ilast] + im_integral(age*1D9,sfh,1D9*age[ilast],1D9*age[iage])
          endif else begin
             if (nb eq 0) then begin
                m100burst = 0D
                m1000burst = 0D
                mtotburst = 0D
             endif else begin
                m100burst = im_integral(age*1D9,sfhburst,1D9*(age[iage]-dt100)>0,1D9*age[iage])
                m1000burst = im_integral(age*1D9,sfhburst,1D9*(age[iage]-dt1000)>0,1D9*age[iage])
                mtotburst = im_integral(age*1D9,sfhburst,0D,1D9*age[iage])
             endelse
             if (tau eq 0D) then begin
                mtot100 = 0D + m100burst
                mtot1000 = 0D + m1000burst
                mgalaxy[iage] = mtau + mtotburst
             endif else begin
                mtot100 = mtau*(exp(-((age[iage]-dt100)>0)/tau)-exp(-age[iage]/tau)) + m100burst
                mtot1000 = mtau*(exp(-((age[iage]-dt1000)>0)/tau)-exp(-age[iage]/tau)) + m1000burst
                mgalaxy[iage] = mtau*(1D0-exp(-age[iage]/tau)) + mtotburst
             endelse
          endelse
          if (age[iage] gt 0D) then begin
             sfr100[iage] = mtot100/(dt100*1D9)
             sfr1000[iage] = mtot1000/(dt1000*1D9)
             b100[iage] = sfr100[iage]/(mgalaxy[iage]/(1D9*age[iage])) ; see Brinchmann+04
             b1000[iage] = sfr1000[iage]/(mgalaxy[iage]/(1D9*age[iage]))
          endif
; compute the SFR/mass-weighted age
          if (tau eq 0D) then sfrage[iage] = age[iage] else begin
             norm = im_integral(age*1D9,sfh,0D,1D9*age[iage])
             if (norm eq 0D) then sfrage[iage] = age[iage] else $
               sfrage[iage] = im_integral(age*1D9,sfh*(age[iage]-age)*1D9,0D,1D9*age[iage])/norm
          endelse
       endfor
    endif
;   splog, systime(1)-t0

; normalize?
    if n_elements(totalmass) ne 0 then begin
       if n_elements(mgalaxy) ne 0L then mgal = mgalaxy[nage-1L] else $
         mgal = im_integral(age*1D9,sfh,0D,1D9*age[nage-1L])
       if mgal eq 0.0 then mgal = 1.0 ; special case for tau=0
       massnorm = totalmass/mgal
       sfh *= massnorm
       sfhtau *= massnorm
       sfhburst *= massnorm
       if n_elements(mgalaxy) ne 0L then mgalaxy *= massnorm
       if n_elements(sfr100) ne 0L then sfr100 *= massnorm
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
    if arg_present(outb1000) then outb1000 = interpolate(b1000,findx)
    if arg_present(outsfrage) then outsfrage = interpolate(sfrage,findx)/1D9 ; [Gyr]
    
; QAplot    
    if keyword_set(debug) then begin
       djs_plot, age, sfh, xlog=0, xsty=3, ysty=3, psym=-6, _extra=extra;, xr=[min(age)>0.01,max(age)]
       djs_oplot, [outage], [outsfh], psym=6, color='orange'
       if dotruncate then djs_oplot, tburst[nb-1]+trunctau*[1,1], !y.crange, color='yellow'
;      for ib = 0, nb-1 do djs_oplot, tburst[ib]*[1,1], !y.crange, color='yellow'
    endif

return, outsfh
end
