function isedfit_reconstruct_sfh, info, age=age, mtau=mtau, $
  aburst=aburst, mburst=mburst, mgalaxy=mgalaxy, sfr100=sfr100, $
  b100=b100, notruncate=notruncate, sfhtau=sfhtau, sfhburst=sfhburst, $
  debug=debug, _extra=extra
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
    
    if (n_elements(age) eq 0) then age = build_isedfit_agegrid(info,$
      debug=0,_extra=extra)
    nage = n_elements(age)
    
; deal with the simple tau model and with bursty SFHs    
    if (n_elements(mtau) eq 0) then mtau = 1D ; normalization
    if (info.tau eq 0.0) then sfhtau = dblarr(nage) else $
      sfhtau = mtau*exp(-age/info.tau)/(info.tau*1D9) 
    
    nb = info.nburst
    if (nb gt 0) then begin
       sfhburst1 = reform(dblarr(nage,nb),nage,nb)
       fburst = 1D*info.fburst[0:nb-1]
       tburst = 1D*info.tburst[0:nb-1]
       dtburst = 1D*info.dtburst[0:nb-1]
; compute the amplitude, and then build each burst as a Gaussian
       aburst = dblarr(nb)
       mburst = dblarr(nb)
       for ib = 0, nb-1 do begin
          if (info.tau eq 0.0) then $
            aburst[ib] = fburst[ib]*mtau/(dtburst[ib]*1D9) else $
              aburst[ib] = fburst[ib]*mtau*(1.0-exp(-tburst[ib]/info.tau))/(dtburst[ib]*1D9)
;         if (info.tau eq 0.0) then $
;           aburst[ib] = fburst[ib]*mtau/(sqrt(2.0*!pi)*dtburst[ib]*1D9) else $
;             aburst[ib] = fburst[ib]*mtau*(1.0-exp(-tburst[ib]/info.tau))/$
;           (sqrt(2.0*!pi)*dtburst[ib]*1D9)

          sfhburst1[*,ib] = aburst[ib]*exp(-0.5*((age-tburst[ib])/dtburst[ib])^2)/sqrt(2.0*!pi) ; [Msun/yr]
;         sfhburst1[*,ib] = aburst[ib]*exp(-0.5*((age-tburst[ib])/dtburst[ib])^2) ; [Msun/yr]
          if arg_present(mburst) then mburst[ib] = im_integral(age*1D9,sfhburst1[*,ib]) ; [Msun]
;         mburst[ib] = aburst[ib]*sqrt(2.0*!pi)*dtburst[ib]*1D9 ; [Msun]
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
                if arg_present(mburst) then mburst[nb-1] = im_integral(age*1D9,sfhburst1[*,nb-1]) ; [Msun]
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
    if arg_present(mgalaxy) or arg_present(sfr100) or arg_present(b100) then begin
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
       endfor
    endif
    
; QAplot    
    if keyword_set(debug) then begin
       djs_plot, age, sfh, /xlog, xsty=3, ysty=3, psym=-6;, xr=[4.0,9]
       djs_oplot, age, sfhtau, color='blue', psym=-6, sym=0.5
       djs_oplot, age, sfhburst, color='red', psym=-6, sym=0.5
       if dotruncate then djs_oplot, tburst[nb-1]+info.tauburst*[1,1], !y.crange, color='yellow'
;      for ib = 0, nb-1 do djs_oplot, tburst[ib]*[1,1], !y.crange, color='yellow'
    endif

return, sfh
end

;; renormalize such that the integrated SFH (from 0-->infinity is 1_Msun)
;    if (nb eq 0) then mgal = mtau else begin ; by definition 
;       mtotburst = total(mburst,/double)
;       mgal = mtau + mtotburst
;       mfracburst = mburst/mgal
;       sfhburst = total(mfracburst)*sfhburst
;    endelse
;
;    mfractau = mtau/mgal
;    sfhtau = mfractau*sfhtau

