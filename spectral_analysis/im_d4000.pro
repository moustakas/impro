function im_d4000, wave, flux, ferr=ferr, d4000_err=d4000_err, $
  bruzual83=bruzual83
; jm09feb27nyu - compute D(4000) [code stolen from SPECTRAL_INDICES];
; compute the Bruzual+83 definition if if /BRUZUAL83, otherwise use
; BALOGH+99 

    if (n_elements(ferr) eq 0L) then ferr = flux*0.05
    light = 2.99792458D18
    fnu = wave*wave*flux/light
    fnu_err = wave*wave*ferr/light

    if keyword_set(bruzual83) then begin
       wbvec = [3750.0,3950.0]
       wrvec = [4050.0,4250.0]
    endif else begin
       wbvec = [3850.0,3950.0]
       wrvec = [4000.0,4100.0]
    endelse

    llimit = total(wbvec)/2.0
    ulimit = total(wrvec)/2.0
    lwidth = wbvec[1]-wbvec[0]
    uwidth = wrvec[1]-wrvec[0]
    midwave = total([llimit,ulimit])/2.0 ; central wavelength

    cbreak = iabslineew(wave,fnu,midwave,ferr=fnu_err,llimit=llimit,$
      lwidth=lwidth,ulimit=ulimit,uwidth=uwidth,label=blabel,/noline,$
      absplot=cbreakplot,debug=0,/fnu,silent=silent,_extra=extra)

    if (cbreak.cratio gt 0.0) then begin
       d4000 = cbreak.cratio
       d4000_err = cbreak.cratio_err
    endif else begin
       d4000 = -1.0
       d4000_err = -1.0
    endelse

    if (n_elements(ferr) eq 0) then d4000_err = -2.0    
    
return, d4000
end
    
