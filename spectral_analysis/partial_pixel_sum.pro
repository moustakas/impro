function partial_pixel_sum, wave, flux, ferr, lowave, hiwave, $
  sum_err=sum_err, indxlohi=indxlohi, nindxlohi=nindxlohi, $
  partial_flag=partial_flag
; jm03oct30uofa
; jm05jul28uofa - added NINDXLOHI and PARTIAL_FLAG optional outputs 
; needs error checking

;  splog, 'THIS ROUTINE NEEDS A REWRITE!!!'

    partial_flag = 0L
    
    flux = double(flux)
    npix = n_elements(wave)

    var = ferr^2.0 ; variance

; find the pixel indices corresponding to the min and max wavelengths 

    get_element, wave, [lowave,hiwave], indxlohi

    if wave[indxlohi[0]] gt lowave then begin
       if ((indxlohi[0]-1L) lt 0L) then partial_flag = 1L
       indxlohi[0] = (indxlohi[0]-1L)>0L
    endif

    if wave[indxlohi[1]] lt hiwave then begin
       if ((indxlohi[1]+1L) gt (npix-1L)) then partial_flag = 1L
       indxlohi[1] = (indxlohi[1]+1L)<(npix-1L)
    endif

    if (indxlohi[0] eq indxlohi[1]) then begin ; return
       sum_err = 0.0
       return, -2.0
    endif

    nindxlohi = (indxlohi[1]-indxlohi[0])+1L
    pixels = lindgen(nindxlohi)+indxlohi[0]
    
    nsamp = 10.0
    dwave = (wave[indxlohi[1]]-wave[indxlohi[0]])/(nindxlohi-1) ; mean dispersion [A/pixel]

    bigwave = findgen((nindxlohi-1.0)*nsamp+1)*dwave/nsamp+wave[pixels[0]]
    bigflux = interpol(flux,wave,bigwave)
    bigvar = interpol(var,wave,bigwave)

; compute the pixel weights; the pixel weight on the endpoints is one
; minus the pixel fraction to be included
    
;   weights = replicate(1.0D,nindxlohi)
;   weights[0L] = 1.0-(lowave-wave[indxlohi[0]])/(wave[indxlohi[0]+1L]-wave[indxlohi[0]])
;   weights[nindxlohi-1L] = 1.0-(wave[indxlohi[1]]-hiwave)/(wave[indxlohi[1]]-wave[indxlohi[1]-1L])

;   print, lowave, wave[indxlohi[0]], lowave-wave[indxlohi[0]], weights[0]
;   print, hiwave, wave[indxlohi[1]], wave[indxlohi[1]]-hiwave, weights[nindxlohi-1]

    sum = im_integral(bigwave,bigflux,lowave,hiwave) / (hiwave-lowave)
;   sum = im_integral(wave,flux,lowave,hiwave) / (hiwave-lowave)
;   sum = tsum(wave[pixels],weights*flux[pixels]) / (hiwave-lowave)
;   sum = int_tabulated(wave[pixels],weights*flux[pixels],/double) / (hiwave-lowave)

    if arg_present(sum_err) then $
      sum_err = sqrt(im_integral(bigwave,bigvar,lowave,hiwave)) / (hiwave-lowave)
;     sum_err = (sqrt(im_integral(wave,var,lowave,hiwave)))[0] / (hiwave-lowave)
;     sum_err = sqrt(tsum(wave[pixels],weights*var[pixels]))
;     sum_err = sqrt(int_tabulated(wave[pixels],weights*var[pixels],/double))

return, sum
end
