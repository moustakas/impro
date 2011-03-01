pro tcorrelate
; jm01jul19uofa
; test out cross-correlation using a couple different methods for my
; edification  
    
; generate two gaussians, separated by 40 pixels.  find the maximum in
; the cross-correlation using IDL's C_CORRELATE (which requires two
; functions with the same number of elements), and Schlegel's
; DJS_CORRELATE, which does not.
    
; note that for both these routines the lags are *negative* and
; correspond to the amount by which you need to shift *x* to overlap
; with *y*

; test #1
; ----------------------------------------------------------------------
    
    npix = findgen(101)
    gausspix1d, npix, [50.0,1.0,5.0], gauss1
    gausspix1d, npix, [10.0,1.0,5.0], gauss2
    
    lags = -lindgen(101)
    result = c_correlate(gauss1,gauss2,lags,/double)

    result = djs_correlate(gauss1,gauss2,lags)
    
    plot, lags, result ; the maximum should be at -40

; ----------------------------------------------------------------------
; test #2
    
    npix1 = findgen(101)
    gausspix1d, npix1, [50.0,1.0,5.0], gauss1
    npix2 = findgen(21)
    gausspix1d, npix, [10.0,1.0,5.0], gauss2
    
    lags = -lindgen(101)
    result = djs_correlate(gauss1,gauss2,lags)

    plot, lags, result ; the maximum should be at -40

return
end
