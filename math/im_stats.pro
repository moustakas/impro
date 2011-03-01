function im_stats, array, verbose=verbose, baremin=baremin, mask=mask, $
  sigrej=sigrej, _extra=extra
; jm03jun2uofa
; jm04sep07uofa - added mask output
; jm05jan20uofa - added MAD output, removed VARIANCE, SKEWNESS, and
;                 KURTOSIS 
; jm05jun14uofa - added BAREMIN; only print the bare minimum
;                 statistics
; jm05aug31uofa - added SIGREJ optional input    

    if (n_elements(sigrej) eq 0L) then sigrej = 4.0
    
    narray = n_elements(array)
    finite = where(finite(array),nfinite)
    
    farray = array[finite]

    djs_iterstat, farray, mean=amean_rej, median=amedian_rej, $
      sigma=asigma_rej, mask=mask, sigrej=sigrej, _extra=extra

    good = where(mask,ngood)
    
    amin = min(array,/nan,_extra=extra)
    amax = max(array,/nan,_extra=extra)
    amean = mean(array,/nan,_extra=extra)
    amedian = median(array,_extra=extra)
    asigma = stddev(array,/nan)

    amad = total(abs(farray-amedian))/float(nfinite)
    
    if (ngood ne 0L) then begin
       aminrej = min(array[good],/nan)
       amaxrej = max(array[good],/nan)
    endif else begin
       aminrej = amin
       amaxrej = amax
    endelse

; compute the lower and upper intervals that include 68% and 95% of
; the data    

    histarray = farray[sort(farray)]
    
    linterp, findgen(nfinite)/float(nfinite), histarray, [0.025,0.16,0.84,0.975], dist

    sig68lo = dist[1]
    sig68up = dist[2]
    sig95lo = dist[0]
    sig95up = dist[3]

    sig68mean = (dist[2]-dist[1]) / 2.0 ; mean 68% interval
    
    stats = {$
      min:        amin,        $
      max:        amax,        $
      mean:       amean,       $
      median:     amedian,     $
      sigma:      asigma,      $
      mad:        amad,        $
      npts:       nfinite,     $
      minrej:     aminrej,     $ ; minimum after iterative rejection
      maxrej:     amaxrej,     $ ; maximum after iterative rejection
      mean_rej:   amean_rej,   $
      median_rej: amedian,     $
      sigma_rej:  asigma_rej,  $
      npts_rej:   ngood,       $
      sig68lo:    sig68lo, $
      sig68up:    sig68up, $
      sig95lo:    sig95lo, $
      sig95up:    sig95up, $
      sig68mean:  sig68mean}
    
    if keyword_set(verbose) then begin
       if keyword_set(baremin) then begin
          struct_print, struct_trimtags(stats,select=['MINREJ','MAXREJ','MEAN_REJ','MEDIAN_REJ','SIGMA_REJ','NPTS']), _extra=extra
;         struct_print, struct_trimtags(stats,select=['MIN','MAX','MEAN','MEDIAN','SIGMA','NPTS']), _extra=extra
;         help, struct_trimtags(stats,select=['MEDIAN','MEAN_REJ','SIG68MEAN']), /structure x[
       endif else begin
          help, stats, /structure
       endelse
    endif
;   if keyword_set(verbose) then struct_print, stats, _extra=extra
    
return, stats
end    
