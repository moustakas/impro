;+
; NAME:
;   IM_STATS()
;
; PURPOSE:
;   Compute some basic statistics and pack into a structure.
;
; INPUTS: 
;   array - input vector
;
; OPTIONAL INPUTS: 
;   sigrej - rejection threshold (default 4)
;   extra - extra keywords for DJS_ITERSTAT
;
; KEYWORD PARAMETERS: 
;   verbose - print to the screen
;   baremin - only print the bare minimum statistics 
;
; OUTPUTS: 
;   stats - output data structure with lots of goodies
;
; OPTIONAL OUTPUTS:
;   mask - bit mask of rejected values
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 Jun 02, U of A
;   jm04sep07uofa - added mask output
;   jm05jan20uofa - added MAD output, removed VARIANCE, SKEWNESS, and
;     KURTOSIS  
;   jm05jun14uofa - added /BAREMIN
;   jm05aug31uofa - added SIGREJ optional input    
;
; Copyright (C) 2003-2005, John Moustakas
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

function im_stats, array, verbose=verbose, baremin=baremin, $
  mask=mask, sigrej=sigrej, _extra=extra

    if (n_elements(array) eq 0L) then begin
       doc_library, 'im_stats'
       return, -1
    endif
    
    if (n_elements(sigrej) eq 0L) then sigrej = 4.0
    
    narray = n_elements(array)
    finite = where(finite(array),nfinite)
    
    farray = array[finite]

    amin = min(farray,_extra=extra)
    amax = max(farray,_extra=extra)
    amean = mean(farray,_extra=extra)
    amedian = djs_median(farray,_extra=extra)
    asigma = djsig(farray)
    amad = total(abs(farray-amedian))/float(nfinite)

    if narray eq 1 then begin
       amean_rej = amean
       amedian_rej = amedian
       asigma_rej = asigma
       good = 0 & ngood = 1
    endif else begin
       djs_iterstat, farray, mean=amean_rej, median=amedian_rej, $
         sigma=asigma_rej, mask=mask, sigrej=sigrej, _extra=extra
       good = where(mask,ngood)
    endelse
    
    
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
    if narray gt 1 then begin
       linterp, findgen(nfinite)/float(nfinite), histarray, [0.025,0.16,0.84,0.975], dist
       sig68lo = dist[1]
       sig68up = dist[2]
       sig95lo = dist[0]
       sig95up = dist[3]
       sig68mean = (dist[2]-dist[1]) / 2.0 ; mean 68% interval
    endif else begin
       sig68lo = 0.0
       sig68up = 0.0
       sig95lo = 0.0
       sig95up = 0.0
       sig68mean = 0.0
    endelse
    
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
          struct_print, struct_trimtags(stats,select=['MINREJ','MAXREJ',$
            'MEAN_REJ','MEDIAN_REJ','SIGMA_REJ','NPTS']), _extra=extra
       endif else begin
          help, stats, /structure
       endelse
    endif
;   if keyword_set(verbose) then struct_print, stats, _extra=extra
    
return, stats
end    
