;+
; NAME:
;   GALEX_STARFINDER()
;
; PURPOSE:
;   Simple wrapper script to run STARFINDER on GALEX images using a
;   standard set of parameters and return the point-source catalog.
;
; INPUTS: 
;   img - input (-int) image [NX,NY]
;   sky - sky background map (-skybg) [NX,NY]
;   psf - point spread function [NXPSF,NYPSF]
;
; OPTIONAL INPUTS: 
;   hdr - FITS header; if present, used to convert the pixel
;     coordinates of each source to celestial coordinates
;   threshold - *absolute* detection threshold; the default is to
;     compute the threshold from the sky-subtracted image
;   zpt - magnitude zeropoint (e.g., zpt_NUV=20.08, zpt_FUV=18.82) 
;   min_correlation, no_slant, four, deblend, deblost - see STARFINDER
;     documentation  
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   catalog - output catalog with the following tags [NSTARS]
;     id - unique identification number
;     [x,y]image - [x,y] centroid
;     flux - total flux (integrated to infinity) [count/s]
;     mag - instrumental magnitude if ZPT is not defined
;       (=-2.5*alog10(FLUX)), otherwise calibrated magnitude 
;     correlation - see STARFINDER documentation
;     ra,dec - celestial coordinates (if HDR is provided) 
;
; OPTIONAL OUTPUTS:
;   modelimg - model image (see STARFINDER documentation for STARS) 
;   params - input/output parameters
;
; COMMENTS:
;   The threshold is computed from the standard deviation of the
;   background-subtracted image using the MMM algorithm, which is
;   optimized for crowded stellar fields.
;   
;   ToDo: interpolate at the exposure-time map and the flag map. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Apr 25, UCSD
;
; Copyright (C) 2010, John Moustakas
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

function galex_starfinder, img, sky, psf, hdr=hdr, threshold=threshold, $
  zpt=zpt, min_correlation=min_correlation, no_slant=no_slant, four=four, $
  deblend=deblend, deblost=deblost, modelimg=modelimg, params=params

    if (n_elements(img) eq 0) or (n_elements(sky) eq 0) or $
      (n_elements(psf) eq 0) then begin
       doc_library, 'galex_starfinder'
       return, -1
    endif

; unless given, compute the threshold for source detection
    if (n_elements(threshold) eq 0) then begin
       splog, 'Computing image sigma'
       good = where(img gt 0)
       mmm, img[good]-sky[good], mode, sigma
       threshold = [1.0]*sigma
;      threshold = [1.0,1.0]*sigma
    endif else begin
       sigma = -1.0
       mode = -1.0
    endelse
    if (n_elements(min_correlation) eq 0) then min_correlation = 0.7
    if (n_elements(zpt) eq 0) then zpt = 0.0
    if (n_elements(no_slant) eq 0) then no_slant = 1 ; always recommended!
    if (n_elements(deblend) eq 0) then deblend = 1   ; always recommended!
;   if (n_elements(deblost) eq 0) then deblost = 1   ; always recommended!

; pack the parameters into a structure
    params = {$
      sigma:           sigma,$
      mode:            mode,$
      threshold:       threshold,$
      min_correlation: min_correlation,$
      zpt:             zpt,$
      no_slant:        keyword_set(no_slant),$
      four:            keyword_set(four),$
      deblend:         keyword_set(deblend),$
      deblost:         keyword_set(deblost),$
      time:            0.0} ; [minutes]

; call starfinder!    
    t0 = systime(1)
    starfinder, img, psf, threshold, min_correlation, $
      xx, yy, ff, xx_err, yy_err, ff_err, cc, background=sky, $
      no_slant=no_slant, four=four, deblend=deblend, $
      deblost=deblost, stars=modelimg
    splog, 'Total time (minutes) = ', (systime(1)-t0)/60.0
    params.time = (systime(1)-t0)/60.0

; pack into a catalog and return
    catalog = {id: -1L, ximage: -1.0, yimage: -1.0, ra: -1.0D, $
      dec: -1.0D, flux: -1.0, mag: -1.0, correlation: -1.0}

    nobj = n_elements(xx)
    if (nobj gt 0L) then begin
       catalog = replicate(catalog,nobj)
       catalog.id = lindgen(nobj)
       catalog.ximage = xx
       catalog.yimage = yy
       catalog.flux = ff
       catalog.mag = -2.5*alog10(ff)+zpt
       catalog.correlation = cc
       if (n_elements(hdr) ne 0) then begin
          extast, hdr, astr
          xy2ad, xx, yy, astr, ra, dec
          catalog.ra = ra
          catalog.dec = dec
       endif
    endif

return, catalog
end
