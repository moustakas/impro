;+
; NAME:
;       IM_NORMALIZE()
;
; PURPOSE:
;       Normalize a spectrum or data array.
;
; CALLING SEQUENCE:
;       normspec = im_normalize(spec,wave,[normwave=,binsize=], $
;          /mean,/median,/max,[const=])
;
; INPUTS:
;       spec - one-dimensional spectrum or flux array
;       wave - corresponding wavelength vector or index number
;
; OPTIONAL INPUTS:
;       normwave - normalize by the median flux in the interval 
;                  NORMWAVE +/- BINSIZE/2.0 (default mean wavelength) 
;       binsize  - see NORMWAVE (default 50.0)
;
; KEYWORD PARAMETERS:
;       mean   - normalize to the mean value in SPEC
;       median - normalize to the median value in SPEC
;       max    - normalize to the max value in SPEC
;
; OUTPUTS:
;       normspec - normalized SPEC
;
; OPTIONAL OUTPUTS:
;       const    - normalization constant
;
; PROCEDURES USED:
;       DJS_MEAN(), DJS_MEDIAN(), GET_ELEMENT
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Feb 20, U of A, written
;       jm03dec30uofa - if NORMWAVE is specified then compute the
;                       mean, median, or max within the specified 
;                       wavelength interval
;
; Copyright (C) 2003, John Moustakas
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

function im_normalize, spec, wave, normwave=normwave, $
  binsize=binsize, mean=mean, median=median, max=max, $
  const=const

    nwave = n_elements(wave)
    nspec = n_elements(spec)
    
    if (nspec eq 0L) then begin
       print, 'Syntax - normspec = im_normalize(spec,wave,normwave=,binsize=,$'
       print, '   const=,/mean,/median,/max)'
       return, -1L
    endif
    
    if (nwave ne 0L) then if (nwave ne nspec) then begin
       print, 'WAVE and SPEC do not have the same number of elements.'
       return, -1L       
    endif

    if (nwave ne 0L) then begin
       
       if n_elements(binsize) eq 0L then binsize = 50.0 ; [Angstrom] (could be smarter)
       if n_elements(normwave) eq 0L then normwave = djs_mean(wave) else normwave = float(normwave) 
       get_element, wave, normwave+[-binsize,+binsize]/2.0, x12

    endif
       
    if keyword_set(mean) then begin
       if keyword_set(normwave) then const = djs_mean(spec[x12[0]:x12[1]]) else const = djs_mean(spec)
       return, spec/const
    endif
    
    if keyword_set(max) then begin
       if keyword_set(normwave) then const = max(spec[x12[0]:x12[1]]) else const = max(spec)
       return, spec/const
    endif
    
    if keyword_set(median) then begin
       if keyword_set(normwave) then const = djs_median(spec[x12[0]:x12[1]]) else const = djs_median(spec)
       return, spec/const
    endif

    if keyword_set(normwave) then const = djs_median(spec[x12[0]:x12[1]]) else const = djs_median(spec)

return, spec/const
end
