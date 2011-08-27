;+
; NAME:
;       IM_GAUSS_BROADEN()
;
; PURPOSE:
;       Degrade the spectral resolution of a spectrum. 
;
; CALLING SEQUENCE:
;       broadflux = im_gauss_broaden(wave,flux,oldres,newres)
;
; INPUTS:
;       wave   - wavelength vector (Angstrom) [NPIX]
;       flux   - corresponding flux vector (either an [NPIX] vector or
;                an [NPIX,NSTAR] array)
;       oldres - original spectral resolution that can be either a
;                scalar or an [NPIX] vector [Angstrom, FWHM]
;       newres - new, output spectral resolution [Angstrom, FWHM]
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       broadflux - FLUX array degraded in spectral resolution
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       GCONVOLVE()
;
; COMMENTS:
;       Note that if either the spectral dispersion, OLDRES, or NEWRES
;       are not constant then this routine is quite slow (because we
;       have to loop on each pixel in the flux array).
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 November 2, U of A
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

function im_gauss_broaden, wave, flux, oldres, newres

    if (n_params() ne 4L) then begin
       print, 'Syntax = broadflux = im_gauss_broaden(wave,flux,oldres,newres)'
       return, -1L
    endif
    
    npix = n_elements(wave)

    ndim = size(flux,/n_dimension)
    dim = size(flux,/dimension)

    if (ndim ne 1L) and (ndim ne 2L) then begin
       print, 'FLUX must be a one- or two-dimensional array.'
       return, flux
    endif
    
    if (ndim eq 2L) then begin

       nfluxpix = dim[0]
       nstar = dim[1]

    endif else begin

       nfluxpix = n_elements(flux)
       nstar = 1L
       flux = reform(flux,nfluxpix,1)

    endelse
    
    if (npix ne nfluxpix) then begin
       print, 'WAVE and FLUX do not have the same number of elements.'
       return, flux
    endif

    newres = float(newres)
    oldres = float(oldres)
    
    if max(newres) lt max(oldres) then begin
       print, 'NEWRES must be larger than OLDRES at all wavelengths.'
       return, flux
    endif

; if NEWRES has only one uniq value then make this routine faster by
; choosing a fixed SIGMA

    if n_elements(uniq(newres)) eq 1L then newres = newres[uniq(newres)]
    
    noldrespix = n_elements(oldres)
    nnewrespix = n_elements(newres)

    if (noldrespix ne 1L) and (noldrespix ne npix) then begin
       print, 'OLDRES must be either a scalar or an NPIX length vector.'
       return, flux
    endif

    if (nnewrespix ne 1L) and (nnewrespix ne npix) then begin
       print, 'NEWRES must be either a scalar or an NPIX length vector.'
       return, flux
    endif

; do not assume the spectral dispersion is constant; pad the
; wavelength compute the sigma width of the Gaussian convolution
; kernel in pixels

    waveone = fltarr(npix+1)
    waveone[1L:npix] = wave
    waveone[0L] = interpol(wave,lindgen(npix),-1)
    waveshift = shift(waveone,1)
    
    disp = wave - waveshift[1:npix]
    nuniqdisp = n_elements(uniq(disp))

    fwhm2sig = 2.0*sqrt(2.0*alog(2.0))                      ; = 2.35
    sigma = sqrt(newres^2.0 - oldres^2.0) / fwhm2sig / disp ; [pixels]

    broadflux = flux*0.0
;   broadflux = reform(flux*0.0,npix,nstar)

    if (noldrespix gt 1L) or (nnewrespix gt 1L) or (nuniqdisp ne 1L) then begin

       splog, 'WARNING: Entering buggy code!'
       for j = 0L, nstar-1L do for k = 0L, npix-1L do $
         broadflux[k,j] = (gconvolve(reform(flux[*,j]),sigma[k],/edge_truncate))[k] 

    endif else begin

      for j = 0L, nstar-1L do broadflux[*,j] = gconvolve(reform(flux[*,j]),sigma[0],/edge_truncate)

   endelse

return, reform(broadflux)
end
