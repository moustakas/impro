;+
; NAME:
;   IM_LOG_REBIN()
;
; PURPOSE:
;   Rebin a spectrum logarithmically (natural log) in wavelength. 
;
; INPUTS: 
;   wave - input wavelength vector with arbitrary pixel spacing
;     [Angstrom, NPIX]
;   flux - corresponding flux vector [NPIX]
;
; OPTIONAL INPUTS: 
;   vsc - velocity scale, i.e. size of each pixel in the log-rebinned
;     spectrum; if not specified, then VSC is computed from the input
;     spectrum  [km/s]
;   ferr - corresponding error spectrum [NPIX]
;   minwave - minimum output wavelength [Angstrom] (default:
;     min(wave)) 
;   maxwave - maximum output wavelength [Angstrom] (default:
;     max(wave)) 
;
; KEYWORD PARAMETERS: 
;   conserve_intensity - by default, this routine conserved *flux*;
;     use this routine to account for the change in the size of the 
;     pixels in order to conserve intensity
;
; OUTPUTS: 
;   outflux - output flux vector [NOUTPIX] 
;
; OPTIONAL OUTPUTS:
;   outwave - output wavelength vector [Angstrom, NOUTPIX] 
;   outferr - output err spectrum [NOUTPIX] 
;
; COMMENTS:
;   Double precision is strongly recommended!
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009-Nov-30, UCSD - code hacked and adapted from
;     original code written by Mina Koleva 
;
; Copyright (C) 2009, John Moustakas
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

function im_log_rebin, wave, flux, vsc=vsc, ferr=ferr, outwave=outwave, $
  outferr=outferr, minwave=minwave, maxwave=maxwave, $
  conserve_intensity=conserve_intensity

    compile_opt idl2
;   on_error, 2

    npix = n_elements(wave)
    if (npix eq 0) or (npix ne n_elements(flux)) then begin
       doc_library, 'uly_log_regin'
       return, -1
    endif

    if (n_elements(minwave) eq 0) then minwave = min(wave)*1.0D
    if (n_elements(maxwave) eq 0) then maxwave = max(wave)*1.0D
    waverange = [minwave,maxwave]
    
    c = 299792.458d             ; Speed of light in km/s

    borders = [(wave + shift(wave, 1)) / 2d, 0d]
    borders[0] = 2.*wave[0] - borders[1]
    borders[npix[0]] = 2.*wave[npix[0]-1] - borders[npix[0]-1]
    bordersLog = alog(borders)

; Determine the velocity scale, vsc, preserving the total number of pixels
    if n_elements(vsc) eq 0 then begin
;       We assume that the wavelength are ordered, and we do not check it
       if n_elements(waverange) gt 0 then begin
          nw = value_locate(wave, waverange)
          if nw[0] eq -1 then nw[0] = 0
          if n_elements(waverange) eq 1 then nw = [nw[0], npix[0]]
          if nw[1] eq -1 then nw[1] = 0 ; this is bad
       endif else nw = [0d, npix[0]-1]
       wr = wave[nw]
       vsc =  alog(wr[1]/wr[0]) / npix[0] * c 
    endif 

    logScale = vsc/c            ; Scaling in the output spectrum

; Determine the start wavelength, logStart, and number of pix, nout, in output
    logRange = alog([wave[0],wave[npix-1]])
    logRange += [0.5d,-0.5d]*logScale
    if n_elements(waverange) gt 0 then begin
;  the 1D-7 is a tolerance for rounding errors
       nshift = ceil(max([0d, (logRange[0] - alog(waverange[0])) / logScale - 1d-7]))
       logRange[0] = alog(waverange[0]) + logScale * nshift ; center of 1st out pix

       if n_elements(waverange) eq 2 then $
         logRange[1] = min([alog(waverange[1]), logRange[1]]) 

       if logRange[1] lt logRange[0] then begin
          message, /INFO, 'waverange is not valid: '+ $
            strjoin(strtrim(waverange,2),',')
          return, 0
       endif
    endif 

    nout = round((logRange[1]-logRange[0]) / logScale + 1)
    logStart = logRange[0]

    if logStart-logScale/2d gt bordersLog[npix[0]] then begin
       message, 'start value is not in the valid range', /INFO
       return, -1
    endif

; define the new borders of the bins
    NewBorders = exp(logStart + (dindgen(nout+1)-0.5d) * logScale)

    dim = size(flux, /DIM)
    n_data = n_elements(flux)
    n_err = n_elements(ferr)
    n_dim = size(flux, /N_DIM)

    if keyword_set(conserve_intensity) then begin
; Determine the conversion factor/vector accounting for pixel size change
       flat = fltarr(dim[0], /NOZERO) ; faster to create an array and
       replicate_inplace, flat, 1B    ; populate with 1s 
       flat = xrebin(borders, flat, NewBorders, /SPLINF) 
       if n_dim gt 1 then flat = rebin(flat, [nout, dim[1:*]])
       outflux = xrebin(borders,flux,NewBorders,/SPLINF) / flat 
    endif else begin
       outflux = xrebin(borders,flux,NewBorders,/cubic)
;      outflux = xrebin(borders,flux,NewBorders,/SPLINF)
    endelse

    outwave = logstart + dindgen(nout)*logscale

; handle the error spectrum    
    if (n_err eq n_data) then begin
       var = xrebin(borders,(1.0d*ferr)^2,NewBorders,/SPLINF)
       if keyword_set(conserve_intensity) then $
         outferr = sqrt(abs(var)) / flat else $
           outferr = sqrt(abs(var))
    endif

return, outflux
end
