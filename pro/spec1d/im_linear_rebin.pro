;+
; NAME:
;   IM_LINEAR_REBIN()
;
; PURPOSE:
;   Rebin a spectrum linearly in wavelength.
;
; INPUTS: 
;   wave - input wavelength vector with arbitrary pixel spacing
;     [Angstrom, NPIX]
;   flux - corresponding flux vector [NPIX]
;
; OPTIONAL INPUTS: 
;   ferr - corresponding error spectrum [NPIX]
;   minwave - minimum output wavelength [Angstrom] (default:
;     min(wave)) 
;   maxwave - maximum output wavelength [Angstrom] (default:
;     max(wave)) 
;   dwave - output pixel size [Angstrom] (defaults to the average
;     pixel size over the whole spectrum)
;
; KEYWORD PARAMETERS: 
;   conserve_intensity - by default, this routine conserves *flux*;
;     use this keyword to account for the change in the size of the 
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

function im_linear_rebin, wave, flux, minwave=minwave, maxwave=maxwave, $
  dwave=dwave, ferr=ferr, outwave=outwave, outferr=outferr, $
  conserve_intensity=conserve_intensity, borders=borders, $
  newborders=newborders

    compile_opt idl2
    on_error, 2

    npix = n_elements(wave)
    if (npix eq 0) or (npix ne n_elements(flux)) then begin
       doc_library, 'uly_log_regin'
       return, -1
    endif

    if (n_elements(minwave) eq 0) then minwave = min(wave)*1.0D
    if (n_elements(maxwave) eq 0) then maxwave = max(wave)*1.0D
    if (n_elements(dwave) eq 0) then dwave = $
      (((max(wave)*1.0D)-min(wave)*1.0D)/(npix-1.0D))[0]

; compute the borders of input bins
    borders = [(wave + shift(wave, 1)) / 2d, 0d]
    borders[0] = 2.*wave[0] - borders[1]
    borders[npix] = 2.*wave[npix-1] - borders[npix-1]

    linrange = [minwave,maxwave]
    nout = round((linRange[1]-linRange[0]) / dwave + 1)
    linStart = linRange[0]

    if linStart-dwave/2d gt borders[npix] then begin
       message, 'start value is not in the valid range', /INFO
       return, 0
    endif

; define the new borders of the bins
    NewBorders = linStart + (dindgen(nout+1)-0.5d) * dwave

; the rest of this routine is exactly the same as uly_spect_logrebin
; (shall we put it in a common low-level routine?)
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
       outflux = xrebin(borders,flux,NewBorders,/SPLINF)
    endelse

    outwave = linstart + dindgen(nout)*dwave

; handle the error spectrum    
    if (n_err eq n_data) then begin
       var = xrebin(borders,ferr^2,NewBorders,/SPLINF) 
       if keyword_set(conserve_intensity) then $
         outferr = sqrt(abs(var)) / flat else $
           outferr = sqrt(abs(var))
    endif

; deal with extrapolation
    mask = ((outwave ge min(wave)) and (outwave le max(wave)))
    outflux *= mask
    if (n_err eq n_data) then outferr *= mask
    
return, outflux
end

