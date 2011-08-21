;+
; NAME:
;       IM_BINSPEC()
;
; PURPOSE:
;       Bin a spectrum.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;       spec - one-dimensional spectrum or array
;       wave - corresponding wavelength vector or index number 
;
; OPTIONAL INPUTS:
;       spec_err  - error spectrum for SPEC.  if it exists then
;                   compute BINSPEC_ERR
;       binsize   - sum the flux in bins of width BINSIZE
;       bandstart - central wavelength or pixel number of the first bin
;       bandend   - central wavelength or pixel number of the last bin
;       finerpix  - over-sampling factor at the bandpass edges
;                   (default 10)
;
; KEYWORD PARAMETERS:
;       abmag  - assume that SPEC is in [erg/s/cm2/Hz] and convert to
;                AB magnitude 
;       doplot - generate a plot of the input and the binned spectrum 
;
; OUTPUTS:
;       binspec - binned spectrum
;
; OPTIONAL OUTPUTS:
;       binwave     - central bandpass wavelengths or pixel numbers 
;       binspec_err - the error spectrum for BINSPEC if SPEC_ERR was
;                     passed (not supported yet!)
;
; COMMENTS:
;       Although this routine was written with spectra in mind, it
;       should be general enough to bin any one-dimensional array. 
;
; PROCEDURES USED:
;       LEGEND
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 13, U of A, based on earlier codes 
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

function im_binspec, spec, wave, spec_err=spec_err, binsize=binsize, $
  bandstart=bandstart, bandend=bandend, finerpix=finerpix, binwave=binwave, $
  binspec_err=binspec_err, abmag=abmag, doplot=doplot

    nspec = n_elements(spec)
    nwave = n_elements(wave)

    if (nspec eq 0L) or (nwave eq 0L) then begin
       print, 'Syntax - '
       return, -1L
    endif

    if (nspec ne nwave) then begin
       print, 'SPEC and WAVE must have the same number of elements.'
       return, -1L       
    endif

    nspec_err = n_elements(spec_err)
    if nspec_err eq 0L then spec_var = spec*0.0+1.0 else begin
       if nspec ne nspec_err then begin
          print, 'SPEC and SPEC_ERR must have the same number of elements.'
          return, -1L       
       endif
       spec_var = 1.0/spec_err^2.0          
    endelse

    if n_elements(finerpix) eq 0L then finerpix = 10.0
    if n_elements(binsize) eq 0L then binsize = 50.0 else binsize = float(binsize)
    if n_elements(bandstart) eq 0L then bandstart = min(wave) + (fix(2*binsize) - (min(wave) mod fix(2*binsize)))
    if n_elements(bandend) eq 0L then bandend = max(wave) - (max(wave) mod fix(2*binsize))

;   if n_elements(bandstart) eq 0L then bandstart = min(wave) + (100 - (min(wave) mod 100))
;   if n_elements(bandend) eq 0L then bandend = max(wave) - (max(wave) mod 100)

    nbins = long((bandend-bandstart)/binsize)
    binwave = findgen(nbins)*binsize+bandstart ;+binsize/2.0

    binspec = fltarr(nbins)
    binspec_err = binspec*0.0

    for j = 0L, nbins-1L do begin
       
       wavearray = bandstart+j*binsize-binsize/2.0+findgen(binsize*finerpix)/finerpix+1.0/finerpix/2.0

       bflux = interpol(spec,wave,wavearray)
       bvar = interpol(spec_var,wave,wavearray)/finerpix

;      binspec[j] = total(bvar*bflux)/binsize/finerpix/total(bvar)
       binspec[j] = total(bflux)/binsize/finerpix
       if nspec_err ne 0L then binspec_err = 1.0/sqrt(total(bvar))

    endfor

    if keyword_set(abmag) then binspec = -2.5*alog10(binwave*binwave*binspec/2.99793D18)-48.59 ; AB magnitude

    if keyword_set(doplot) then begin

       normspec = im_normalize(spec,wave,/max,const=norm)
       plot, wave, normspec, ps=10, xsty=3, ysty=3, $
         ytitle=''+string(1.0/norm,format='(G0)')+' x SPEC', xtitle='WAVE', $
         charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, thick=2.0
       oplot, binwave, binspec/norm, ps=4, syms=2.0, thick=2.0
       legend, ['BINSIZE = '+string(binsize,format='(G0)'),$
         'BANDSTART = '+string(bandstart,format='(G0)'),$
         'BANDEND   = '+string(bandend,format='(G0)')], /right, /top, $
         box=0, charsize=2.0, charthick=2.0
       
    endif
    
return, binspec
end
