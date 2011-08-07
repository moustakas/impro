;+
; NAME:
;       IM_TOTAL_FIR()
;
; PURPOSE:
;       Compute the FAR-infrared luminosity given the IRAS fluxes at
;       12, 25, 60, and 100 microns. 
;
; INPUTS:
;       iras_flux - IRAS flux densities at 12, 25, 60, and 100 microns
;                   in Jansky [4,NGALAXY]
;
; OPTIONAL INPUTS:
;       iras_ferr - errors corresponding to IRAS_FLUX in Jansky
;                   [4,NGALAXY] 
;       alpha     - extend the data beyond 100 microns assuming a
;                   lambda^(-alpha) modified blackbody curve 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       result - output data structure (see documentation below)
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       DJS_PLANCK(), QPINT1D()
;
; COMMENTS:
;       L(FIR) is defined between 40 and 500 microns.  Upper limits are not treated.
;
;       DJS_PLANCK() is available in David Schlegel's IDLUTILS library
;       and QPINT1D() is Craig Markwardt's very robust numerical
;       integration routine. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       S. Juneau, 2008, U of A, modified from im_total_ir.pro to get
;       FIR flux instead of IR flux.
;
;       J. Moustakas, 2004 Feb 25, U of A, based in large part on code  
;          originally written by Karl Gordon
;
; Copyright (C) 2004, John Moustakas
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

function qpint1d_func, x, wave=wave, flux=flux
return, interpol(flux,wave,x)
end

function im_total_fir, iras_flux, iras_ferr=iras_ferr, alpha=alpha

    if (n_elements(iras_flux) eq 0L) then begin
       doc_library, 'im_total_ir'
       return, -1L
    endif

    ndim = size(iras_flux,/n_dimension)
    dims = size(iras_flux,/dimension)

    if (n_elements(iras_ferr) eq 0L) then iras_ferr = iras_flux*0.0

    ndim_err = size(iras_ferr,/n_dimension)
    dims_err = size(iras_ferr,/dimension)

    if (ndim_err ne ndim) then begin
       splog, 'IRAS_FLUX and IRAS_FERR have incompatible dimensions.'
       return, -1L
    endif
    
; if there are multiple objects then call this routine recursively 
    
    if (ndim eq 2L) then begin

       ngalaxy = dims[1]
       for igalaxy = 0L, ngalaxy-1L do begin
          result1 = im_total_fir(iras_flux[*,igalaxy],iras_ferr=iras_ferr[*,igalaxy],alpha=alpha)
          if (igalaxy eq 0L) then result = result1 else result = [ [result], [result1] ]
       endfor
       return, reform(result)
       
    endif
    
    if (n_elements(alpha) eq 0L) then alpha = 1.0

; initialize the output data structure

    result = {$
      iras_flux: iras_flux, $ ; input fluxes [Jy]
      iras_ferr: iras_ferr, $ ; input flux errors [Jy]
      dust_temp:      -1.0, $ ; dust temperature corresponding to f(60)/f(100) ratio
      dust_temp_min:  -1.0, $ ; minimum possible dust temperature
      dust_temp_max:  -1.0, $ ; maximum possible dust temperature
      fir_flux:       -1.0, $ ; far-IR flux from Helou et al. (1988) [erg/s/cm2]
      fir_ferr:       -1.0, $ ; error in above
      fir40_flux:     -1.0, $ ; far-IR (40-500 micron) flux [erg/s/cm2]
      fir40_ferr:     -1.0, $ ; error in above
      ir_flux_short:  -1.0, $ ; 40-100 micron flux [erg/s/cm2]
      ir_ferr_short:  -1.0, $ ; error in above
      ir_flux_long:   -1.0, $ ; 100-500 micron flux [erg/s/cm2]
      ir_ferr_long:   -1.0}   ; error in above
    
; compute the FIR flux according to the formula in Helou et al. (1988) 

    result.fir_flux = 1D3*1.26D-14*(2.58*iras_flux[2]+iras_flux[3])           ; [erg/s/cm2]
    result.fir_ferr = 1D3*1.26D-14*sqrt((2.58*iras_ferr[2])^2+iras_ferr[3]^2) ; [erg/s/cm2]

; initialize parameters of the IRAS satellite and some constants 

    ir_wave = [12.0,25.0,60.0,100.0]      ; central wavelength [micron]
    ir_delta = [13.0,24.0,44.0,40.0]      ; bandpass widths [micron]

    light = 2.99792458D14            ; speed of light [micron/s]
    cf = 1D-26*1D3*light/(ir_wave^2) ; [Jy] --> [erg/s/cm2/micron]
    ir_flux = cf*iras_flux           ; [erg/s/cm2/micron]
    ir_ferr = cf*iras_ferr           ; [erg/s/cm2/micron]

    ir_wave_extended = ir_wave
    ir_flux_extended = ir_flux
    ir_ferr_extended = ir_ferr
    
; determine the temperature of a blackbody corresponding to the
; observed f(60)/f(100) ratio; since this is an implicit equation, get
; the temperature by numerical interpolation; also figure out the
; minimum and maximum dust temperatures allowed by the data

    ratio = ir_flux[2]/ir_flux[3]
    ratio_min = (ir_flux[2]-ir_ferr[2])/(ir_flux[3]+ir_ferr[3]) ; minimum ratio (temperature)
    ratio_max = (ir_flux[2]+ir_ferr[2])/(ir_flux[3]-ir_ferr[3]) ; maximum ratio (temperature)

    temperature = findgen((100.0-10.0)/1.0)*1.0+10.0 ; [K]

    BB_60  = djs_planck(temperature,ir_wave[2],/MJy)*cf[2]*ir_wave[2]^(-alpha)
    BB_100 = djs_planck(temperature,ir_wave[3],/MJy)*cf[3]*ir_wave[3]^(-alpha)

    result.dust_temp = interpol(temperature,BB_60/BB_100,ratio)
    result.dust_temp_min = interpol(temperature,BB_60/BB_100,ratio_min)
    result.dust_temp_max = interpol(temperature,BB_60/BB_100,ratio_max)
    
; integrate the flux density between 40 and 100 microns

    result.ir_flux_short = qpint1d('qpint1d_func',40.,$
      ir_wave_extended[3],functargs={wave:ir_wave_extended,flux:ir_flux_extended})
    result.ir_ferr_short = sqrt(qpint1d('qpint1d_func',40.,$
      ir_wave_extended[3],functargs={wave:ir_wave_extended,flux:ir_ferr_extended^2}))

; to get the total flux longward of 100 microns, construct a blackbody
; at the appropriate dust temperature, normalized at the observed
; 100-micron flux density; then simply integrate

    minwave = 40.0    ; [micron]
    maxwave = 500.0 ; [micron]
    dlogwave = 0.001
    
    nwave = round(1.0+(alog10(maxwave)-alog10(minwave))/dlogwave)
    logwave = alog10(minwave) + dlogwave*findgen(nwave)
    wave = float(10^logwave) ; [micron]

    cf2 = 1D-26*1D3*1D6*light/wave^2 ; [MJy/sr] --> [erg/s/cm2/micron]

    BB = djs_planck(result.dust_temp,wave,/MJy)*cf2*wave^(-alpha) ; [erg/s/cm2]
    BB = ir_flux[3]*BB/interpol(BB,wave,ir_wave[3])               ; [erg/s/cm2]

; generate blackbodies corresponding to the minimum and maximum
; allowable temperatures; remember! a higher temperature means less
; flux longward of 100 microns!
    
    BB_min = djs_planck(result.dust_temp_min,wave,/MJy)*cf2*wave^(-alpha)    ; [erg/s/cm2]
    BB_min = (ir_flux[3]+ir_ferr[3])*BB_min/interpol(BB_min,wave,ir_wave[3]) ; [erg/s/cm2]
    
    BB_max = djs_planck(result.dust_temp_max,wave,/MJy)*cf2*wave^(-alpha)    ; [erg/s/cm2]
    BB_max = (ir_flux[3]-ir_ferr[3])*BB_max/interpol(BB_max,wave,ir_wave[3]) ; [erg/s/cm2]

; integrate

    ir_flux_long_min = qpint1d('qpint1d_func',100.0,500.0,functargs={wave:wave,flux:BB_min})
    ir_flux_long_max = qpint1d('qpint1d_func',100.0,500.0,functargs={wave:wave,flux:BB_max})

    result.ir_flux_long = (ir_flux_long_min + ir_flux_long_max) / 2.0
    result.ir_ferr_long = (ir_flux_long_min - ir_flux_long_max) / 2.0

    result.fir40_flux = result.ir_flux_short + result.ir_flux_long
    result.fir40_ferr = sqrt(result.ir_ferr_short^2 + result.ir_ferr_long^2)
    
return, result
end
