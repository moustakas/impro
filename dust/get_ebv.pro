;+
; NAME:
;       GET_EBV()
;
; PURPOSE:
;       Get E(B-V) from the Balmer decrement and an extinction curve. 
;
; INPUTS:
;       decrement - the observed Balmer decrement, for example,
;                   F(H-alpha)/F(H-beta)  
;
; OPTIONAL INPUTS:
;       decrement_err - error in DECREMENT
;       temperature   - temperature used to derive the intrinsic
;                       Balmer decrement [K] (default 10000)
;       use_hahb      - instead of the default 2.86, use this
;                       intrinsic Ha/Hb ratio; useful for, e.g., AGN,
;                       for which Ha/Hb is 3.1
;
; KEYWORD PARAMETERS:
;       HaHb  - de-redden using the H-alpha/H-beta line ratio
;               (default) 
;       HbHg  - de-redden using the H-beta/H-gamma line ratio
;       HaHg  - de-redden using the H-alpha/H-gamma line ratio
;       extra - keywords for K_LAMBDA()
;
; OUTPUTS:
;       ebv - color excess E(B-V) [mag]
;
; OPTIONAL OUTPUTS:
;       ebv_err   - error in EBV [mag]
;       color     - color excess = -2.5 log (R_int/R_obs), computed
;                   directly from the observed and intrinsic decrement
;                   [mag] 
;       err_color - error in COLOR [mag]
;
; COMMENTS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 May 28, U of A - written
;       jm02nov6uofa  - added K_LAMBDA() compatibility
;       jm03jun13uofa - added HGAMMA keyword
;       jm03aug06uofa - added temperature dependence on the Balmer
;                       decrement and HBETA keyword
;       jm04nov03uofa - removed HGAMMA and HBETA keywords in favor of
;                       HaHb, HaHg, and HbHg keywords
;       jm05feb03uofa - if no temperature is given, then compute the
;                       error in R_int, otherwise assume the error is
;                       zero; return COLOR and ERR_COLOR keywords 
;       jm08feb06nyu  - added USE_HAHB optional input
;
; Copyright (C) 2002-2005, 2008, John Moustakas
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

function get_ebv, decrement, decrement_err=decrement_err, temperature=temperature1, $
  ebv_err=ebv_err, use_hahb=use_hahb, color=color, err_color=err_color, HaHb=HaHb, $
  HaHg=HaHg, HbHg=HbHg, _extra=extra

    ndec = n_elements(decrement)
    if ndec eq 0L then begin
       doc_library, 'get_ebv'
       return, -1
    endif

    ndecerr = n_elements(decrement_err)
    if ndecerr eq 0L then decrement_err = decrement*0.0 else begin
       if ndecerr ne ndec then begin
          print, 'DECREMENT and DECREMENT_ERR have incompatible dimensions.'
          return, -1L
       endif
    endelse

    if (n_elements(temperature1) eq 0L) then temperature = 10000.0 else temperature = temperature1

    mintemp = 5E3 ; minimum temperature
    maxtemp = 2E4 ; maximum temperature
    
; intrinsic H-alpha/H-beta ratio for type B recombination (default) 

    if (not keyword_set(HbHg)) and (not keyword_set(HaHg)) then HaHb = 1L
    
    if keyword_set(HaHb) then begin

       wave1 = 4861.33
       wave2 = 6562.80
       if (n_elements(use_hahb) eq 0L) then R_int = 2.86 else R_int = use_hahb
       R_int_err = R_int*0.0
;      R_int_err = R_int*0.05

;      R_int = return_tbalmer(temperature,/HaHb)
;      R_int_err = mean([R_int-return_tbalmer(maxtemp,/Hahb), return_tbalmer(mintemp,/Hahb)-R_int])

    endif
       
; intrinsic H-beta/H-gamma ratio for type B recombination

    if keyword_set(HbHg) then begin

       wave1 = 4340.464
       wave2 = 4861.33

       R_int = 2.14
       R_int_err = R_int*0.03

;      R_int = return_tbalmer(temperature,/HbHg)
;      R_int_err = mean([R_int-return_tbalmer(maxtemp,/HbHg), return_tbalmer(mintemp,/HbHg)-R_int])

    endif

; intrinsic H-alpha/H-gamma ratio for type B recombination

    if keyword_set(HaHg) then begin

       wave1 = 4340.464
       wave2 = 6562.80

       R_int = 6.11
       R_int_err = R_int*0.07

;      R_int = return_tbalmer(temperature,/HaHg)
;      R_int_err = mean([R_int-return_tbalmer(maxtemp,/HaHg), return_tbalmer(mintemp,/HaHg)-R_int])

    endif
    
    if (n_elements(temperature1) ne 0L) then R_int_err = 0.0 ; assume the temperature is known perfectly

    rcurve = k_lambda(wave1,_extra=extra)-k_lambda(wave2,_extra=extra)
    rcurve_err = 0.0

    ratio = R_int/decrement
    ratio_err = sqrt( (R_int_err/decrement)^2.0 + ( (R_int/decrement^2.0)*decrement_err)^2 )
    
    color = -2.5 * alog10(ratio)
    err_color = 2.5 * ratio_err / ratio  / alog(10.0)

    ebv = color / rcurve
    ebv_err = sqrt( (err_color/rcurve)^2.0 + ( (color/rcurve^2.0)*rcurve_err)^2 )

return, ebv
end
