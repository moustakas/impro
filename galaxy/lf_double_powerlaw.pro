;+
; NAME:
;   LF_DOUBLE_POWERLAW()
;
; PURPOSE:
;   Generate a double power-law luminosity function given the relevant
;   parameters. 
;
; INPUTS: 
;   loglum - log-base-10 luminosity at which to evaluate the model
;     [NLUM] 
;   phistar - number density at the inflection point of the LF 
;   loglstar - luminosity at the inflection point of the LF
;   alpha - first power law slope
;   beta - second power law slope
;
; OPTIONAL INPUTS: 
;   powerlaw - LF_FIT_DOUBLE_POWERLAW() style structure with PHISTAR, 
;     LOGLSTAR, ALPHA, and BETA (takes precedence over PHISTAR,
;     LOGLSTAR, ALPHA, and BETA)
;
; OUTPUTS: 
;   model - output model [NLUM]
; 
; COMMENTS:
;   This type of model is frequently used to fit the X-ray or mid-IR
;   luminosity function; see, for example, Yahil+91; Aird+08; or
;   Rujopakarn+10.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Apr 14, NYU
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

function lf_double_powerlaw, loglum, phistar, loglstar, $
  alpha, beta, powerlaw=powerlaw

    if (n_elements(loglum) eq 0L) then begin
       doc_library, 'lf_double_powerlaw'
       return, -1
    endif

    if (n_elements(powerlaw) eq 0) then begin
       if (n_elements(phistar) eq 0) or (n_elements(loglstar) eq 0) or $
         (n_elements(alpha) eq 0) or (n_elements(beta) eq 0) then begin
          splog, 'PHISTAR, LOGLSTAR, ALPHA, and BETA inputs required'
          return, -1
       endif
       use_phistar = phistar
       use_loglstar = loglstar
       use_alpha = alpha
       use_beta = beta
    endif else begin
       if (tag_exist(powerlaw,'phistar') eq 0) or $
         (tag_exist(powerlaw,'loglstar') eq 0) or $
         (tag_exist(powerlaw,'alpha') eq 0) or $
         (tag_exist(powerlaw,'beta') eq 0) then begin
          splog, 'Improper POWERLAW structure!'
          return, -1
       endif
       use_phistar = powerlaw.phistar
       use_loglstar = powerlaw.loglstar
       use_alpha = powerlaw.alpha
       use_beta = powerlaw.beta
    endelse

; integrated and differential LF (see Yahil+91)
    intlf = use_phistar*(10D^(loglum-use_loglstar))^(-use_alpha)*(1.0+10D^(loglum-use_loglstar))^(-use_beta)
    model = alog(10)*10D^loglum*(use_alpha/10D^loglum+use_beta/(10D^loglum+10D^use_loglstar))*intlf
    
return, model
end
