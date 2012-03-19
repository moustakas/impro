;+
; NAME:
;   MF_SCHECHTER_PLUS()
;
; PURPOSE:
;   Generate a Schechter plus additional low-mass component model
;   given the relevant parameters. 
;
; INPUTS: 
;   logmass - log-base-10 stellar masses (in the same units as MSTAR,
;     presumably Msun) at which to evaluate the model
;   phistar - number density at the 'knee' of the MF
;   logmstar - log-base-10 stellar mass at the 'knee' of the MF 
;   alpha - low-mass slope
;   phiplus - second normalization at low mass
;   alphaplus - second slope at low mass
;
; OPTIONAL INPUTS: 
;   schechter_plus - MF_SCHECHTER() style structure with PHISTAR,
;     LOGMSTAR, ALPHA, PHIPLUS, and ALPHAPLUS (takes precedence over
;     PHISTAR, LOGMSTAR, ALPHA, PHIPLUS, and ALPHAPLUS) 
;
; OUTPUTS: 
;   model - logarithmic double Schechter model Phi(log M) 
; 
; COMMENTS:
;   Note that Phi(log M) = ln(10)*M*Phi(M)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 11, UCSD - based largely on
;     M. Blanton's LF_FIT_SCHECHTER_PLUS
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

function mf_schechter_plus, logmass, phistar, logmstar, alpha, $
  phiplus, alphaplus, schechter_plus=schechter_plus

    if (n_elements(logmass) eq 0L) then begin
       doc_library, 'mf_schechter_plus'
       return, -1
    endif

    if (n_elements(schechter_plus) eq 0) then begin
       if (n_elements(phistar) eq 0) or (n_elements(logmstar) eq 0) or $
         (n_elements(alpha) eq 0) or (n_elements(phiplus) eq 0) or $
         (n_elements(alphaplus) eq 0) then begin
          splog, 'PHISTAR, LOGMSTAR, ALPHA, PHIPLUS, and ALPHAPLUS inputs required'
          return, -1
       endif
       use_phistar = phistar
       use_logmstar = logmstar
       use_alpha = alpha
       use_phiplus = phiplus
       use_alphaplus = alphaplus
    endif else begin
       if (tag_exist(schechter_plus,'phistar') eq 0) or $
         (tag_exist(schechter_plus,'logmstar') eq 0) or $
         (tag_exist(schechter_plus,'alpha') eq 0) or $
         (tag_exist(schechter_plus,'phiplus') eq 0) or $
         (tag_exist(schechter_plus,'alphaplus') eq 0) then begin
          splog, 'Improper SCHECHTER_PLUS structure!'
          return, -1
       endif
       use_phistar = schechter_plus.phistar
       use_logmstar = schechter_plus.logmstar
       use_alpha = schechter_plus.alpha
       use_phiplus = schechter_plus.phiplus
       use_alphaplus = schechter_plus.alphaplus
    endelse

    mratio = 10D^(logmass-use_logmstar)
    model = alog(10)*exp(-mratio)*(use_phistar*mratio^(use_alpha+1) + $
      use_phiplus*mratio^(use_alphaplus+1))

return, model
end
