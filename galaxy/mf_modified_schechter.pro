;+
; NAME:
;   MF_MODIFIED_SCHECHTER()
;
; PURPOSE:
;   Generate a modified Schechter model given the relevant parameters.
;
; INPUTS: 
;   logmass - log-base-10 stellar masses (in the same units as MSTAR,
;     presumably Msun) at which to evaluate the model
;   phistar - number density at the 'knee' of the MF
;   logmstar - log-base-10 stellar mass at the 'knee' of the MF 
;   alpha - low-mass slope
;   beta - second power-law parameter
;
; OPTIONAL INPUTS: 
;   modified_schechter - MF_SCHECHTER() style structure with PHISTAR,
;     LOGMSTAR, ALPHA, and BETA (takes precedence over PHISTAR,
;     LOGMSTAR, ALPHA, and BETA)
;
; OUTPUTS: 
;   model - logarithmic double Schechter model Phi(log M) 
; 
; COMMENTS:
;   Note that Phi(log M) = ln(10)*M*Phi(M).  This model is used by
;   Bernardi+10, MNRAS.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Aug 26, UCSD
;
; Copyright (C) 2011, John Moustakas
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

function mf_modified_schechter, logmass, phistar, logmstar, alpha, $
  beta, modified_schechter=modified_schechter

    if (n_elements(logmass) eq 0L) then begin
       doc_library, 'mf_modified_schechter'
       return, -1
    endif

    if (n_elements(modified_schechter) eq 0) then begin
       if (n_elements(phistar) eq 0) or (n_elements(logmstar) eq 0) or $
         (n_elements(alpha) eq 0) or (n_elements(beta) eq 0) then begin
          splog, 'PHISTAR, LOGMSTAR, ALPHA, and BETA inputs required'
          return, -1
       endif
       use_phistar = phistar
       use_logmstar = logmstar
       use_alpha = alpha
       use_beta = beta
    endif else begin
       if (tag_exist(modified_schechter,'phistar') eq 0) or $
         (tag_exist(modified_schechter,'logmstar') eq 0) or $
         (tag_exist(modified_schechter,'alpha') eq 0) or $
         (tag_exist(modified_schechter,'beta') eq 0) then begin
          splog, 'Improper MODIFIED_SCHECHTER structure!'
          return, -1
       endif
       use_phistar = modified_schechter.phistar
       use_logmstar = modified_schechter.logmstar
       use_alpha = modified_schechter.alpha
       use_beta = modified_schechter.beta
    endelse

    mratio = 10D^(logmass-use_logmstar)
    model = alog(10)*use_phistar*mratio^use_alpha*$
      exp(-mratio^use_beta)*use_beta/gamma(use_alpha/use_beta)

return, model
end
