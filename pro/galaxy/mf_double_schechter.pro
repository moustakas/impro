;+
; NAME:
;   MF_DOUBLE_SCHECHTER()
;
; PURPOSE:
;   Generate a double Schechter model given the relevant parameters.
;
; INPUTS: 
;   logmass - log-base-10 stellar masses (in the same units as MSTAR,
;     presumably Msun) at which to evaluate the model
;   phistar - number density at the 'knee' of the first Schechter
;     model 
;   logmstar - log-base-10 stellar mass at the 'knee' of the first
;     Schechter model 
;   alpha - low-mass slope of the first Schechter model
;   phistar2 - number density at the 'knee' of the second Schechter
;     model 
;   logmstar2 - log-base-10 stellar mass at the 'knee' of the second
;     Schechter model  
;   alpha2 - low-mass slope of the second Schechter model 
;
; OPTIONAL INPUTS: 
;   schechter - structure with PHISTAR, LOGMSTAR, ALPHA, PHISTAR2,
;     LOGMSTAR2, and ALPHA2 (takes precedence over PHISTAR, LOGMSTAR, 
;     ALPHA, PHISTAR2, LOGMSTAR2, and ALPHA)  
;
; OUTPUTS: 
;   model - logarithmic model Phi(log M) 
; 
; COMMENTS:
;   Note that Phi(log M) = ln(10)*M*Phi(M)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 201 Aug 23, UCSD
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

function mf_double_schechter, logmass, phistar, logmstar, alpha, phistar2, $
  logmstar2, alpha2, double_schechter=double_schechter

    if (n_elements(logmass) eq 0L) then begin
       doc_library, 'mf_schechter'
       return, -1
    endif

    if (n_elements(double_schechter) eq 0) then begin
       if (n_elements(phistar) eq 0) or (n_elements(logmstar) eq 0) or $
         (n_elements(alpha) eq 0) or (n_elements(phistar2) eq 0) or $
         (n_elements(logmstar2) eq 0) or (n_elements(alpha2) eq 0) then begin
          splog, 'PHISTAR, LOGMSTAR, ALPHA, PHISTAR2, LOGMSTAR2, and ALPHA2 inputs required'
          return, -1
       endif
       use_phistar = phistar
       use_logmstar = logmstar
       use_alpha = alpha
       use_phistar2 = phistar2
       use_logmstar2 = logmstar2
       use_alpha2 = alpha2
    endif else begin
       if (tag_exist(double_schechter,'phistar') eq 0) or $
         (tag_exist(double_schechter,'logmstar') eq 0) or $
         (tag_exist(double_schechter,'alpha') eq 0) or $
         (tag_exist(double_schechter,'phistar2') eq 0) or $
         (tag_exist(double_schechter,'logmstar2') eq 0) or $
         (tag_exist(double_schechter,'alpha2') eq 0) then begin
          splog, 'Improper SCHECHTER structure!'
          return, -1
       endif
       use_phistar = double_schechter.phistar
       use_logmstar = double_schechter.logmstar
       use_alpha = double_schechter.alpha
       use_phistar2 = double_schechter.phistar2
       use_logmstar2 = double_schechter.logmstar2
       use_alpha2 = double_schechter.alpha2
    endelse

    mratio = 10D^(logmass-use_logmstar)
    mratio2 = 10D^(logmass-use_logmstar2)
    model = alog(10)*(exp(-mratio)*use_phistar*mratio^(use_alpha+1) + $
      exp(-mratio2)*use_phistar2*mratio2^(use_alpha2+1))

return, model
end
