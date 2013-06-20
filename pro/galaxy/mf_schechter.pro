;+
; NAME:
;   MF_SCHECHTER()
;
; PURPOSE:
;   Generate a Schechter model given the relevant parameters.
;
; INPUTS: 
;   logmass - log-base-10 stellar masses (in the same units as MSTAR,
;     presumably Msun) at which to evaluate the model
;   phistar - number density at the 'knee' of the MF
;   logmstar - log-base-10 stellar mass at the 'knee' of the MF 
;   alpha - low-mass slope
;
; OPTIONAL INPUTS: 
;   schechter - FIT_MF_SCHECHTER() style structure with PHISTAR,
;     LOGMSTAR, and ALPHA (takes precedence over PHISTAR, LOGMSTAR,
;     and ALPHA)  
;
; OUTPUTS: 
;   model - logarithmic Schechter model Phi(log M) 
; 
; COMMENTS:
;   Note that Phi(log M) = ln(10)*M*Phi(M)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Apr 14, NYU - based largely on
;     M. Blanton's LF_FIT_SCHECHTER 
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

function mf_schechter, logmass, phistar, logmstar, alpha, $
  schechter=schechter, numden=numden, rhoden=rhoden
    
    if (n_elements(logmass) eq 0L) then begin
       doc_library, 'mf_schechter'
       return, -1
    endif

    if (n_elements(schechter) eq 0) then begin
       if (n_elements(phistar) eq 0) or (n_elements(logmstar) eq 0) or $
         (n_elements(alpha) eq 0) then begin
          splog, 'PHISTAR, LOGMSTAR, and ALPHA inputs required'
          return, -1
       endif
       use_phistar = phistar
       use_logmstar = logmstar
       use_alpha = alpha
    endif else begin
       if (tag_exist(schechter,'phistar') eq 0) or $
         (tag_exist(schechter,'logmstar') eq 0) or $
         (tag_exist(schechter,'alpha') eq 0) then begin
          splog, 'Improper SCHECHTER structure!'
          return, -1
       endif
       use_phistar = schechter.phistar
       use_logmstar = schechter.logmstar
       use_alpha = schechter.alpha
    endelse

    mratio = 10D^(logmass-use_logmstar)
    model = alog(10)*exp(-mratio)*use_phistar*mratio^(use_alpha+1)

; integrate the number and mass density
;   if arg_present(numden) then numden = use_phistar*igamma(use_alpha+1,10D^((7D)-use_logmstar))
;   if arg_present(numden) then numden = use_phistar*gamma(use_alpha+1)*$
;     (1-igamma(use_alpha+1,10D^((5D)-use_logmstar)))
;   if arg_present(rhoden) then rhoden = use_phistar*gamma(use_alpha+2)*$
;     (1-igamma(use_alpha+2,10D^((5D)-use_logmstar)))
    
return, model
end
