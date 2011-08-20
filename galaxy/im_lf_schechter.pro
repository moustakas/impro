;+
; NAME:
;   IM_LF_SCHECHTER()
;
; PURPOSE:
;   Generate a Schechter model given the relevant parameters.
;
; INPUTS: 
;   absmag - absolute magnitudes at which to evaluate the model 
;   phistar - number density at the 'knee' of the LF
;   mstar - absolute magnitude at the 'knee' of the LF
;   alpha - low-mass slope
;
; OPTIONAL INPUTS: 
;   schechter - IM_FIT_LF_SCHECHTER() style structure with PHISTAR,
;     MSTAR,  and ALPHA (takes precedence over PHISTAR, MSTAR, and
;     ALPHA)  
;
; OUTPUTS: 
;   model - absolute magnitude (not luminosity) version of the
;     Schechter model Phi(M)  
; 
; COMMENTS:
;   Note that Phi(M) = ln(10)*L*Phi(L) where L is luminosity. 
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

function im_lf_schechter, absmag, phistar, mstar, alpha, schechter=schechter

    if (n_elements(absmag) eq 0L) then begin
       doc_library, 'im_lf_schechter'
       return, -1
    endif

    if (n_elements(schechter) eq 0) then begin
       if (n_elements(phistar) eq 0) or (n_elements(mstar) eq 0) or $
         (n_elements(alpha) eq 0) then begin
          splog, 'PHISTAR, MSTAR, and ALPHA inputs required'
          return, -1
       endif
       use_phistar = phistar
       use_mstar = mstar
       use_alpha = alpha
    endif else begin
       if (tag_exist(schechter,'phistar') eq 0) or $
         (tag_exist(schechter,'mstar') eq 0) or $
         (tag_exist(schechter,'alpha') eq 0) then begin
          splog, 'Improper SCHECHTER structure!'
          return, -1
       endif
       use_phistar = schechter.phistar
       use_mstar = schechter.mstar
       use_alpha = schechter.alpha
    endelse

    model = 0.4*alog(10)*use_phistar*10.0^(-0.4*(absmag-use_mstar)*(use_alpha+1))*$
      exp(-10.0^(-0.4*(absmag-use_mstar)))

return, model
end
