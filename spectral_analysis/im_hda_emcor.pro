;+
; NAME:
;       IM_HDA_EMCOR()
;
; PURPOSE:
;       Given an H-alpha flux and a Ha/Hb Balmer decrement, compute
;       the appropriate correction to apply to the Lick H-delta_A
;       index to account for emission-line filling.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
;
; KEYWORD PARAMETERS: 
;
;
; OUTPUTS: 
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;
; Copyright (C) 2008-May-20, John Moustakas, NYU
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

function im_hda_emcor, linenodust

    nspec = n_elements(linenodust) 
    if (nspec eq 0L) then begin
       splog, 'Does not compute.'
       return, linenodust
    endif

; NEEDS ERROR CHECKING!    
    
    kl = k_lambda(4101.734,_extra=extra)
    R_HaHd = 11.1 ; intrinsic Ha/Hd ratio

    hdcor = {h_delta_dust_predict: [0.0,-2.0], $
      h_delta_nodust_predict: [0.0,-2.0], $
      lick_hd_a_cor: [0.0,-2.0]}
    hdcor = replicate(hdcor,nspec)

; predict the amount of H-delta *emission* given the H-alpha flux and
; the measured Balmer decrement
    
    good = where((linenodust.h_alpha[1] gt 0.0) and $
      (linenodust.ebv_hahb_err gt 0.0) and $
      (linenodust.h_delta_continuum[1] gt 0.0) and $
      (linenodust.h_delta_dust_predict[1] gt 0.0),ngood)
    if (ngood ne 0L) then begin

       ebv = linenodust[good].ebv_hahb
       ebv_err = linenodust[good].ebv_hahb_err

       flux = linenodust[good].h_alpha[0] / R_HaHd
       ferr = linenodust[good].h_alpha[1] / R_HaHd
       
       hdcor[good].h_delta_nodust_predict[0] = flux
       hdcor[good].h_delta_nodust_predict[1] = ferr

       hdcor[good].h_delta_dust_predict[0] = flux * 10^(-0.4*ebv*kl)
       hdcor[good].h_delta_dust_predict[1] = sqrt( (ferr*10.0^(0.4*ebv*kl))^2.0 + $
         (flux*0.4*kl*alog(10)*10.0^(0.4*ebv*kl)*ebv_err)^2.0 )

; now correct the LICK_HD_A index for emission-line contamination 

       hd_dust_predict_ew = linenodust[good].h_delta_dust_predict[0]/$
         linenodust[good].h_delta_continuum[0]
       hd_dust_predict_ew_err = im_compute_error(linenodust[good].h_delta_dust_predict[0],$
         linenodust[good].h_delta_dust_predict[1],linenodust[good].h_delta_continuum[0],$
         linenodust[good].h_delta_continuum[1],/quotient)

; *add* the correction to make the absorption indices bigger
       
       hdcor[good].lick_hd_a_cor[0] = linenodust[good].lick_hd_a[0] + hd_dust_predict_ew
       hdcor[good].lick_hd_a_cor[1] = sqrt(linenodust[good].lick_hd_a[1]^2 + $
         hd_dust_predict_ew_err^2)

    endif

return, hdcor
end
