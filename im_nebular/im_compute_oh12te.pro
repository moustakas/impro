;+
; NAME:
;   IM_COMPUTE_OH12TE()
;
; PURPOSE:
;   Compute electron-temperature abundances.
;
; CALLING SEQUENCE:
;   result = im_compute_oh12te(data)
;
; INPUTS:
;   data - output from COMPUTE_TE() 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   result - results
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;   IDL_FIVEL(), IM_COMPUTE_ERROR()
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 July 5, U of A - written
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

function im_compute_oh12te, data

    nobject = n_elements(data)
    if (nobject eq 0L) then begin
       print, 'Syntax - result = im_compute_oh12te(data)'
       return, -1L
    endif

    oiii_indx = 5L     ; index of the [O III] 5007 transition
    oii_3726_indx = 0L ; [O II] 3726
    oii_3728_indx = 1L ; [O II] 3728
    nii_indx = 5L      ; [N II] 6584

    result = {$
      ZT_oh_pp:            -999.0D, $ ; O++ abundance
;     ZT_oh_pp_lower:      -999.0D, $
;     ZT_oh_pp_upper:      -999.0D, $
      ZT_oh_pp_err:        -999.0D, $
      ZT_oh_p:             -999.0D, $ ; O+ abundance
;     ZT_oh_p_lower:       -999.0D, $
;     ZT_oh_p_upper:       -999.0D, $
      ZT_oh_p_err:         -999.0D, $
      ZT_log12oh:          -999.0D, $ ; total oxygen abundance
;     ZT_log12oh_lower:    -999.0D, $
;     ZT_log12oh_upper:    -999.0D, $
      ZT_log12oh_err:      -999.0D, $
      ZT_nh_p:             -999.0D, $ ; N+ abundance
      ZT_nh_p_err:         -999.0D, $
      ZT_log_no:           -999.0D, $ ; log N/O ratio
      ZT_log_no_err:       -999.0D}
      
    result = replicate(result,nobject)

    if tag_exist(data,'NII_6584') then nitrogen = 1L else nitrogen = 0L

    for iobj = 0L, nobject-1L do begin

       if (data[iobj].ZT_T_oiii gt -900.0) and (data[iobj].oii_3727[1] gt 0.0) and $
         (data[iobj].oiii_5007[1] gt 0.0) then begin
                    
; compute the emissivities of the [O III] lines and the abundance of
; the O++ zone

          fivel = idl_fivel(2,7,temperature=data[iobj].ZT_T_oiii,$
            min_temperature=data[iobj].ZT_T_oiii-data[iobj].ZT_T_oiii_err,$
            max_temperature=data[iobj].ZT_T_oiii+data[iobj].ZT_T_oiii_err,$
            density=data[iobj].ZT_density,/montecarlo,/silent)

          numer = data[iobj].oiii_5007[0]
          numer_err = data[iobj].oiii_5007[1]

          denom = fivel.emissivity[oiii_indx] / fivel.jhb
          denom_err = im_compute_error(fivel.emissivity[oiii_indx],fivel.emissivity_err[oiii_indx],$
            fivel.jhb,fivel.jhb_err,/quotient)
          
          result[iobj].ZT_oh_pp = numer / denom
          result[iobj].ZT_oh_pp_err = im_compute_error(numer,numer_err,denom,denom_err,/quotient)

; compute the emissivities of the [O II] lines and the abundance of
; the O+ zone; we have to predict the intensity of one of the oxygen
; lines -- choose [O II] 3728.794 (assuming no emissivity error)

          fivel = idl_fivel(2,6,temperature=data[iobj].ZT_T_oii,$
            min_temperature=data[iobj].ZT_T_oii-data[iobj].ZT_T_oii_err,$
            max_temperature=data[iobj].ZT_T_oii+data[iobj].ZT_T_oii_err,$
            density=data[iobj].ZT_density,/montecarlo,/silent)
          
          numer = data[iobj].oii_3727[0] / (1.0 + fivel.emissivity[oii_3726_indx]/fivel.emissivity[oii_3728_indx])
          numer_err = data[iobj].oii_3727[1]
          
          denom = fivel.emissivity[oii_3728_indx] / fivel.jhb
          denom_err = im_compute_error(fivel.emissivity[oii_3728_indx],fivel.emissivity_err[oii_3728_indx],$
            fivel.jhb,fivel.jhb_err,/quotient)
          
          result[iobj].ZT_oh_p = numer / denom
          result[iobj].ZT_oh_p_err = im_compute_error(numer,numer_err,denom,denom_err,/quotient)

; compute the total abundance
          sum = result[iobj].ZT_oh_pp + result[iobj].ZT_oh_p
          sum_err = sqrt(result[iobj].ZT_oh_pp_err^2.0 + result[iobj].ZT_oh_p_err^2.0)
          
          result[iobj].ZT_log12oh = 12.0 + alog10(sum)
          result[iobj].ZT_log12oh_err = sum_err/sum/alog(10.0)

; finally compute the N/O ratio assuming N+/O+ = N/O and T(NII) = T(OII)
          if nitrogen then begin

             if (data[iobj].nii_6584[1] gt 0.0) then begin

                fivel = idl_fivel(2,3,temperature=data[iobj].ZT_T_oii,$
                  min_temperature=data[iobj].ZT_T_oii-data[iobj].ZT_T_oii_err,$
                  max_temperature=data[iobj].ZT_T_oii+data[iobj].ZT_T_oii_err,$
                  density=data[iobj].ZT_density,/montecarlo,/silent)

                numer = data[iobj].nii_6584[0]
                numer_err = data[iobj].nii_6584[1]
                
                denom = fivel.emissivity[nii_indx] / fivel.jhb
                denom_err = im_compute_error(fivel.emissivity[nii_indx],fivel.emissivity_err[nii_indx],$
                  fivel.jhb,fivel.jhb_err,/quotient)

                result[iobj].ZT_nh_p = numer / denom
                result[iobj].ZT_nh_p_err = im_compute_error(numer,numer_err,denom,denom_err,/quotient)

                no = result[iobj].ZT_nh_p / result[iobj].ZT_oh_p
                no_err = im_compute_error(result[iobj].ZT_nh_p,result[iobj].ZT_nh_p_err,$
                  result[iobj].ZT_oh_p,result[iobj].ZT_oh_p_err,/quotient)
                
                result[iobj].ZT_log_no = alog10(no)
                result[iobj].ZT_log_no_err = no_err/no/alog(10.0)
                
             endif ; close [N II] condition
          endif ; close NITROGEN condition
       endif ; close T(OIII), [O II] and [O III] conditions
    endfor 

return, result
end
