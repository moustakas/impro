;+
; NAME:
;       IM_IZOTOV_TE()
;
; PURPOSE:
;       Compute the nebular electron temperature from the [O III] flux 
;       ratio using the equations in Izotov et al. (2006).
;
; CALLING SEQUENCE:
;       result = im_compute_te(data,snrcut=)
;
; INPUTS:
;       data - linefit data structure
;
; OPTIONAL INPUTS:
;       snrcut - S/N cut on the emission lines (default 3)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       result - results
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       SPLOG, IDL_FIVEL(), IM_COMPUTE_ERROR(), IM_PREDICT_TOII() 
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 July 5, U of A - written
;       jm06may01uofa - added SNRCUT
;
; Copyright (C) 2004, 2006, John Moustakas
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

function im_izotov_oh12te, data, snrcut=snrcut

    nobject = n_elements(data)
    if (nobject eq 0L) then begin
       print, 'Syntax - result = im_izotov_oh12te(data,snrcut=)'
       return, -1L
    endif

    if (n_elements(snrcut) eq 0L) then snrcut = 3.0
    
    result = {$
      ff_r3:                        -999.0, $ ; [O III] 4959,5007 / H-beta
      ff_r3_err:                    -999.0, $
      

      ZT_oiii_ratio:       -999.0D, $ ; 4959 + 5007 / 4363
      ZT_oiii_ratio_err:   -999.0D, $
      ZT_sii_ratio:        -999.0D, $ ; 6716 / 6731
      ZT_sii_ratio_err:    -999.0D, $
      ZT_density:          -999.0D, $ ; density [cm-3]
      ZT_density_err:      -999.0D, $
      ZT_T_oiii:           -999.0D, $ ; T(O++) [K]
      ZT_T_oiii_lower:     -999.0D, $
      ZT_T_oiii_upper:     -999.0D, $
      ZT_T_oiii_err:       -999.0D, $
      ZT_T_oii:            -999.0D, $ ; T(O+) [K]
      ZT_T_oii_lower:      -999.0D, $
      ZT_T_oii_upper:      -999.0D, $
      ZT_T_oii_err:        -999.0D}
    result = replicate(result,nobject)

    if (tag_exist(data,'OIII_4363') eq 0L) or (tag_exist(data,'OIII_5007') eq 0L) then begin
       splog, 'Insufficient information to compute the electron temperature.'
       return, result
    endif

    if tag_exist(data,'SII_6716') and tag_exist(data,'SII_6731') then sulfur = 1L else sulfur = 0L

; loop on every object

    niter = -1L ; starting value
    
    for iobj = 0L, nobject-1L do begin

; compute the [S II] line ratio
       
       if sulfur then begin

          if (data[iobj].sii_6716[0]/data[iobj].sii_6716[1] gt snrcut) and $
            (data[iobj].sii_6731[0]/data[iobj].sii_6731[1] gt snrcut) then begin
          
             sii_ratio = data[iobj].sii_6716[0]/data[iobj].sii_6731[0]
             sii_ratio_err = im_compute_error(data[iobj].sii_6716[0],data[iobj].sii_6716[1],$
               data[iobj].sii_6731[0],data[iobj].sii_6731[1],/quotient)

             result[iobj].ZT_sii_ratio = sii_ratio
             result[iobj].ZT_sii_ratio_err = sii_ratio_err
             
          endif else sii_ratio_err = -2.0
       
       endif else sii_ratio_err = -2.0
       
; compute the [O III] line ratio

       if (data[iobj].oiii_5007[0]/data[iobj].oiii_5007[1] gt snrcut) and $
         (data[iobj].oiii_4363[0]/data[iobj].oiii_4363[1] gt snrcut) then begin

          oiii_ratio = data[iobj].oiii_5007[0]*ocor/data[iobj].oiii_4363[0]
          oiii_ratio_err = im_compute_error(data[iobj].oiii_5007[0]*ocor,data[iobj].oiii_5007[1],$
            data[iobj].oiii_4363[0],data[iobj].oiii_4363[1],/quotient)

          result[iobj].ZT_oiii_ratio = oiii_ratio
          result[iobj].ZT_oiii_ratio_err = oiii_ratio_err
       
       endif else oiii_ratio_err = -2.0
       
; iterate on the density and temperature

       if (sii_ratio_err gt 0.0) and (oiii_ratio_err gt 0.0) then niter = 2L

       if ((sii_ratio_err lt 0.0) and (oiii_ratio_err gt 0.0)) or $
         ((sii_ratio_err gt 0.0) and (oiii_ratio_err lt 0.0)) then niter = 1L
       
       for iter = 0L, niter-1L do begin
       
          if (sii_ratio_err gt 0.0) then begin
             fivel = idl_fivel(1,11,lineratio=sii_ratio,err_lineratio=sii_ratio_err,$
               temperature=temperature,density=density,/montecarlo,/silent)
             result[iobj].ZT_density = fivel.density
             result[iobj].ZT_density_err = djs_mean([fivel.density_lower,fivel.density_upper])
             density = result[iobj].ZT_density
          endif else density = 1D2

          if (oiii_ratio_err gt 0.0) then begin
             fivel = idl_fivel(1,6,lineratio=oiii_ratio,err_lineratio=oiii_ratio_err,$
               temperature=temperature,density=density,/montecarlo,/silent)
             result[iobj].ZT_T_oiii = fivel.temperature
             result[iobj].ZT_T_oiii_lower = fivel.temperature_lower
             result[iobj].ZT_T_oiii_upper = fivel.temperature_upper
             result[iobj].ZT_T_oiii_err = djs_mean([fivel.temperature_lower,fivel.temperature_upper])
             temperature = result[iobj].ZT_T_oiii
          endif else temperature = 1D4

       endfor

; if no density calculation was possible then assume low density

       if (sii_ratio_err lt 0.0) then result[iobj].ZT_density = 1D2
       
; adopt the Garnett (1992) or the Stasinska (1990) photoionization
; relation to predict the temperature of the O+ zone from the
; temperature of the O++ zone
       
       if (oiii_ratio_err gt 0.0) then begin

          T_oii = im_predict_toii(result[iobj].ZT_T_oiii,toiii_lower=result[iobj].ZT_T_oiii_lower,$
            toiii_upper=result[iobj].ZT_T_oiii_upper,stasinska=stasinska)

          result[iobj].ZT_T_oii       = T_oii.ZT_T_oii
          result[iobj].ZT_T_oii_lower = T_oii.ZT_T_oii_lower
          result[iobj].ZT_T_oii_upper = T_oii.ZT_T_oii_upper
          result[iobj].ZT_T_oii_err   = T_oii.ZT_T_oii_err

       endif
       
    endfor

return, result
end
