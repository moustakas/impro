;+
; NAME:
;       IM_COMPUTE_FIVEL_TE()
;
; PURPOSE:
;       Compute the nebular electron temperature from the [O III] flux 
;       ratio and, if available, the electron density from [S II].
;
; INPUTS:
;       data - linefit data structure
;
; OPTIONAL INPUTS:
;       snrcut - S/N cut on the emission lines (default 3)
;
; KEYWORD PARAMETERS:
;       fivel     - use the old FIVEL code rather than IM_TEMDEN() 
;       stasinska - use the Stasinska relation (default is to use
;                   Garnett) 
;
; OUTPUTS:
;       result - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 July 5, U of A - written
;       jm06may01uofa - added SNRCUT
;       jm07nov26nyu - rewritten to use IM_TEMDEN()
;
; Copyright (C) 2004, 2006-2007, John Moustakas
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

function im_compute_fivel_te, data, nmonte=nmonte, snrcut=snrcut, $
  fivel=fivel, stasinska=stasinska

    nobject = n_elements(data)
    if (nobject eq 0L) then begin
       doc_library, 'im_compute_fivel_te'
       return, -1L
    endif

    if (n_elements(snrcut) eq 0L) then snrcut = 1.0

    ocor = 1.0+1.0/(im_branch_ratios()).o_iii

    result = {$
      zt_oiii_ratio:       -999.0, $ ; 4959 + 5007 / 4363
      zt_oiii_ratio_err:   -999.0, $
      zt_sii_ratio:        -999.0, $ ; 6716 / 6731
      zt_sii_ratio_err:    -999.0, $
      zt_density:          -999.0, $ ; density [cm-3]
      zt_density_err:      -999.0, $
      zt_toiii:           -999.0, $ ; t(O++) [K]
      zt_toiii_err:       -999.0, $
      zt_toii:            -999.0, $ ; t(O+) [K]
      zt_toii_err:        -999.0, $
      zt_s_oiii:           -999.0, $ ; t(O+) [K]
      zt_s_oiii_err:       -999.0}
    result = replicate(result,nobject)

    for iobj = 0L, nobject-1L do begin

; compute the [S II] line ratio
       
       if tag_exist(data,'SII_6716') and tag_exist(data,'SII_6731') then begin

          if (data[iobj].sii_6716[0]/data[iobj].sii_6716[1] gt snrcut) and $
            (data[iobj].sii_6731[0]/data[iobj].sii_6731[1] gt snrcut) then begin
          
             sii_ratio = data[iobj].sii_6716[0]/data[iobj].sii_6731[0]
             sii_ratio_err = im_compute_error(data[iobj].sii_6716[0],data[iobj].sii_6716[1],$
               data[iobj].sii_6731[0],data[iobj].sii_6731[1],/quotient)

             result[iobj].zt_sii_ratio = sii_ratio
             result[iobj].zt_sii_ratio_err = sii_ratio_err
             
          endif else sii_ratio_err = -2.0
       
       endif else sii_ratio_err = -2.0
       
; compute the [O III] line ratio

       if tag_exist(data,'OIII_4363') and tag_exist(data,'OIII_5007') then begin

          if (data[iobj].oiii_5007[0]/data[iobj].oiii_5007[1] gt snrcut) and $
            (data[iobj].oiii_4363[0]/data[iobj].oiii_4363[1] gt snrcut) then begin

             oiii_ratio = data[iobj].oiii_5007[0]*ocor/data[iobj].oiii_4363[0]
             oiii_ratio_err = im_compute_error(data[iobj].oiii_5007[0]*ocor,data[iobj].oiii_5007[1],$
               data[iobj].oiii_4363[0],data[iobj].oiii_4363[1],/quotient)

             result[iobj].zt_oiii_ratio = oiii_ratio
             result[iobj].zt_oiii_ratio_err = oiii_ratio_err
       
          endif else oiii_ratio_err = -2.0

       endif else oiii_ratio_err = -2.0
       
; iterate on the density and temperature

       if keyword_set(fivel) then begin ; use FIVEL
          
          for iter = 0L, niter-1L do begin
             
             if (sii_ratio_err gt 0.0) then begin
                fivel = idl_fivel(1,11,lineratio=sii_ratio,err_lineratio=sii_ratio_err,$
                  temperature=temperature,density=density,/montecarlo,/silent)
                result[iobj].zt_density = fivel.density
                density = result[iobj].zt_density
             endif else density = 1D2

             if (oiii_ratio_err gt 0.0) then begin
                fivel = idl_fivel(1,6,lineratio=oiii_ratio,err_lineratio=oiii_ratio_err,$
                  temperature=temperature,density=density,/montecarlo,/silent)
                result[iobj].zt_toiii = fivel.temperature
                temperature = result[iobj].zt_toiii
             endif else temperature = 1D4

          endfor

       endif else begin ; use IM_TEMDEN()

          if (oiii_ratio_err gt 0.0) and (sii_ratio_err gt 0.0) then begin
             temden = im_temden(['o_iii','s_ii'],[oiii_ratio,sii_ratio],$
               ratio_err=[oiii_ratio_err,sii_ratio_err],nmonte=nmonte)
             result[iobj].zt_toiii      = temden[0].temp
             result[iobj].zt_toiii_err  = temden[0].temp_err
             result[iobj].zt_density     = temden[1].dens
             result[iobj].zt_density_err = temden[1].dens_err
          endif
          if (oiii_ratio_err lt 0.0) and (sii_ratio_err gt 0.0) then begin
             temden = im_temden('s_ii',sii_ratio,ratio_err=sii_ratio_err,$
               nmonte=nmonte)
             result[iobj].zt_density     = temden[0].dens
             result[iobj].zt_density_err = temden[0].dens_err
          endif
          if (oiii_ratio_err gt 0.0) and (sii_ratio_err lt 0.0) then begin
             temden = im_temden('o_iii',oiii_ratio,ratio_err=oiii_ratio_err,$
               nmonte=nmonte)
             result[iobj].zt_toiii      = temden[0].temp
             result[iobj].zt_toiii_err  = temden[0].temp_err
          endif

       endelse
          
; predict the temperature of the O+ zone from the temperature of the
; O++ zone; adopt the Garnett (1992) or the Stasinska (1990)
; photoionization relation
       
       if (result[iobj].zt_toiii gt -900.0) then begin

          toii = im_predict_toii(result[iobj].zt_toiii,$
            toiii_err=result[iobj].zt_toiii_err,$
            toii_err=toii_err)
          result[iobj].zt_toii       = toii
          result[iobj].zt_toii_err   = toii_err

       endif
       
    endfor

return, result
end
