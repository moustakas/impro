;+
; NAME:
;   IM_COMPUTE_TE()
;
; PURPOSE:
;   Compute the nebular electron temperature based on various
;   temperature-sensitive line-ratios and, if available, the
;   electron density from [S II]. 
;
; INPUTS:
;   data - linefit data structure
;
; OPTIONAL INPUTS:
;   snrcut - S/N cut on the emission lines (default 3)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   result - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 July 5, U of A - written
;   jm06may01uofa - added SNRCUT
;   jm07nov26nyu - rewritten to use IM_TEMDEN()
;   jm07dec11nyu - additional ions included 
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

function im_compute_te, data, nmonte=nmonte, snrcut=snrcut

    nobject = n_elements(data)
    if (nobject eq 0L) then begin
       doc_library, 'im_compute_te'
       return, -1L
    endif

    if (n_elements(snrcut) eq 0L) then snrcut = 1.0

    ocor = 1.0+1.0/(im_branch_ratios()).o_iii

    result = {$
      zt_sii_ratio:          -999.0, $ ; 6716/6731
      zt_sii_ratio_err:      -999.0, $
      zt_oii_dens_ratio:     -999.0, $ ; 3726/3729
      zt_oii_dens_ratio_err: -999.0, $
      zt_oii_temp_ratio:     -999.0, $ ; (3726+3729)/7325 = 3727/7325
      zt_oii_temp_ratio_err: -999.0, $
      zt_nii_ratio:          -999.0, $ ; (6548+6584)/5755
      zt_nii_ratio_err:      -999.0, $
      zt_siii_ratio:         -999.0, $ ; (9069+9532)/6312
      zt_siii_ratio_err:     -999.0, $
      zt_oiii_ratio:         -999.0, $ ; (4959+5007)/4363
      zt_oiii_ratio_err:     -999.0, $

      zt_sii_dens:           -999.0, $ ; [S II] density [cm-3]
      zt_sii_dens_err:       -999.0, $
      zt_oii_dens:           -999.0, $ ; [O II] density [cm-3]
      zt_oii_dens_err:       -999.0, $

      zt_t7325:              -999.0, $ ; T[7325] [K]
      zt_t7325_err:          -999.0, $
      zt_t5755:              -999.0, $ ; T[5755] [K]
      zt_t5755_err:          -999.0, $
      zt_t6312:              -999.0, $ ; T[6312] [K]
      zt_t6312_err:          -999.0, $
      zt_t4363:              -999.0, $ ; T[4363] [K]
      zt_t4363_err:          -999.0, $

      zt_toiii:              -999.0, $ ; final T(O++)=T(Ne++) [K]
      zt_toiii_err:          -999.0, $
      zt_tsiii:              -999.0, $ ; final T(S++)=T(Ar++) [K]
      zt_tsiii_err:          -999.0, $
      zt_toii:               -999.0, $ ; final T(O+)=T(N+)=T(S+) [K]
      zt_toii_err:           -999.0}

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
       
; compute the (density-dependent) [O II] line ratio
       if tag_exist(data,'OII_3726') and tag_exist(data,'OII_3729') then begin

          if (data[iobj].oii_3726[0]/data[iobj].oii_3726[1] gt snrcut) and $
            (data[iobj].oii_3729[0]/data[iobj].oii_3729[1] gt snrcut) then begin
          
             oii_dens_ratio = data[iobj].oii_3726[0]/data[iobj].oii_3729[0]
             oii_dens_ratio_err = im_compute_error(data[iobj].oii_3726[0],data[iobj].oii_3726[1],$
               data[iobj].oii_3729[0],data[iobj].oii_3729[1],/quotient)

             result[iobj].zt_oii_dens_ratio     = oii_dens_ratio
             result[iobj].zt_oii_dens_ratio_err = oii_dens_ratio_err
             
          endif else oii_dens_ratio_err = -2.0
       
       endif else oii_dens_ratio_err = -2.0
       
; compute the (temperature-dependent) [O II] line ratio
       
       if tag_exist(data,'OII_3727') and tag_exist(data,'OII_7325') then begin

          if (data[iobj].oii_3727[0]/data[iobj].oii_3727[1] gt snrcut) and $
            (data[iobj].oii_7325[0]/data[iobj].oii_7325[1] gt snrcut) then begin
          
             oii_temp_ratio = data[iobj].oii_3727[0]/data[iobj].oii_7325[0]
             oii_temp_ratio_err = im_compute_error(data[iobj].oii_3727[0],data[iobj].oii_3727[1],$
               data[iobj].oii_7325[0],data[iobj].oii_7325[1],/quotient)

             result[iobj].zt_oii_temp_ratio     = oii_temp_ratio
             result[iobj].zt_oii_temp_ratio_err = oii_temp_ratio_err
             
          endif else oii_temp_ratio_err = -2.0
       
       endif else oii_temp_ratio_err = -2.0
       
; compute the [N II] line ratio

       if tag_exist(data,'NII_5755') and tag_exist(data,'NII_6548') and tag_exist(data,'NII_6584') then begin

          if (data[iobj].nii_6584[0]/data[iobj].nii_6584[1] gt snrcut) and $
            (data[iobj].nii_6548[0]/data[iobj].nii_6548[1] gt snrcut) and $
            (data[iobj].nii_5755[0]/data[iobj].nii_5755[1] gt snrcut) then begin

             nii     = data[iobj].nii_6548[0]+data[iobj].nii_6584[0]
             nii_err = sqrt(data[iobj].nii_6548[1]^2+data[iobj].nii_6584[1]^2)
             
             nii_ratio = nii/data[iobj].nii_5755[0]
             nii_ratio_err = im_compute_error(nii,nii_err,data[iobj].nii_5755[0],$
               data[iobj].nii_5755[1],/quotient)

             result[iobj].zt_nii_ratio = nii_ratio
             result[iobj].zt_nii_ratio_err = nii_ratio_err
       
          endif else nii_ratio_err = -2.0

       endif else nii_ratio_err = -2.0
       
; compute the [S III] line ratio

       if tag_exist(data,'SIII_6312') and tag_exist(data,'SIII_9069') and tag_exist(data,'SIII_9532') then begin

          if (data[iobj].siii_9532[0]/data[iobj].siii_9532[1] gt snrcut) and $
            (data[iobj].siii_9069[0]/data[iobj].siii_9069[1] gt snrcut) and $
            (data[iobj].siii_6312[0]/data[iobj].siii_6312[1] gt snrcut) then begin

             siii     = data[iobj].siii_9069[0]+data[iobj].siii_9532[0]
             siii_err = sqrt(data[iobj].siii_9069[1]^2+data[iobj].siii_9532[1]^2)
             
             siii_ratio = siii/data[iobj].siii_6312[0]
             siii_ratio_err = im_compute_error(siii,siii_err,data[iobj].siii_6312[0],$
               data[iobj].siii_6312[1],/quotient)

             result[iobj].zt_siii_ratio = siii_ratio
             result[iobj].zt_siii_ratio_err = siii_ratio_err
       
          endif else siii_ratio_err = -2.0

       endif else siii_ratio_err = -2.0
       
; compute the [O III] line ratio

       if tag_exist(data,'OIII_4363') and tag_exist(data,'OIII_4959') and tag_exist(data,'OIII_5007') then begin

          if (data[iobj].oiii_5007[0]/data[iobj].oiii_5007[1] gt snrcut) and $
            (data[iobj].oiii_4959[0]/data[iobj].oiii_4959[1] gt snrcut) and $
            (data[iobj].oiii_4363[0]/data[iobj].oiii_4363[1] gt snrcut) then begin

             oiii     = data[iobj].oiii_4959[0]+data[iobj].oiii_5007[0]
             oiii_err = sqrt(data[iobj].oiii_4959[1]^2+data[iobj].oiii_5007[1]^2)
             
             oiii_ratio = oiii/data[iobj].oiii_4363[0]
             oiii_ratio_err = im_compute_error(oiii,oiii_err,data[iobj].oiii_4363[0],$
               data[iobj].oiii_4363[1],/quotient)

             result[iobj].zt_oiii_ratio = oiii_ratio
             result[iobj].zt_oiii_ratio_err = oiii_ratio_err
       
          endif else oiii_ratio_err = -2.0

       endif else oiii_ratio_err = -2.0
       
; iterate on the density and temperature

;      if (oiii_ratio_err gt 0.0) and (sii_ratio_err gt 0.0) then begin
;     temden = im_temden(['o_iii','s_ii'],[oiii_ratio,sii_ratio],$
;       ratio_err=[oiii_ratio_err,sii_ratio_err],nmonte=nmonte)
;     result[iobj].zt_t4363        = temden[0].temp
;     result[iobj].zt_t4363_err    = temden[0].temp_err
;     result[iobj].zt_sii_dens     = temden[1].dens
;     result[iobj].zt_sii_dens_err = temden[1].dens_err
;      endif
;      if (oiii_ratio_err lt 0.0) and (sii_ratio_err gt 0.0) then begin
;     temden = im_temden('s_ii',sii_ratio,ratio_err=sii_ratio_err,$
;       nmonte=nmonte)
;     result[iobj].zt_sii_dens     = temden[0].dens
;     result[iobj].zt_sii_dens_err = temden[0].dens_err
;      endif
;      if (oiii_ratio_err gt 0.0) and (sii_ratio_err lt 0.0) then begin
;     temden = im_temden('o_iii',oiii_ratio,ratio_err=oiii_ratio_err,$
;       nmonte=nmonte)
;     result[iobj].zt_t4363      = temden[0].temp
;     result[iobj].zt_t4363_err  = temden[0].temp_err
;      endif

       if (sii_ratio_err gt 0.0) then begin
          temden = im_temden('s_ii',sii_ratio,$
            ratio_err=sii_ratio_err,nmonte=nmonte)
          result[iobj].zt_sii_dens     = temden[0].dens
          result[iobj].zt_sii_dens_err = temden[0].dens_err
       endif
       if (oii_dens_ratio_err gt 0.0) then begin
          temden = im_temden('o_ii_dens',oii_dens_ratio,$
            ratio_err=oii_dens_ratio_err,nmonte=nmonte)
          result[iobj].zt_oii_dens     = temden[0].dens
          result[iobj].zt_oii_dens_err = temden[0].dens_err
       endif
       if (oii_temp_ratio_err gt 0.0) then begin
          temden = im_temden('o_ii_temp',oii_temp_ratio,$
            ratio_err=oii_temp_ratio_err,nmonte=nmonte)
          result[iobj].zt_t7325      = temden[0].temp
          result[iobj].zt_t7325_err  = temden[0].temp_err
       endif
       if (nii_ratio_err gt 0.0) then begin
          temden = im_temden('n_ii',nii_ratio,$
            ratio_err=nii_ratio_err,nmonte=nmonte)
          result[iobj].zt_t5755      = temden[0].temp
          result[iobj].zt_t5755_err  = temden[0].temp_err
       endif
       if (siii_ratio_err gt 0.0) then begin
          temden = im_temden('s_iii',siii_ratio,$
            ratio_err=siii_ratio_err,nmonte=nmonte)
          result[iobj].zt_t6312      = temden[0].temp
          result[iobj].zt_t6312_err  = temden[0].temp_err
       endif
       if (oiii_ratio_err gt 0.0) then begin
          temden = im_temden('o_iii',oiii_ratio,$
            ratio_err=oiii_ratio_err,nmonte=nmonte)
          result[iobj].zt_t4363      = temden[0].temp
          result[iobj].zt_t4363_err  = temden[0].temp_err
       endif

    endfor
    
; assign final ionic temperatures here, mostly following Bresolin et
; al. (2005): 

; T[OIII] - (1) if T4363 exists, then set T[OIII]=T4363; (2) if T4363
;       does not exist *and both* T6312 and T5755 exist, then
;       invert the Garnett (1992) relations to get two estimates
;       of T[OIII], combining them using a weighted mean (or an
;       unweighted mean if NMONTE=0); (3) if only one of T6312 or
;       T5755 exist then do the same as in step (2) with just the
;       single value

    tt = where((result.zt_toiii lt -900.0) and (result.zt_t4363 gt -900.0),ntt)
    if (ntt ne 0L) then begin
       result[tt].zt_toiii     = result[tt].zt_t4363
       result[tt].zt_toiii_err = result[tt].zt_t4363_err
    endif

    tt = where((result.zt_toiii lt -900.0) and (result.zt_t4363 lt -900.0) and $
      (result.zt_t6312 gt -900.0) and (result.zt_t5755 gt -900.0),ntt)
    if (ntt ne 0L) then begin
       toiii_6312 = (result[tt].zt_t6312-1700.0)/0.83 ; invert Garnett's relations
       toiii_5755 = (result[tt].zt_t5755-3000.0)/0.70
       toiii_6312_err = (result[tt].zt_t6312_err/0.83)>0.0
       toiii_5755_err = (result[tt].zt_t5755_err/0.70)>0.0
       for ii = 0L, ntt-1L do begin
          if (toiii_6312_err[ii] gt 0.0) and (toiii_5755_err[ii] gt 0.0) then begin ; weighted mean
             result[tt[ii]].zt_toiii     = total([toiii_6312[ii],toiii_5755[ii]]/$
               [toiii_6312_err[ii],toiii_5755_err[ii]]^2.0) / $
               total(1.0/[toiii_6312_err[ii],toiii_5755_err[ii]]^2)
             result[tt[ii]].zt_toiii_err = 1.0/sqrt(total(1.0/[toiii_6312_err[ii],toiii_5755_err[ii]]^2))
          endif else begin ; straight average
             result[tt[ii]].zt_toiii     = djs_mean([toiii_6312[ii],toiii_5755[ii]])
             result[tt[ii]].zt_toiii_err = djsig([toiii_6312[ii],toiii_5755[ii]])
          endelse
       endfor
    endif

    tt = where((result.zt_toiii lt -900.0) and (result.zt_t4363 lt -900.0) and $
      (result.zt_t6312 lt -900.0) and (result.zt_t5755 gt -900.0),ntt)
    if (ntt ne 0L) then begin
       result[tt].zt_toiii     = (result[tt].zt_t5755-3000.0)/0.7 ; invert Garnett's relations
       result[tt].zt_toiii_err = result[tt].zt_t5755_err/0.7
    endif

    tt = where((result.zt_toiii lt -900.0) and (result.zt_t4363 lt -900.0) and $
      (result.zt_t6312 gt -900.0) and (result.zt_t5755 lt -900.0),ntt)
    if (ntt ne 0L) then begin
       result[tt].zt_toiii     = (result[tt].zt_t6312-1700.0)/0.83 ; invert Garnett's relations
       result[tt].zt_toiii_err = result[tt].zt_t6312_err/0.83
    endif

; predict the temperature of the O+ zone from the temperature of the
; O++ zone; adopt the Garnett (1992) or the Stasinska (1990)
; photoionization relation

     toii = im_predict_toii(result.zt_toiii,toii_err=toii_err,$
       toiii_err=result.zt_toiii_err)
     result.zt_toii = toii
     result.zt_toii_err = toii_err

return, result
end
