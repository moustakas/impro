;+
; NAME:
;       IM_DIRECT_ABUNDANCE()
;
; PURPOSE:
;       Given the relevant physical conditions (temperature and
;       density), and the appropriate emission-line ratio (relative to
;       H-beta), compute the total elemental abundance.
;
; INPUTS:
;       data - ispec1d type structure with the necessary emission-line
;              ratios (assumed relative to H-beta) [see
;              INIT_HII_LINEFIT_STRUCTURE] 
;
; OPTIONAL INPUTS:
;       oii_ion  - (default 3727)
;       oiii_ion - (default 5007)
;       nii_ion  - (default 6584)
;       nmonte   - number of Monte Carlo iterations; if NMONTE=0
;                  then don't do the Monte Carlo bit
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       result - data structure with lots of goodies
;
; COMMENTS:
;       Currently only computes the oxygen abundance and the N/O
;       ratio. 
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Dec 4, NYU - written
;
; Copyright (C) 2007, John Moustakas
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

function im_direct_abundance, data, snrcut=snrcut, oii_temp=oii_temp, oiii_temp=oiii_temp, $
  err_oii_temp=err_oii_temp, err_oiii_temp=err_oiii_temp, dens=dens, err_dens=err_dens, $
  oii_ion=oii_ion, oiii_ion=oiii_ion, nii_ion=nii_ion, nmonte=nmonte

    nobj = n_elements(data)    
    if (nobj eq 0L) then begin
       doc_library, 'im_direct_abundance'
       return, -1L
    endif

    if (n_elements(oii_temp) eq 0L) then oii_temp = replicate(1E4,nobj)
    if (n_elements(oiii_temp) eq 0L) then oiii_temp = replicate(1E4,nobj)
    if (n_elements(dens) eq 0L) then dens = replicate(1E2,nobj)

    if (n_elements(err_oii_temp) eq 0L) then err_oii_temp = oii_temp*0.0
    if (n_elements(err_oiii_temp) eq 0L) then err_oiii_temp = oiii_temp*0.0
    if (n_elements(err_dens) eq 0L) then err_dens = dens*0.0

    if (n_elements(oii_temp) ne nobj) then begin
       splog, 'Dimensions of DATA and OII_TEMP must match.'
       return, -1L
    endif
    if (n_elements(oiii_temp) ne nobj) then begin
       splog, 'Dimensions of DATA and OIII_TEMP must match.'
       return, -1L
    endif
    if (n_elements(dens) ne nobj) then begin
       splog, 'Dimensions of DATA and DENS must match.'
       return, -1L
    endif
    
    if (n_elements(err_oii_temp) ne nobj) then begin
       splog, 'Dimensions of DATA and ERR_OII_TEMP must match.'
       return, -1L
    endif
    if (n_elements(err_oiii_temp) ne nobj) then begin
       splog, 'Dimensions of DATA and ERR_OIII_TEMP must match.'
       return, -1L
    endif
    if (n_elements(err_dens) ne nobj) then begin
       splog, 'Dimensions of DATA and ERR_DENS must match.'
       return, -1L
    endif
    
    if (n_elements(snrcut) eq 0L) then snrcut = 1.0

    if (n_elements(oii_ion) eq 0L) then oii_ion = '3727'
    if (n_elements(oiii_ion) eq 0L) then oiii_ion = '5007'
    if (n_elements(nii_ion) eq 0L) then nii_ion = '6584'

    if (n_elements(nmonte) eq 0L) then nmonte = 500L

; initialize the output data structure

    result = {$
      ZT_oh_p:             -999.0, $ ; O+ abundance
      ZT_oh_p_err:         -999.0, $
      ZT_oh_pp:            -999.0, $ ; O++ abundance
      ZT_oh_pp_err:        -999.0, $
      ZT_log12oh_te:       -999.0, $ ; total oxygen abundance
      ZT_log12oh_te_err:   -999.0, $
      ZT_nh_p:             -999.0, $ ; N+ abundance
      ZT_nh_p_err:         -999.0, $
      ZT_log_no:           -999.0, $ ; log N/O ratio
      ZT_log_no_err:       -999.0}
    result = replicate(result,nobj)

; check which line-ratio data are available, and store the necessary
; fluxes
    
    datatags = tag_names(data[0])
    oii_indx = where(strmatch(datatags,'*'+oii_ion+'*',/fold),noii_indx)
    oiii_indx = where(strmatch(datatags,'*'+oiii_ion+'*',/fold),noiii_indx)
    nii_indx = where(strmatch(datatags,'*'+oii_ion+'*',/fold),nii_indx)

; compute O+/H    
    if (noii_indx eq 1L) then begin
       good = where((((data.(oii_indx))[0,*])/((data.(oii_indx))[1,*]) gt snrcut) and $
         (oii_temp gt 0.0) and (dens gt 0.0),ngood)
       if (ngood ne 0L) then begin
          
          oii_flux = reform((data[good].(oii_indx))[0,*])
          oii_ferr = reform((data[good].(oii_indx))[1,*])

          for ii = 0L, ngood-1L do begin
             
             level = im_nlevel('o_ii',dens=dens[good[ii]],temp=oii_temp[good[ii]])
             case oii_ion of
                '3727': oii_emissivity = level.emissivity[2,0]+level.emissivity[1,0]     ; 3726+3729
                '7325': oii_emissivity = (level.emissivity[4,1]+level.emissivity[3,1])+$ ; 7319,7320+7330,7331
                  (level.emissivity[4,2]+level.emissivity[3,2])                          ;   
                else: message, '[O II] ion '+oii_ion+' not supported.'
             endcase

             result[good[ii]].zt_oh_p = oii_flux[ii]/oii_emissivity
          endfor
       endif
    endif

; compute O++/H    
    
    if (noiii_indx eq 1L) then begin

       good = where((((data.(oiii_indx))[0,*])/((data.(oiii_indx))[1,*]) gt snrcut) and $
         (oiii_temp gt 0.0) and (dens gt 0.0),ngood)
       if (ngood ne 0L) then begin

          oiii_flux = reform((data[good].(oiii_indx))[0,*])
          oiii_ferr = reform((data[good].(oiii_indx))[1,*])

          for ii = 0L, ngood-1L do begin

             level = im_nlevel('o_iii',dens=dens[good[ii]],temp=oiii_temp[good[ii]])
             case oiii_ion of
                '4959': oiii_emissivity = level.emissivity[3,1] ; 4959
                '5007': oiii_emissivity = level.emissivity[3,2] ; 5007
                else: message, '[O III] ion '+oiii_ion+' not supported.'
             endcase

             result[good[ii]].zt_oh_pp = oiii_flux[ii]/oiii_emissivity

          endfor
             
       endif
          
    endif

; compute the total oxygen abundance    
    
    total_oh = where((result.zt_oh_p gt -900.0) and (result.zt_oh_pp gt -900.0),ntotal_oh)
    if (ntotal_oh ne 0L) then begin
       result[total_oh].zt_log12oh_te = 12.0 + alog10(result[total_oh].zt_oh_p + result[total_oh].zt_oh_pp)
    endif

; compute the uncertainties on the abundances using a Monte Carlo
; method

    if (nmonte gt 0L) then begin
       
; perturb the line-fluxes, temperature, and density by the
; corresponding errors        

       oii_temp_monte = fltarr(nobj,nmonte)-999.0
       oiii_temp_monte = fltarr(nobj,nmonte)-999.0
       dens_monte = (cmreplicate(dens,nmonte) + $
         randomn(seed,nobj,nmonte)*cmreplicate(err_dens,nmonte))>1.0 ; make sure the density is never -999!!!
       
       data_monte = cmreplicate(data,nmonte)
       if (noii_indx ne 0L) then begin
          good = where(((data.(oii_indx))[1,*] gt 0.0) and (oii_temp gt 0.0) and (err_oii_temp ge 0.0),ngood)
          if (ngood ne 0L) then begin
             data_monte[good,*].(oii_indx)[0,*] = cmreplicate(data[good].(oii_indx)[0,*],nmonte) + $
               randomn(seed,ngood,nmonte)*cmreplicate(data[good].(oii_indx)[1,*],nmonte)
             oii_temp_monte[good,*] = cmreplicate(oii_temp[good],nmonte) + $
               randomn(seed,ngood,nmonte)*cmreplicate(err_oii_temp[good],nmonte)
          endif
       endif

       if (noiii_indx ne 0L) then begin
          good = where(((data.(oiii_indx))[1,*] gt 0.0) and (oiii_temp gt 0.0) and (err_oiii_temp ge 0.0),ngood)
          if (ngood ne 0L) then begin
             data_monte[good,*].(oiii_indx)[0,*] = cmreplicate(data[good].(oiii_indx)[0,*],nmonte) + $
               randomn(seed,ngood,nmonte)*cmreplicate(data[good].(oiii_indx)[1,*],nmonte)
             oiii_temp_monte[good,*] = cmreplicate(oiii_temp[good],nmonte) + $
               randomn(seed,ngood,nmonte)*cmreplicate(err_oiii_temp[good],nmonte)
          endif
       endif

; now loop on each Monte Carlo iteration       
       
       for imonte = 0L, nmonte-1L do begin
          resmonte1 = im_direct_abundance(data_monte[*,imonte],oii_ion=oii_ion,oiii_ion=oiii_ion,$
            oii_temp=oii_temp_monte[*,imonte],oiii_temp=oiii_temp_monte[*,imonte],dens=dens_monte[*,imonte],$
            nmonte=0L)
          if (imonte eq 0L) then resmonte = resmonte1 else $
            resmonte = [[resmonte],[resmonte1]]
       endfor

       for iobj = 0L, nobj-1L do begin
          good = where((resmonte[iobj,*].zt_oh_p gt 0.0),ngood)
          if (ngood ne 0L) then begin
             result[iobj].zt_oh_p_err = djsig(resmonte[iobj,good].zt_oh_p)
          endif
          good = where((resmonte[iobj,*].zt_oh_pp gt 0.0),ngood)
          if (ngood ne 0L) then begin
             result[iobj].zt_oh_pp_err = djsig(resmonte[iobj,good].zt_oh_pp)
          endif
          good = where((resmonte[iobj,*].zt_log12oh_te gt 0.0),ngood)
          if (ngood ne 0L) then begin
             result[iobj].zt_log12oh_te_err = djsig(resmonte[iobj,good].zt_log12oh_te)
          endif
       endfor
    endif

return, result
end
