;+
; NAME:
;   IM_TEMDEN()
;
; PURPOSE:
;   Given one or more appropriate line-ratios, compute the
;   nebular electron temperature and/or density.
;
; INPUTS:
;   ionlist   - string array of density- and temperature-sensitive 
;               ion names corresponding to RATIO
;   ratio - density- and/or temperature-sensitive line-ratios
;
; OPTIONAL INPUTS:
;   ratio_err - uncertainty on RATIO
;   frac_tol  - minimum fractional tolerance to consider a good
;               measurement (see below)
;   nmonte    - number of Monte Carlo iterations (if RATIO_ERR is
;               passed); if NMONTE=0 then don't do the Monte
;               Carlo part
;   temp_guess - initial guess at the temperature
;   max_temp   - maximum temperature to consider (default 25,000 K)
;   min_temp   - minimum temperature to consider (default 5,000 K)
;   step_temp  - temperature step size
;   max_dens   - maximum density to consider (default 10,000 cm^-3)
;   min_dens   - minimum density to consider (default 1 cm^-3)
;   step_dens  - density step size
;
; KEYWORD PARAMETERS:
;   noiter - in the special case that [O III] is passed as the
;            temperature diagnostic and [S II] is passed as the
;            density diagnostic, the default is to iterate the
;            solution so that the values converge
;
; OUTPUTS:
;   result - data structure with lots of goodies
;
; COMMENTS:
;   Temperature-sensitive line-ratios:
;      O_III = (4959+5007)/4363
; 
;   Density-sensitive line-ratios:
;      S_II = 6716/6731
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2007 Nov 25, NYU - written
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

function im_temden, ionlist, ratio, ratio_err=ratio_err, $
  frac_tol=frac_tol, nmonte=nmonte, temp_guess=temp_guess1, $
  dens_guess=dens_guess1, noiter=noiter

    common temden_table, temden_table
    
    nion = n_elements(ionlist)
    nratio = n_elements(ratio)
    
    if (nion eq 0L) or (nratio eq 0L) then begin
       doc_library, 'im_temden'
       return, -1L
    endif

    if (n_elements(ratio_err) ne 0L) then begin
       if (nratio ne n_elements(ratio_err)) then begin
          print, 'Dimensions of RATIO and RATIO_ERR must agree.'
          return, -1L
       endif
    endif

; if NION=1 and NRATIO>1 then call this routine recursively

    if (nion eq 1L) and (nratio gt 1L) then begin
       for ir = 0L, nratio-1L do begin
          if (n_elements(ratio_err) ne 0L) then thisratio_err = ratio_err[ir]
          result1 = im_temden(ionlist[0],ratio[ir],ratio_err=thisratio_err,$
            frac_tol=frac_tol,nmonte=nmonte,temp_guess=temp_guess1,$
            dens_guess=dens_guess1,noiter=noiter)
          if (ir eq 0L) then result = result1 else result = [result,result1]
       endfor
       return, reform(result)
    endif
    
; read the line-ratio look-up table, but do this just once, unless it
; has changed on disk; see WRITE_TEMDEN_LOOKUP_TABLE for more details 

    if (n_elements(temden_table) eq 0L) then begin
       temden_table_file = getenv('IMPRO_DIR')+'/etc/temden_table.fits'
       if (file_test(temden_table_file,/regular) eq 0L) then begin
          splog, 'TEMDEN lookup table '+temden_table_file+' not found.'
          return, -1L
       endif
       temden_table = mrdfits(temden_table_file,1,/silent)
    endif
    
    ngrid_temp = n_elements(temden_table.grid_temp)
    ngrid_dens = n_elements(temden_table.grid_dens)

; initial guesses

    if (n_elements(dens_guess1) eq 0L) then dens_guess = 1D2 else dens_guess = dens_guess1 ; [cm^-3]
    if (n_elements(temp_guess1) eq 0L) then temp_guess = 1D4 else temp_guess = temp_guess1 ; [K]
    
    if (n_elements(nmonte) eq 0L) then nmonte = 500L
    if (n_elements(frac_tol) eq 0L) then frac_tol = 1D-4

; initialize the output data structure; when the uncertainty is given,
; compute the [5,16,25,50,75,84,95] percentiles of the Monte Carlo
; distribution 

    quant = [0.05,0.16,0.25,0.50,0.75,0.84,0.95]
    nquant = n_elements(quant)
    
    result = {$
      ion:                        '', $
      ratio:                  -999.0, $
      ratio_err:              -999.0, $
      model_ratio:            -999.0, $
      ratio_frac_err:         -999.0, $
      temp_ion:                   0B, $ ; temperature-sensitive ion (Boolean) 
      temp_edge:                  0B, $ ; implied temperature is on the edge of the grid
      temp:                   -999.0, $
      temp_err:               -999.0, $
      temp_quant:      quant*0.0-1.0, $ 
      dens_ion:                    0, $ ; density-sensitive ion (Boolean) 
      dens_edge:                   0, $ ; implied density is on the edge of the grid
      dens:                   -999.0, $
      dens_err:               -999.0, $
      dens_quant:      quant*0.0-1.0}
    result = replicate(result,nion)

    result.ion   = ionlist
    result.temp  = temp_guess
    result.dens  = dens_guess
    result.ratio = ratio
    if (n_elements(ratio_err) ne 0L) then result.ratio_err = ratio_err
    
    dens_ions = ['S_II','O_II_DENS']
    temp_ions = ['O_II_TEMP','N_II','S_III','O_III'] 

    ndens = n_elements(dens_ions)
    ntemp = n_elements(temp_ions)

; loop on each ION/RATIO; note that I'm using INTERPOLATE in
; combination with FINDEX because when the ratio hits against the
; limits of the density and temperature grid (see
; WRITE_TEMDEN_LOOKUP_TABLE) the last physically meaningful
; temperature and/or density is returned

;   chi2 = exp(-[0.5*(ratio[0]-temden_table.oiii_ratio)^2.0/ratio_err[0]^2.0 + $
;     0.5*(ratio[1]-temden_table.sii_ratio)^2.0/ratio_err[1]^2.0])

;   chi2 = exp(-0.5*(ratio-temden_table.oiii_ratio[*,40])^2.0/ratio_err^2.0)
;   plot, temden_table.grid_temp, chi2, charsize=2

    for ii = 0L, nion-1L do begin
; loop on the known density diagnostics
       for idens = 0L, ndens-1L do begin
          if (strupcase(strtrim(ionlist[ii],2)) eq dens_ions[idens]) then begin
             result[ii].dens_ion = 1B
             case dens_ions[idens] of
                'S_II': ion_ratio = transpose(temden_table.sii_ratio)
                'O_II_DENS': ion_ratio = transpose(temden_table.oii_dens_ratio)
             endcase
             dens_model_ratio = interpolate(ion_ratio,findex(temden_table.grid_temp,temp_guess))
             dens = interpolate(temden_table.grid_dens,findex(dens_model_ratio,ratio[ii]))
             result[ii].dens = dens
             result[ii].model_ratio = interpolate(dens_model_ratio,findex(temden_table.grid_dens,dens))
             result[ii].ratio_frac_err = (result[ii].ratio-result[ii].model_ratio)/$
               result[ii].model_ratio
             if (dens ge max(temden_table.grid_dens)) or (dens le min(temden_table.grid_dens)) then $
               result[ii].dens_edge = 1B
          endif
       endfor
; loop on the known temperature diagnostics
       for itemp = 0L, ntemp-1L do begin
          if (strupcase(strtrim(ionlist[ii],2)) eq temp_ions[itemp]) then begin
             result[ii].temp_ion = 1B
             case temp_ions[itemp] of
                'O_II_TEMP': temp_model_ratio = interpolate(temden_table.oii_temp_ratio,$
                  findex(temden_table.grid_dens,dens_guess))
                'N_II': temp_model_ratio = interpolate(temden_table.nii_ratio,$
                  findex(temden_table.grid_dens,dens_guess))
                'S_III': temp_model_ratio = interpolate(temden_table.siii_ratio,$
                  findex(temden_table.grid_dens,dens_guess))
                'O_III': temp_model_ratio = interpolate(temden_table.oiii_ratio,$
                  findex(temden_table.grid_dens,dens_guess))
             endcase
             temp = interpolate(temden_table.grid_temp,findex(temp_model_ratio,ratio[ii]))
             result[ii].temp = temp
             result[ii].model_ratio = interpolate(temp_model_ratio,findex(temden_table.grid_temp,temp))
             result[ii].ratio_frac_err = (result[ii].ratio-result[ii].model_ratio)/$
               result[ii].model_ratio
             if (temp ge max(temden_table.grid_temp)) or (temp le min(temden_table.grid_temp)) then $
               result[ii].temp_edge = 1B
          endif
       endfor
    endfor

;; special case for most of what I'm doing right now (2007); iterate on
;; the temperature and density if [O III] and [S II] have been
;; computed; this is important at the ~100 K level
;
;    sii_indx = where(strlowcase(strtrim(result.ion,2)) eq 's_ii',nsii_indx)
;    oiii_indx = where(strlowcase(strtrim(result.ion,2)) eq 'o_iii',noiii_indx)
;
;    if (nsii_indx eq 1L) and (noiii_indx eq 1L) and (not keyword_set(noiter)) then begin
;   niter = 2L
;   resiter = result
;   for iter = 0L, niter do begin
;      dens_guess = resiter[sii_indx].dens
;      temp_guess = resiter[oiii_indx].temp
;      resiter = im_temden(['o_iii','s_ii'],resiter[[oiii_indx,sii_indx]].ratio,$
;        dens_guess=dens_guess,temp_guess=temp_guess,/noiter)
;   endfor
;   result[[oiii_indx,sii_indx]].dens = resiter.dens
;   result[[oiii_indx,sii_indx]].temp = resiter.temp
;    endif

; compute the uncertainties on the temperature and density using a
; Monte Carlo method

    if (nmonte gt 0L) and (n_elements(ratio_err) ne 0L) then begin

       for ii = 0L, nion-1L do begin
          
          ratio_monte = result[ii].ratio + randomn(seed,nmonte)*result[ii].ratio_err

          for imonte = 0L, nmonte-1L do begin
             if (result[ii].temp_ion eq 1B) and (result[ii].temp_edge eq 0B) then begin ; temperature
                resmonte1 = im_temden(ionlist[ii],ratio_monte[imonte],/noiter,$
                  temp_guess=result[ii].temp,dens_guess=result[ii].dens)
                if (imonte eq 0L) then resmonte = resmonte1 else $
                  resmonte = [resmonte,resmonte1]
             endif
             if (result[ii].dens_ion eq 1B) and (result[ii].dens_edge eq 0B) then begin ; density
                resmonte1 = im_temden(ionlist[ii],ratio_monte[imonte],/noiter,$
                  dens_guess=result[ii].dens,temp_guess=result[ii].temp)
                if (imonte eq 0L) then resmonte = resmonte1 else $
                  resmonte = [resmonte,resmonte1]
             endif
          endfor

; discard Monte Carlo realizations that result in large fractional
; differences between the measured and the model line-ratios as being
; unphysical 
          if result[ii].temp_ion and (n_elements(resmonte) ne 0L) then begin
             temp_good = where((abs(resmonte.ratio_frac_err) lt frac_tol),ngood)
             if (ngood gt 3L) then begin
                result[ii].temp_err = djsig(resmonte[temp_good].temp)
                result[ii].temp_quant = weighted_quantile(resmonte[temp_good].temp,quant=quant)
             endif
          endif
          if result[ii].dens_ion and (n_elements(resmonte) ne 0L) then begin
             dens_good = where((abs(resmonte.ratio_frac_err) lt frac_tol),ngood)
             if (ngood gt 3L) then begin
                result[ii].dens_err = djsig(resmonte[dens_good].dens)
                result[ii].dens_quant = weighted_quantile(resmonte[dens_good].dens,quant=quant)
             endif
          endif
       endfor 
    
    endif

return, result
end
