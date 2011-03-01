;+
; NAME:
;       IDL_FIVEL()
;
; PURPOSE:
;       Compute nebular physical conditions or emission-line
;       emissivities.
;
; INPUTS:
;       ictg - scalar long integer: 1 - compute the density and
;              temperature given the appropriate line ratio; 2 -
;              compute emission-line emissivities given the density
;              and temperature
;       ion  - atomic ion number (see COMMENTS)
;
; OPTIONAL INPUTS:
;       lineratio       - emission-line ratio (see EXAMPLES)
;       err_lineratio   - error in LINERATIO
;       temperature     - if ICTG=1 then this variable is the initial 
;                         temperature guess (default 10,000); if
;                         ICTG=2 then this variable is the nebular
;                         temperature 
;       min_temperature - minimum temperature
;       max_temperature - maximum temperature
;       density         - if ICTG=1 then this variable is the initial
;                         density guess (default 100); if ICTG=2 then
;                         this variable is the nebular density  
;
; KEYWORD PARAMETERS:
;       montecarlo - Monte Carlo the bad boy
;       silent     - suppress messages to STDOUT
;
; OUTPUTS:
;       result - 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Jun 30, U of A - written
;       jm05jan20uofa - fixed an error in how the upper and lower
;                       standard deviations were being calculated 
;
; Copyright (C) 2004-2005, John Moustakas
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

function compute_emissivity, t, jl, density, lambda=lambda
; compute the emission-line wavelengths and emissivities; if the
; emissivity is less than 1D-35 then set it to zero

    lambda = dblarr(10)
    emissivity = dblarr(10)
    
    indx = 0L
    for i = 1L, 4L do begin
       for j = 0L, 3L do begin
          if (j lt i) then begin
             lambda[indx] = 1.0D/(t[i]-t[j])
             emissivity[indx] = jl[i,j]/density
             if (emissivity[indx] lt 1D-35) then emissivity[indx] = 0.0D
             indx = indx+1L
          endif
       endfor
    endfor

;   vactoair, lambda ; convert to wavelengths in air

return, emissivity
end    

function idl_fivel, ictg, ion, lineratio=lineratio, err_lineratio=err_lineratio, $
  temperature=temperature, min_temperature=min_temperature, max_temperature=max_temperature, $
  density=density, montecarlo=montecarlo, silent=silent, nmonte=nmonte

    nictg = n_elements(ictg)
    nion = n_elements(ion)
    
    if (nictg ne 1L) or (nion ne 1L) then begin
       doc_library, 'idl_fivel'
       return, -1L
    endif

    if (n_elements(nmonte) eq 0L) then nmonte = 500L

    nlineratio = n_elements(lineratio)
    nerr_lineratio = n_elements(err_lineratio)
    ntemperature = n_elements(temperature)
    ndensity = n_elements(density)

    ictg = long(ictg)
    ion = long(ion)
    
    if (ictg ne 1L) and (ictg ne 2L) then begin
       splog, 'ICTG must be either 1 or 2.'
       return, -1L
    endif

    if (ictg eq 1L) then begin

;      if (n_elements(silent) eq 0L) then splog, 'Computing nebular physical conditions.'

       if (nlineratio eq 0L) then begin
          splog, 'LINERATIO must be specified.'
          return, -1L
       endif else lineratio = double(lineratio)
       
       if (nerr_lineratio eq 0L) then err_lineratio = lineratio*0.0 else begin
          if (nerr_lineratio ne nlineratio) then begin
             splog, 'LINERATIO and ERR_LINERATIO must have the same number of elements.'
             return, -1L
          endif else err_lineratio = double(err_lineratio)
       endelse
       
       if (ntemperature eq 0L) then temperature = lineratio*0.0+1D4 else begin
          if (ntemperature ne nlineratio) then begin
             splog, 'LINERATIO and TEMPERATURE must have the same number of elements.'
             return, -1L
          endif else temperature = double(temperature)
       endelse
       
       if (ndensity eq 0L) then density = lineratio*0.0+1D2 else begin
          if (ndensity ne nlineratio) then begin
             splog, 'LINERATIO and DENSITY must have the same number of elements.'
             return, -1L
          endif else density = double(density)
       endelse

       nerr_lineratio = n_elements(err_lineratio)
       ntemperature = n_elements(temperature)
       ndensity = n_elements(density)
       
       case ion of
          0L: begin
             ionum = 6021L
             if (n_elements(silent) eq 0L) then splog, 'Computing the density from C III] - I(1909)/(1907).'
          end
          1L: begin
             ionum = 7001L
             if (n_elements(silent) eq 0L) then splog, 'Computing the density from [N I] - I(5200)/(I5198).'
          end
          2L:  begin
             ionum = 7012L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [N II] - I(6548+6583)/I(5755).'
          end
          3L:  begin
             ionum = 8002L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [O I] - ???'
          end
          4L:  begin
             ionum = 8011L
             if (n_elements(silent) eq 0L) then splog, 'Computing the density from [O II] - I(3726)/I(3729).'
          end
          5L:  begin
             ionum = 8012L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [O II] - I(3727)/I(7324).'
          end
          6L:  begin
             ionum = 8022L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [O III] - I(4959+5007)/I(4363).'
          end
          7L:  begin
             ionum = 10022L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [Ne III] - I(3869+3968)/I(3343).'
          end
          8L:  begin
             ionum = 10031L
             if (n_elements(silent) eq 0L) then splog, 'Computing the density from [Ne IV] - I(2423)/I(2425).'
          end
          9L:  begin
             ionum = 10042L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [Ne V] - I(3426+3346)/I(2975).'
          end
          10L: begin
             ionum = 14021L
             if (n_elements(silent) eq 0L) then splog, 'Computing the density from [Si III] - I(1883)/I(1893).'
          end
          11L: begin
             ionum = 16011L
             if (n_elements(silent) eq 0L) then splog, 'Computing the density from [S II] - I(6717)/I(6731).'
          end
          12L: begin
             ionum = 16012L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [S II] - I(6717+6731)/I(4069+4076).'
          end
          13L: begin
             ionum = 16022L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [S III] - I(9069+9532)/I(6312).'
          end
          14L: begin
             ionum = 17021L
             if (n_elements(silent) eq 0L) then splog, 'Computing the density from [Cl III] - I(5518)/I(5538).'
          end
          15L: begin
             ionum = 17032L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [Cl IV] - I(7530+8045)/I(5323).'
          end
          16L: begin
             ionum = 18022L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [Ar III] - I(7136+7751)/I(5192).'
          end
          17L: begin
             ionum = 18031L
             if (n_elements(silent) eq 0L) then splog, 'Computing the density from [Ar IV] - I(4711)/I(4740).'
          end
          18L: begin
             ionum = 18032L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [Ar IV] - I(4711+4740)/I(2854+2868).'
          end
          19L: begin
             ionum = 18042L
             if (n_elements(silent) eq 0L) then splog, 'Computing the temperature from [Ar V] - I(6435+7006)/I(4626).'
          end
          else: begin
             splog, 'Unrecognized ion number '+string(ion,format='(I0)')+'.'
             return, -1L
          end
       endcase

    endif

    if (ictg eq 2L) then begin

;      if (n_elements(silent) eq 0L) then splog, 'Computing emissivities.'

       if (ntemperature eq 0L) or (ndensity eq 0L) then begin
          splog, 'TEMPERATURE and DENSITY must be specified.'
          return, -1L
       endif

       if (ntemperature ne ndensity) then begin
          splog, 'TEMPERATURE and DENSITY must have the same number of elements.'
          return, -1L
       endif

       lineratio = temperature*0.0-999.0D     ; internal variable
       err_lineratio = temperature*0.0-999.0D ; internal variable
       
       case ion of
          0L:  begin
             ionum = 6010L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [C II].'
          end
          1L:  begin
             ionum = 6020L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for C III].'
          end
          2L:  begin
             ionum = 7000L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [N I].'
          end
          3L:  begin
             ionum = 7010L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [N II].'
          end
          4L:  begin
             ionum = 7020L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [N III].'
          end
          5L:  begin
             ionum = 8000L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [O I].'
          end
          6L:  begin
             ionum = 8010L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [O II].'
          end
          7L:  begin
             ionum = 8020L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [O III].'
          end
          8L:  begin
             ionum = 10020L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [Ne III].'
          end
          9L:  begin
             ionum = 10030L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [Ne IV].'
          end
          10L:  begin
             ionum = 10040L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [Ne V].'
          end
          11L:  begin
             ionum = 14020L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for Si III].'
          end
          12L:  begin
             ionum = 16010L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [S II].'
          end
          13L:  begin
             ionum = 16020L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [S III].'
          end
          14L:  begin
             ionum = 16030L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [S IV].'
          end
          15L:  begin
             ionum = 17010L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [Cl II].'
          end
          16L:  begin
             ionum = 17020L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [Cl III].'
          end
          17L:  begin
             ionum = 17030L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [Cl IV].'
          end
          18L:  begin
             ionum = 18020L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [Ar III].'
          end
          19L:  begin
             ionum = 18030L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [Ar IV].'
          end
          20L:  begin
             ionum = 18040L
             if (n_elements(silent) eq 0L) then splog, 'Computing emissivities for [Ar V].'
          end
          else: begin
             splog, 'Unrecognized ion number '+string(ion,format='(I0)')+'.'
             return, -1L
          end
       endcase

    endif

; for multiple elements call this routine recursively

    if (ntemperature gt 1L) then begin

       for i = 0L, ntemperature-1L do begin

          result1 = idl_fivel(ictg,ion,lineratio=lineratio[i],err_lineratio=err_lineratio[i],$
            temperature=temperature[i],density=density[i],montecarlo=montecarlo,/silent)
          if (i eq 0L) then result = result1 else result = struct_append(result,result1)

       endfor
       return, reform(result)

    endif

; initialize variables for FIVEL.F and the output structure

    a = dblarr(5,5)
    c = dblarr(5,5)
    q = dblarr(5,5)
    jl = dblarr(5,5)
    e = dblarr(5,5)
    t = dblarr(5)
    pop = dblarr(5)
    ncrit = dblarr(5)
    weight = lonarr(5)
    temden_error = 0L

    result = {$
      ion:                ion,   $
      ionum:              ionum, $
      lineratio:          0.0D,  $
      lineratio_err:      0.0D,  $
      temperature:        0.0D,  $
      temperature_lower:  0.0D,  $
      temperature_upper:  0.0D,  $
      density:            0.0D,  $
      density_lower:      0.0D,  $
      density_upper:      0.0D,  $
      a:                  a,     $
      c:                  c,     $
      q:                  q,     $
      jl:                 jl,    $
      e:                  e,     $
      t:                  t,     $
      pop:                pop,   $
      ncrit:              ncrit, $
      weight:             weight,$
      temden_error:         0L,  $
      logx:               0.0D,  $
      jhb:                0.0D,  $
      jhb_err:            0.0D,  $
      lambda:         dblarr(10),$
      emissivity:     dblarr(10),$
      emissivity_err: dblarr(10)}

    srcdir = filepath('',root_dir=getenv('IMPRO_DIR'),subdirectory='nebular')

    flag = call_external(srcdir+'idl_fivel.so', 'idl_fivel_wrapper_', $
      long(ictg),long(ionum),double(temperature),double(density), $
      double(lineratio),double(a),double(c),double(q),double(jl), $
      double(e),double(t),double(pop),double(ncrit),long(weight), $
      long(temden_error))
    if (temden_error ne 0L) and (n_elements(silent) eq 0L) then $
      splog, 'Ratio predicts unreasonable conditions - restoring initial values.'

; store the results

    result.temperature = temperature
    result.density = density
    result.lineratio = lineratio
    result.lineratio_err = err_lineratio
    result.a = a
    result.c = c
    result.q = q
    result.jl = jl
    result.e = e
    result.t = t
    result.pop = pop
    result.ncrit = ncrit
    result.weight = weight
    result.temden_error = temden_error
    result.logx = alog10(density/sqrt(temperature))
    result.jhb = 1.387D-25 / (temperature/1D4)^0.983 / 10D^(0.0424D4/temperature)

; compute the emission-line wavelengths and emissivities 

    result.emissivity = compute_emissivity(result.t,result.jl,$
      result.density,lambda=lambda)
    result.lambda = lambda

; to get the median temperature and a proper error  Monte Carlo the
; line ratio NMONTE times
    
    if (n_elements(montecarlo) ne 0L) then begin

       if (n_elements(min_temperature) eq 0L) then min_temperature = temperature
       if (n_elements(max_temperature) eq 0L) then max_temperature = temperature
       temperature_range = (max_temperature - min_temperature)/2.0
       
       temperature_array = dblarr(nmonte)
       density_array = dblarr(nmonte)
       lineratio_array = dblarr(nmonte)
       emissivity_array = dblarr(10,nmonte)
       jhb_array = dblarr(nmonte)
       temden_array = lonarr(nmonte)
       
       for imonte = 0L, nmonte-1L do begin

; if ICTG=1 then Gaussian perturb the line ratio to compute the error
; in the temperature and density; if ICTG=2 then choose a flat prior
; in the temperature to compute the error in the emissivity
          
          if (ictg eq 1L) then begin
          
             lineratio_monte = lineratio + (randomn(seed,1))[0]*err_lineratio
             temperature_monte = temperature

          endif

          if (ictg eq 2L) then begin
          
             lineratio_monte = lineratio
;            temperature_monte = temperature + (randomu(seed,1)-0.5)[0]*temperature_range
             temperature_monte = temperature + (randomn(seed,1))[0]*temperature_range

          endif

          density_monte = density

          flag = call_external(srcdir+'idl_fivel.so', 'idl_fivel_wrapper_', $
            long(ictg),long(ionum),double(temperature_monte),double(density_monte), $
            double(lineratio_monte),double(a),double(c),double(q),double(jl), $
            double(e),double(t),double(pop),double(ncrit),long(weight), $
            long(temden_error))

          jhb = 1.387D-25 / (temperature_monte/1D4)^0.983 / 10D^(0.0424D4/temperature_monte)
          
          temperature_array[imonte] = temperature_monte
          density_array[imonte] = density_monte
          lineratio_array[imonte] = lineratio_monte
          emissivity_array[*,imonte] = compute_emissivity(t,jl,density_monte)
          jhb_array[imonte] = jhb

          temden_array[imonte] = temden_error

       endfor
       
       onesigma = errorf(1/sqrt(2.0))-0.5

       good = where(temden_array eq 0L,ngood) ; physical conditions predicted
       if (ngood gt 2L) then begin

          temperature_array = temperature_array[good]
          density_array = density_array[good]
          lineratio_array = lineratio_array[good]
          emissivity_array = emissivity_array[*,good]
          jhb_array = jhb_array[good]
          
; temperature

          stats = im_stats(temperature_array)
          result.temperature = stats.median
          result.temperature_lower = stats.median-stats.sig68lo
          result.temperature_upper = stats.sig68up-stats.median
          
; density
          
          stats = im_stats(density_array)
          result.density = stats.median
          result.density_lower = stats.median-stats.sig68lo
          result.density_upper = stats.sig68up-stats.median
          
; emissivity

          result.jhb_err = stddev(jhb_array)
          for jindx = 0L, 9L do result.emissivity_err[jindx] = stddev(emissivity_array[jindx,*])
          
       endif
          
    endif
    
return, result
end
