function compute_inclination, diam, twomass=twomass
; jm06feb16uofa - compute the inclination angle; the input should be
;                 the output from NED_WEBGET_DIAMETERS, or similar
;                 structure 

    ndiam = n_elements(diam)
    if (ndiam eq 0L) then begin
       print, 'Syntax - incl = compute_inclination(diam,/twomass)'
       return, -1L
    endif

    incl = fltarr(ndiam)-999.0
    
    if keyword_set(twomass) then begin

       good = where((diam.twomass_k20_major_axis gt -900.0) and (diam.twomass_k20_minor_axis gt -900.0),ngood)
       if (ngood ne 0L) then begin

          d25_maj = float(diam[good].twomass_k20_major_axis)
          d25_min = float(diam[good].twomass_k20_minor_axis)

          ratio = d25_min/d25_maj ; "b" divided by "a"
          quantity = sqrt((1.0-ratio^2)/0.96)

          incl[good] = asin(quantity<1)*!radeg

       endif

    endif else begin

       good = where((diam.rc3_major_axis gt -900.0) and (diam.rc3_minor_axis gt -900.0),ngood)
       if (ngood ne 0L) then begin

          d25_maj = float(diam[good].rc3_major_axis)
          d25_min = float(diam[good].rc3_minor_axis)

          ratio = d25_min/d25_maj ; "b" divided by "a"
          quantity = sqrt((1.0-ratio^2)/0.96)

          incl[good] = asin(quantity<1)*!radeg
       
       endif

    endelse
       
return, incl
end
