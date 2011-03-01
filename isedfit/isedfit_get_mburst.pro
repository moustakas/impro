function isedfit_get_mburst, ised
; jm10feb01ucsd - given an isedfit structure, get the burst mass

    ngal = n_elements(ised)
    mburst = fltarr(ngal)
    for ii = 0L, ngal-1 do begin
       delvarx, age
       if (ised[ii].nburst ge 1) then begin
          sfh = isedfit_reconstruct_sfh(ised[ii],age=age)
stop          
       endif
    endfor
       
return, mburst
end
