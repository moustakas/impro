function abundance_catalogs_log12oh, data, ewbranch=ewbranch, branch=branch
; jm06apr11uofa - assign branches and compute final abundances

    niihacut = -1.05
    
; ---------------------------------------------------------------------------
; EW abundances    
; ---------------------------------------------------------------------------

    if (n_elements(ewbranch) ne 0L) then begin ; use assigned branches

       if (n_elements(ewbranch) ne n_elements(data)) then message, 'Incompatible dimensions.'
    
       good = where((ewbranch eq 'U') and (data.zstrong_ew_12oh_kk04_r23_upper gt -900) and (data.zstrong_ew_12oh_kk04_r23_lower gt -900),ngood)
       if (ngood ne 0L) then begin
          data[good].r23branch_kk04_ew        = 'U'
          data[good].zstrong_ew_12oh_kk04     = data[good].zstrong_ew_12oh_kk04_r23_upper
          data[good].zstrong_ew_12oh_kk04_err = data[good].zstrong_ew_12oh_kk04_r23_upper_err
       endif
   
       good = where((ewbranch eq 'T') and (data.zstrong_ew_12oh_kk04_r23_upper gt -900) and (data.zstrong_ew_12oh_kk04_r23_lower gt -900),ngood)
       if (ngood ne 0L) then begin
          data[good].r23branch_kk04_ew = 'T'
          for i = 0L, ngood-1L do begin
             oh = [data[good[i]].zstrong_ew_12oh_kk04_r23_upper,data[good[i]].zstrong_ew_12oh_kk04_r23_lower]
             oh_err = [data[good[i]].zstrong_ew_12oh_kk04_r23_upper_err,data[good[i]].zstrong_ew_12oh_kk04_r23_lower_err]
             data[good[i]].zstrong_ew_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
             data[good[i]].zstrong_ew_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
          endfor
       endif
   
       good = where((ewbranch eq 'L') and (data.zstrong_ew_12oh_kk04_r23_upper gt -900) and (data.zstrong_ew_12oh_kk04_r23_lower gt -900),ngood)
       if (ngood ne 0L) then begin
          data[good].r23branch_kk04_ew        = 'L'
          data[good].zstrong_ew_12oh_kk04     = data[good].zstrong_ew_12oh_kk04_r23_lower
          data[good].zstrong_ew_12oh_kk04_err = data[good].zstrong_ew_12oh_kk04_r23_lower_err
       endif

    endif else begin            ; determine branches
    
       good = where((data.zstrong_ew_12oh_kk04_r23_upper gt -900) and (data.zstrong_ew_12oh_kk04_r23_lower gt -900),ngood)
       if (ngood ne 0L) then begin
          data[good].r23branch_kk04_ew        = 'U'
          data[good].zstrong_ew_12oh_kk04     = data[good].zstrong_ew_12oh_kk04_r23_upper
          data[good].zstrong_ew_12oh_kk04_err = data[good].zstrong_ew_12oh_kk04_r23_upper_err
       endif

;      good = where((data.zstrong_ew_12oh_kk04_r23_upper gt -900) and (data.zstrong_ew_12oh_kk04_r23_lower gt -900),ngood)
;      if (ngood ne 0L) then begin
;         turn = where((data[good].zstrong_ew_12oh_kk04_r23_upper-data[good].zstrong_ew_12oh_kk04_r23_upper_err) lt $
;                      (data[good].zstrong_ew_12oh_kk04_r23_lower+data[good].zstrong_ew_12oh_kk04_r23_lower_err),nturn)
;         if (nturn ne 0L) then begin
;            data[good[turn]].r23branch_kk04_ew = 'T'
;            for i = 0L, nturn-1L do begin
;               oh = [data[good[turn[i]]].zstrong_ew_12oh_kk04_r23_upper,data[good[turn[i]]].zstrong_ew_12oh_kk04_r23_lower]
;               oh_err = [data[good[turn[i]]].zstrong_ew_12oh_kk04_r23_upper_err,data[good[turn[i]]].zstrong_ew_12oh_kk04_r23_lower_err]
;               data[good[turn[i]]].zstrong_ew_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
;               data[good[turn[i]]].zstrong_ew_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
;            endfor
;         endif
;      endif
    
       good = where((data.zstrong_ew_niiha gt -900.0),ngood)
;      good = where((data.h_alpha_ew[1] gt 0.0) and (data.nii_6584_ew[1] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          lo = where(alog10(data[good].nii_6584_ew[0]/data[good].h_alpha_ew[0]) lt niihacut,nlo)
          if (nlo ne 0L) then begin
             data[good[lo]].r23branch_kk04_ew        = 'L'
             data[good[lo]].zstrong_ew_12oh_kk04     = data[good[lo]].zstrong_ew_12oh_kk04_r23_lower
             data[good[lo]].zstrong_ew_12oh_kk04_err = data[good[lo]].zstrong_ew_12oh_kk04_r23_lower_err
          endif
       endif
       
    endelse

; reject    
    
    reject = where((data.zstrong_ew_12oh_kk04_r23_upper gt -900) and (data.zstrong_ew_12oh_kk04_r23_lower gt -900) and $
      ((data.zstrong_ew_12oh_kk04_r23_upper lt data.zstrong_ew_12oh_kk04_r23_lower) and $
       ((data.zstrong_ew_12oh_kk04_r23_upper+data.zstrong_ew_12oh_kk04_r23_upper_err) lt $ ; note!
        (data.zstrong_ew_12oh_kk04_r23_lower-data.zstrong_ew_12oh_kk04_r23_lower_err))),nreject)
    if (nreject ne 0L) then begin
       data[reject].r23branch_kk04_ew        = 'Rejected'
       data[reject].zstrong_ew_12oh_kk04     = -999.0
       data[reject].zstrong_ew_12oh_kk04_err = -999.0
    endif

; if the lower branch metallicity is larger than the upper branch
; metallicity, but the object was not rejected from the above
; algorithm (because the measurements were within one-sigma of one
; another; e.g., GDDS/SA02-0756, or CFRS/03.0599), then average them!
    
    average = where((data.zstrong_ew_12oh_kk04_r23_upper gt -900) and (data.zstrong_ew_12oh_kk04_r23_lower gt -900) and $
      ((data.zstrong_ew_12oh_kk04_r23_upper lt data.zstrong_ew_12oh_kk04_r23_lower) and $
       ((data.zstrong_ew_12oh_kk04_r23_upper+data.zstrong_ew_12oh_kk04_r23_upper_err) gt $
        (data.zstrong_ew_12oh_kk04_r23_lower-data.zstrong_ew_12oh_kk04_r23_lower_err))),naverage)
    if (naverage ne 0L) then begin
;      niceprint, data[average].zstrong_ew_12oh_kk04_r23_upper, data[average].zstrong_ew_12oh_kk04_r23_upper_err, $
;        data[average].zstrong_ew_12oh_kk04_r23_lower, data[average].zstrong_ew_12oh_kk04_r23_lower_err
       data[average].r23branch_kk04_ew        = 'T'
       for i = 0L, naverage-1L do begin
          oh = [data[average[i]].zstrong_ew_12oh_kk04_r23_upper,data[average[i]].zstrong_ew_12oh_kk04_r23_lower]
          oh_err = [data[average[i]].zstrong_ew_12oh_kk04_r23_upper_err,data[average[i]].zstrong_ew_12oh_kk04_r23_lower_err]
          data[average[i]].zstrong_ew_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          data[average[i]].zstrong_ew_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

; ---------------------------------------------------------------------------
; flux abundances    
; ---------------------------------------------------------------------------
    
    if (n_elements(ewbranch) ne 0L) then begin ; use assigned branches

       if (n_elements(ewbranch) ne n_elements(data)) then message, 'Incompatible dimensions.'

       good = where((branch eq 'U') and (data.zstrong_12oh_kk04_r23_upper gt -900) and (data.zstrong_12oh_kk04_r23_lower gt -900),ngood)
       if (ngood ne 0L) then begin
          data[good].r23branch_kk04        = 'U'
          data[good].zstrong_12oh_kk04     = data[good].zstrong_12oh_kk04_r23_upper
          data[good].zstrong_12oh_kk04_err = data[good].zstrong_12oh_kk04_r23_upper_err
       endif
 
       good = where((branch eq 'T') and (data.zstrong_12oh_kk04_r23_upper gt -900) and (data.zstrong_12oh_kk04_r23_lower gt -900),ngood)
       if (ngood ne 0L) then begin
          data[good].r23branch_kk04 = 'T'
          for i = 0L, ngood-1L do begin
             oh = [data[good[i]].zstrong_12oh_kk04_r23_upper,data[good[i]].zstrong_12oh_kk04_r23_lower]
             oh_err = [data[good[i]].zstrong_12oh_kk04_r23_upper_err,data[good[i]].zstrong_12oh_kk04_r23_lower_err]
             data[good[i]].zstrong_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
             data[good[i]].zstrong_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
          endfor
       endif
   
       good = where((branch eq 'L') and (data.zstrong_12oh_kk04_r23_upper gt -900) and (data.zstrong_12oh_kk04_r23_lower gt -900),ngood)
       if (ngood ne 0L) then begin
          data[good].r23branch_kk04        = 'L'
          data[good].zstrong_12oh_kk04     = data[good].zstrong_12oh_kk04_r23_lower
          data[good].zstrong_12oh_kk04_err = data[good].zstrong_12oh_kk04_r23_lower_err
       endif
 
    endif else begin ; determine branches
    
       good = where((data.zstrong_12oh_kk04_r23_upper gt -900) and (data.zstrong_12oh_kk04_r23_lower gt -900),ngood)
       if (ngood ne 0L) then begin
          data[good].r23branch_kk04        = 'U'
          data[good].zstrong_12oh_kk04     = data[good].zstrong_12oh_kk04_r23_upper
          data[good].zstrong_12oh_kk04_err = data[good].zstrong_12oh_kk04_r23_upper_err
       endif
 
;      good = where((data.zstrong_12oh_kk04_r23_upper gt -900) and (data.zstrong_12oh_kk04_r23_lower gt -900),ngood)
;      if (ngood ne 0L) then begin
;         turn = where((data[good].zstrong_12oh_kk04_r23_upper-data[good].zstrong_12oh_kk04_r23_upper_err) lt $
;                      (data[good].zstrong_12oh_kk04_r23_lower+data[good].zstrong_12oh_kk04_r23_lower_err),nturn)
;         if (nturn ne 0L) then begin
;            data[good[turn]].r23branch_kk04 = 'T'
;            for i = 0L, nturn-1L do begin
;               oh = [data[good[turn[i]]].zstrong_12oh_kk04_r23_upper,data[good[turn[i]]].zstrong_12oh_kk04_r23_lower]
;               oh_err = [data[good[turn[i]]].zstrong_12oh_kk04_r23_upper_err,data[good[turn[i]]].zstrong_12oh_kk04_r23_lower_err]
;               data[good[turn[i]]].zstrong_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
;               data[good[turn[i]]].zstrong_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
;            endfor
;         endif
;      endif
 
       good = where((data.zstrong_niiha gt -900.0),ngood)
;      good = where((data.h_alpha[1] gt 0.0) and (data.nii_6584[1] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          lo = where(alog10(data[good].nii_6584[0]/data[good].h_alpha[0]) lt niihacut,nlo)
          if (nlo ne 0L) then begin
             data[good[lo]].r23branch_kk04        = 'L'
             data[good[lo]].zstrong_12oh_kk04     = data[good[lo]].zstrong_12oh_kk04_r23_lower
             data[good[lo]].zstrong_12oh_kk04_err = data[good[lo]].zstrong_12oh_kk04_r23_lower_err
          endif

       endif

    endelse
       
; reject    
    
    reject = where((data.zstrong_12oh_kk04_r23_upper gt -900) and (data.zstrong_12oh_kk04_r23_lower gt -900) and $
      ((data.zstrong_12oh_kk04_r23_upper lt data.zstrong_12oh_kk04_r23_lower) and $
       ((data.zstrong_12oh_kk04_r23_upper+data.zstrong_12oh_kk04_r23_upper_err) lt $ ; note!
        (data.zstrong_12oh_kk04_r23_lower-data.zstrong_12oh_kk04_r23_lower_err))),nreject)
    if (nreject ne 0L) then begin
       data[reject].r23branch_kk04        = 'Rejected'
       data[reject].zstrong_12oh_kk04     = -999.0
       data[reject].zstrong_12oh_kk04_err = -999.0
    endif

; if the lower branch metallicity is larger than the upper branch
; metallicity, but the object was not rejected from the above
; algorithm (because the measurements were within one-sigma of one
; another; e.g., GDDS/SA02-0756, or CFRS/03.0599), then average them!
    
    average = where((data.zstrong_12oh_kk04_r23_upper gt -900) and (data.zstrong_12oh_kk04_r23_lower gt -900) and $
      ((data.zstrong_12oh_kk04_r23_upper lt data.zstrong_12oh_kk04_r23_lower) and $
       ((data.zstrong_12oh_kk04_r23_upper+data.zstrong_12oh_kk04_r23_upper_err) gt $
        (data.zstrong_12oh_kk04_r23_lower-data.zstrong_12oh_kk04_r23_lower_err))),naverage)
    if (naverage ne 0L) then begin
;      niceprint, data[average].zstrong_12oh_kk04_r23_upper, data[average].zstrong_12oh_kk04_r23_upper_err, $
;        data[average].zstrong_12oh_kk04_r23_lower, data[average].zstrong_12oh_kk04_r23_lower_err
       data[average].r23branch_kk04        = 'T'
       for i = 0L, naverage-1L do begin
          oh = [data[average[i]].zstrong_12oh_kk04_r23_upper,data[average[i]].zstrong_12oh_kk04_r23_lower]
          oh_err = [data[average[i]].zstrong_12oh_kk04_r23_upper_err,data[average[i]].zstrong_12oh_kk04_r23_lower_err]
          data[average[i]].zstrong_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          data[average[i]].zstrong_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

return, data
end
