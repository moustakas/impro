;+
; NAME:
;   DEEP2_FILTERLIST()
; PURPOSE:
;   Return the DEEP2 filters in the extended photometric catalog of
;   Matthews et al. (2013).
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Jul 14, Siena
;-

function deep2_filterlist, field24=field24
    if keyword_set(field24) then begin
       filt = [['deep_'+['B','R','I']]+'.par',sdss_filterlist()]
    endif else begin
       filt = ['deep_'+['B','R','I']+'.par',cfhtls_filterlist()]
    endelse
return, filt
end
