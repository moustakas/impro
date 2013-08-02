;+
; NAME:
;   DEEP2_FILTERLIST()
; PURPOSE:
;   Return the DEEP2 filters in the extended photometric catalog of
;   Matthews et al. (2013).
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Jul 14, Siena
;-

function deep2_filterlist
    return, [['deep_'+['B','R','I']]+'.par',sdss_filterlist()]
end
