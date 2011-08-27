;+
; NAME:
;   WISE_FILTERLIST()
; PURPOSE:
;   Return the WISE filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Apr 19, UCSD
;-

function wise_filterlist
    return, 'wise_w'+['1','2','3','4']+'.par'
end
