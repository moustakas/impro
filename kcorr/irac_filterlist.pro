;+
; NAME:
;   IRAC_FILTERLIST()
; PURPOSE:
;   Return the IRAC filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Nov 13, UCSD
;-

function irac_filterlist
    return, 'spitzer_irac_'+['ch1','ch2','ch3','ch4']+'.par'
end
