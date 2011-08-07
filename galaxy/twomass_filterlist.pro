;+
; NAME:
;   TWOMASS_FILTERLIST()
; PURPOSE:
;   Return the 2MASS filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Aug 14, UCSD
;-

function twomass_filterlist
    return, 'twomass_'+['J','H','Ks']+'.par'
end
