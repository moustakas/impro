;+
; NAME:
;   MIPS_FILTERLIST()
; PURPOSE:
;   Return the MIPS filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Feb 15, Siena
;-

function mips_filterlist
    filt = 'spitzer_mips_'+['24','70','160']+'.par'
return, filt
end
