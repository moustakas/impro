;+
; NAME:
;   GALEX_FILTERLIST()
; PURPOSE:
;   Return the GALEX filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Oct 30, UCSD
;-

function galex_filterlist
    return, 'galex_'+['FUV','NUV']+'.par'
end
