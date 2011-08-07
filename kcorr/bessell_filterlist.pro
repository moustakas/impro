;+
; NAME:
;   BESSELL_FILTERLIST()
; PURPOSE:
;   Return the Bessell filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Aug 14, UCSD
;-

function bessell_filterlist
    return, 'bessell_'+['U','B','V','R','I']+'.par'
end
