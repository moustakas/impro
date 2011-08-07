;+
; NAME:
;   SDSS_FILTERLIST()
; PURPOSE:
;   Return the SDSS filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Aug 14, UCSD
;-

function sdss_filterlist
    return, 'sdss_'+['u0','g0','r0','i0','z0']+'.par'
end
