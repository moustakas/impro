;+
; NAME:
;   DECAM_FILTERLIST()
; PURPOSE:
;   Return the DECam filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Jun 8, Siena
;-

function decam_filterlist
    return, 'decam_'+['u','g','r','i','z','Y']+'.par'
end
