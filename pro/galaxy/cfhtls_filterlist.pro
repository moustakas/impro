;+
; NAME:
;   CFHTLS_FILTERLIST()
; PURPOSE:
;   Return the CFHTLS filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Jun 24, Siena
;-

function cfhtls_filterlist
    return, 'capak_cfht_megaprime_sagem_'+['u','g','r','i','z']+'.par'
end
