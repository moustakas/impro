;+
; NAME:
;   WISE_FILTERLIST()
; PURPOSE:
;   Return the WISE filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Apr 19, UCSD
;-

function wise_filterlist, longwavelength=longwavelength, $
  shortwavelength=shortwavelength
    filt = 'wise_w'+['1','2','3','4']+'.par'
    if keyword_set(longwavelength) then filt = filt[2:3]
    if keyword_set(shortwavelength) then filt = filt[0:1]
return, filt
end
