;+
; NAME:
;   GET_MJDSTR()
; PURPOSE:
;   Return the current modified Julian date (MJD).
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Aug 11, NYU
;-

function get_mjdstr
    get_juldate, jd
return, string(long(jd-2400000L),format='(I5)')
end
