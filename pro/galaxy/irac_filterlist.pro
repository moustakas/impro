;+
; NAME:
;   IRAC_FILTERLIST()
; PURPOSE:
;   Return the IRAC filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Nov 13, UCSD
;-

function irac_filterlist, warm_mission=warm_mission
    filt = 'spitzer_irac_'+['ch1','ch2','ch3','ch4']+'.par'
    if keyword_set(warm_mission) then filt = filt[0:1]
return, filt
end
