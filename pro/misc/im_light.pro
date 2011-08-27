;+
; NAME:
;   IM_LIGHT()
; PURPOSE:
;   Return the speed of light in various units.
; KEYWORD PARAMETERS:
;   meters - 
;   cms - 
;   micron - 
;   angstrom - 
;   kms - (default)
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Jan 01, UCSD
;-

function im_light, meters=meters, cms=cms, micron=micron, angstrom=angstrom, kms=kms
    light = 2.99792458D5        ; [km/s]
    if keyword_set(meters) then light = light*1D3
    if keyword_set(cms) then light = light*1D5
    if keyword_set(micron) then light = light*1D9
    if keyword_set(angstrom) then light = light*1D13
return, light
end    
