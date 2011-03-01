function im_light, meters=meters, cm=cm, micron=micron, angstrom=angstrom, kms=kms
; jm10jan01ucsd - return the speed of light
    light = 2.99792458D5        ; [km/s]
    if keyword_set(meters) then light = light*1D3
    if keyword_set(cm) then light = light*1D5
    if keyword_set(micron) then light = light*1D9
    if keyword_set(angstrom) then light = light*1D13
return, light
end    
