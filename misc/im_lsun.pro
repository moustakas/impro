function im_lsun, watts=watts
; jm10feb17ucsd - return the solar luminosity 
    lsun = 3.826D33 ; [erg/s]
    if keyword_set(watts) then lsun = 3.826D26 ; [W]
return, lsun
end    
