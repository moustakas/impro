;+
; NAME:
;   IM_LSUN()
; PURPOSE:
;   Return the solar luminosity in erg/s.
; KEYWORD PARAMETERS:
;   watts - 
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 17, UCSD
;-

function im_lsun, watts=watts
    lsun = 3.826D33 ; [erg/s]
    if keyword_set(watts) then lsun = 3.826D26 ; [W]
return, lsun
end    
