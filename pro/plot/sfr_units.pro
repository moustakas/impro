;+
; NAME:
;   SFR_UNITS()
; PURPOSE:
;   Return a clean title for star formation rate.
; MODIFICATION HISTORY:
;   J. Moustakas ???
;-

function sfr_units
return, textoidl('M'+sunsymbol()+' yr^{-1}')
end
