;+
; NAME:
;   FNU_UNITS()
; PURPOSE:
;   Return a clean title for flux density.
; MODIFICATION HISTORY:
;   J. Moustakas ???
;-

function fnu_units
return, textoidl('erg s^{-1} cm^{-2} Hz^{-1}')
end
