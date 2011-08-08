;+
; NAME:
;   FLAM_UNITS()
; PURPOSE:
;   Return a clean title for flux density.
; MODIFICATION HISTORY:
;   J. Moustakas ???
;-

function flam_units
return, textoidl('erg s^{-1} cm^{-2} \AA^{-1}')
end
