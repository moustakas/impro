;+
; NAME:
;   Z2STRING()
; PURPOSE:
;   Convert a metallicity value to a convenient string format.
; COMMENTS:
;   ISEDFIT support routine.
; MODIFICATION HISTORY:
;   J. Moustakas, ???
;-

function z2string, Z
return, 'Z'+string(Z,format='(G0.0)')
;return, 'Z'+string(Z,format='(F6.4)')
end
