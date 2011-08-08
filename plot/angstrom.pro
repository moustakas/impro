;+
; NAME:
;   ANGSTROM()
; PURPOSE:
;   Render the Angstrom symbol on a plot.
; MODIFICATION HISTORY:
;   J. Moustakas ???
;-

function angstrom
;return, '!6!sA!r!u!9 %!6!n'
;return, string(byte(197))
;return, '!3' + STRING(197B) + '!X'
return, string('305'OB)
end
