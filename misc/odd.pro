;+
; NAME:
;   ODD()
; PURPOSE:
;   Returns 0 if number if even, 1 if odd.
; MODIFICATION HISTORY:
;   8/1/95 MCL
;-

function odd,num
return,(long(num)/2 ne num/2.0)
end
