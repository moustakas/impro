;+
; NAME:
;   L
; PURPOSE:
;   Spawn the unix command 'ls -l | more'.
; MODIFICATION HISTORY:
;   J. Moustakas, 2000 Sep 14, U of A
;-

pro l
spawn, 'ls -lG | grep -v "~" | more', /sh
return
end
