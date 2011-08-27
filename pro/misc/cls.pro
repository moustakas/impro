;+
; NAME:
;   CLS
; PURPOSE:
;   Spawn the unix 'clear' command.
; MODIFICATION HISTORY:
;   J. Moustakas, 2000 Sep 14, U of A
;-

pro cls
spawn, 'clear', /sh
return
end
