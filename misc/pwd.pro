;+
; NAME:
;   PWD()
; PURPOSE:
;   Return the present working directory.
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 Jan 21, U of A
;-

pro pwd
    spawn, 'pwd', /sh
return
end
