;+
; NAME:
;   CWD()
; PURPOSE:
;   Return the current working directory.
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 Jan 21, U of A
;-

function cwd
    spawn, 'pwd', d, /sh
    d = d[0]+'/'
return, d
end
