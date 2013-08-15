;+
; NAME:
;   GET_PWD()
; PURPOSE:
;   Return the present working directory.
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Aug 11, Siena
;-

function get_pwd
    cd, current=pwd
return, pwd+'/'
end
